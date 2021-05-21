// Scorer for ClusteredDNADamage
//
//**************************************************************************************************
// Author: Logan Montgomery
// Based on code written by:
//      - J Schuemann et al. (2019). DOI:10.1667/RR15226.1
//		- A McNamara et al. (2018). DOI:10.1088/1361-6560/aad8eb
//      - E Delage et al. (2015). DOI:10.1016/j.cpc.2015.02.026
//
// This class is a custom Topas ntuple scorer that records clustered DNA damage induced in a DNA
// model (VoxelizedNuclearDNA).
//**************************************************************************************************

#include "ScoreClusteredDNADamage.hh"
#include "TsTrackInformation.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <algorithm>

#include <map>
#include "G4RunManager.hh"

#include "G4Molecule.hh"
#include "G4MoleculeTable.hh"
#include "G4MolecularConfiguration.hh"

//--------------------------------------------------------------------------------------------------
// Struct used to hold parameters of interest for a single cluster of DNA damage.
// Used in RecordClusteredDNADamage().
//--------------------------------------------------------------------------------------------------
struct DamageCluster {
	DamageCluster() : numSSB(0), numBD(0), numDSB(0), start(0), end(0), size(0){}

	G4int numSSB; // # of SSB in cluster
	G4int numBD; // # of base damages in cluster
	G4int numDSB; // # of DSB in cluster

	G4int start; // starting bp index of cluster
	G4int end; // ending bp index of cluster
	G4int size; // bp range of cluster (end-start+1)
};


//--------------------------------------------------------------------------------------------------
// Constructor. Initialize member variables using a variety of methods. Specify which data is output
// to the main data output file.
//--------------------------------------------------------------------------------------------------
ScoreClusteredDNADamage::ScoreClusteredDNADamage(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
											   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
:TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)

{
	//----------------------------------------------------------------------------------------------
	// Initialize member variables
	//----------------------------------------------------------------------------------------------
	ResolveParams(); // initialize some member variables using Topas parameter file

	// Damage yields
	fTotalSSB = 0;
	fTotalBD = 0;
	fTotalDSB = 0;
	fTotalComplexDSB = 0;
	fTotalNonDSBCluster = 0;

	fDoubleCountsDI = 0;
	fDoubleCountsII = 0;
	fNumProcessHitsCalls = 0;

	// Run tracking
	fNumEvents = 0;
	fTotalEdep = 0.;
	fFiberID = 0;
	fVoxelID = 0;

	// Variables used when generating output files
	fDelimiter = ",";
	fOutFileExtension = ".csv";
	fOutHeaderExtension = ".header";

	// Erase contents of existing output files. Must be done here, at start of run, in case doing
	// event-by-event scoring (i.e. need to write to same file many times).
	ClearOutputFiles();

	// Total cubic volume of the DNA component
	fComponentVolume = CalculateComponentVolume();

	// Number of fibers per voxel
	if (fBuildNucleus) {
		fNumFibers = 20;
	}
	else {
		fNumFibers = 1;
	}

	// If using a dose threshold to end simulation, convert to an energy threshold, which is
	// compared against after every event
	if (fUseDoseThreshold) {
		fEnergyThreshold = ConvertDoseThresholdToEnergy();
	}
	else {
		fEnergyThreshold = -1.;
	}

	// Determine order of magnitude of number of base pairs.
	// Set fParser parameters accordingly for use in ProcessHits()
	fNumBpMagnitude = CalculateIntegerMagnitude(fNumNucleosomePerFiber*fNumBpPerNucleosome);
	// fParserResidue = fNumBpMagnitude*10;
	// fParserStrand = fNumBpMagnitude*100;

	fThreadID = 0;
	fEventID = 0;

	//----------------------------------------------------------------------------------------------
	// Assign member variables to columns in the main output file. Contains DNA damage yields.
	//----------------------------------------------------------------------------------------------
	fNtuple->RegisterColumnI(&fThreadID, "Thread ID"); // Unique thread ID
	fNtuple->RegisterColumnI(&fEventID, "Event ID"); // Unique ID of primary particle / event / history
	fNtuple->RegisterColumnI(&fFiberID, "Fiber ID"); // Unique fiber ID
	fNtuple->RegisterColumnI(&fTotalSSB, "Single strand breaks"); // Number of SSB caused by this primary particle
	fNtuple->RegisterColumnI(&fTotalDSB, "Double strand breaks"); // Number of simple DSB caused by this primary particle
	fNtuple->RegisterColumnI(&fTotalBD, "Base damages"); // Number of BD caused by this primary particle
	if (fScoreClusters) {
		fNtuple->RegisterColumnI(&fTotalComplexDSB, "Complex DSBs"); // Number of Complex DSB caused by this primary particle
		fNtuple->RegisterColumnI(&fTotalNonDSBCluster, "Non-DSB clusters"); // Number of Non-DSB Clusters caused by this primary particle
	}
}


//--------------------------------------------------------------------------------------------------
// Destructor
//--------------------------------------------------------------------------------------------------
ScoreClusteredDNADamage::~ScoreClusteredDNADamage() {
}


//--------------------------------------------------------------------------------------------------
// Read in parameters from the Topas parameter file & save values in some member variables. Default
// values are provided in case the user does not specify a value.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::ResolveParams() {
	//----------------------------------------------------------------------------------------------
	// Scoring parameters
	//----------------------------------------------------------------------------------------------
	if ( fPm->ParameterExists(GetFullParmName("BasePairDistanceForDefiningDSB")) )
		fThresDistForDSB = fPm->GetIntegerParameter(GetFullParmName("BasePairDistanceForDefiningDSB"));
	else
		fThresDistForDSB = 10;

	if ( fPm->ParameterExists(GetFullParmName("EnergyThresholdForHavingSSB")) )
		fThresEdepForSSB = fPm->GetDoubleParameter(GetFullParmName("EnergyThresholdForHavingSSB"), "Energy");
	else
		fThresEdepForSSB = 17.5 * eV;

	if ( fPm->ParameterExists(GetFullParmName("EnergyThresholdForHavingBD")) )
		fThresEdepForBD = fPm->GetDoubleParameter(GetFullParmName("EnergyThresholdForHavingBD"), "Energy");
	else
		fThresEdepForBD = 17.5 * eV;

	if ( fPm->ParameterExists(GetFullParmName("BasePairDistanceForDefiningCluster")) )
		fThresDistForCluster = fPm->GetIntegerParameter(GetFullParmName("BasePairDistanceForDefiningCluster"));
	else
		fThresDistForCluster = 40;

	//----------------------------------------------------------------------------------------------
	// Score direct or indirect damage
	//----------------------------------------------------------------------------------------------
	if ( fPm->ParameterExists(GetFullParmName("IncludeDirectDamage")))
		fIncludeDirectDamage = fPm->GetBooleanParameter(GetFullParmName("IncludeDirectDamage"));
	else
		fIncludeDirectDamage = true;
	if ( fPm->ParameterExists(GetFullParmName("IncludeIndirectDamage")))
		fIncludeIndirectDamage = fPm->GetBooleanParameter(GetFullParmName("IncludeIndirectDamage"));
	else
		fIncludeIndirectDamage = true;

	//----------------------------------------------------------------------------------------------
	// Score damage clusters or not
	//----------------------------------------------------------------------------------------------
	if ( fPm->ParameterExists(GetFullParmName("ScoreClusters")))
		fScoreClusters = fPm->GetBooleanParameter(GetFullParmName("ScoreClusters"));
	else
		fScoreClusters = true;

	//----------------------------------------------------------------------------------------------
	// Specify whether to record damage on a per-event basis (vs. per-run basis)
	//----------------------------------------------------------------------------------------------
	if ( fPm->ParameterExists(GetFullParmName("RecordDamagePerEvent")))
		fRecordDamagePerEvent = fPm->GetBooleanParameter(GetFullParmName("RecordDamagePerEvent"));
	else
		fRecordDamagePerEvent = false;

	//----------------------------------------------------------------------------------------------
	// Specify whether to report damage on a per-fibre basis (vs. grouping together). Note that
	// damage in separate fibers are never considered together when clustering damage.
	//----------------------------------------------------------------------------------------------
	if ( fPm->ParameterExists(GetFullParmName("RecordDamagePerFiber")))
		fRecordDamagePerFiber = fPm->GetBooleanParameter(GetFullParmName("RecordDamagePerFiber"));
	else
		fRecordDamagePerFiber = false;

	//----------------------------------------------------------------------------------------------
	// Files for outputting run information & clustered damage details
	//----------------------------------------------------------------------------------------------
	if ( fPm->ParameterExists(GetFullParmName("FileRunSummary")))
		fFileRunSummary = fPm->GetStringParameter(GetFullParmName("FileRunSummary"));
	else
		fFileRunSummary = "output_run_summary";

	if ( fPm->ParameterExists(GetFullParmName("FileComplexDSB")))
		fFileComplexDSB = fPm->GetStringParameter(GetFullParmName("FileComplexDSB"));
	else
		fFileComplexDSB = "output_complex_dsb_specs";

	if ( fPm->ParameterExists(GetFullParmName("FileNonDSBCluster")))
		fFileNonDSBCluster = fPm->GetStringParameter(GetFullParmName("FileNonDSBCluster"));
	else
		fFileNonDSBCluster = "output_non_dsb_cluster_specs";

	//----------------------------------------------------------------------------------------------
	// Specify whether to output headers for data files or not
	//----------------------------------------------------------------------------------------------
	if ( fPm->ParameterExists(GetFullParmName("OutputHeaders")))
		fOutputHeaders = fPm->GetBooleanParameter(GetFullParmName("OutputHeaders"));
	else
		fOutputHeaders = true;

	//----------------------------------------------------------------------------------------------
	// Parameters to handle stopping simulation & scoring when dose threshold is met
	//----------------------------------------------------------------------------------------------
	if (fPm->ParameterExists(GetFullParmName("UseDoseThreshold")))
		fUseDoseThreshold = fPm->GetBooleanParameter(GetFullParmName("UseDoseThreshold"));
	else
		fUseDoseThreshold = false;

	if (fPm->ParameterExists(GetFullParmName("DoseThreshold")))
		fDoseThreshold = fPm->GetDoubleParameter(GetFullParmName("DoseThreshold"), "Dose");
	else
		fDoseThreshold = -1.;

	//----------------------------------------------------------------------------------------------
	// Geometry parameters
	//----------------------------------------------------------------------------------------------
	if (fPm->ParameterExists(GetFullParmName("BuildNucleus")))
		fBuildNucleus = fPm->GetBooleanParameter(GetFullParmName("BuildNucleus"));
	else
		fBuildNucleus = false;

	if (fPm->ParameterExists(GetFullParmName("NumVoxelsPerSide")))
		fNumVoxelsPerSide = fPm->GetIntegerParameter(GetFullParmName("NumVoxelsPerSide"));
	else
		fNumVoxelsPerSide = 1;

	if (fPm->ParameterExists(GetFullParmName("VoxelSideLength")))
		fVoxelSideLength = fPm->GetDoubleParameter(GetFullParmName("VoxelSideLength"),"Length");
	else
		fVoxelSideLength = 150*nm;

	if ( fPm->ParameterExists(GetFullParmName("DnaNumNucleosomePerFiber")) )
		fNumNucleosomePerFiber = fPm->GetIntegerParameter(GetFullParmName("DnaNumNucleosomePerFiber"));
	else
		fNumNucleosomePerFiber = 0;

	if ( fPm->ParameterExists(GetFullParmName("DnaNumBpPerNucleosome")) )
		fNumBpPerNucleosome = fPm->GetIntegerParameter(GetFullParmName("DnaNumBpPerNucleosome"));
	else
		fNumBpPerNucleosome = 0;

	//----------------------------------------------------------------------------------------------
	// Chemistry parameters
	//----------------------------------------------------------------------------------------------
	// MoleculeIDs
	fMoleculeID_OH = G4MoleculeTable::Instance()->GetConfiguration("OH")->GetMoleculeID();
	fMoleculeID_OHm = G4MoleculeTable::Instance()->GetConfiguration("OHm")->GetMoleculeID();
	fMoleculeID_e_aq = G4MoleculeTable::Instance()->GetConfiguration("e_aq")->GetMoleculeID();
	fMoleculeID_H2 = G4MoleculeTable::Instance()->GetConfiguration("H2")->GetMoleculeID();
	fMoleculeID_H3Op = G4MoleculeTable::Instance()->GetConfiguration("H3Op")->GetMoleculeID();
	fMoleculeID_H = G4MoleculeTable::Instance()->GetConfiguration("H")->GetMoleculeID();
	fMoleculeID_H2O2 = G4MoleculeTable::Instance()->GetConfiguration("H2O2")->GetMoleculeID();

	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/SSB/OH")))
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_OH, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/SSB/OH")));
	else
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_OH, 0.4);
	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/BD/OH")))
		fMoleculeDamageProb_BD.emplace(fMoleculeID_OH, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/BD/OH")));
	else
		fMoleculeDamageProb_BD.emplace(fMoleculeID_OH, 0.0);

	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/SSB/OHm")))
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_OHm, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/SSB/OHm")));
	else
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_OHm, 0.0);
	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/BD/OHm")))
		fMoleculeDamageProb_BD.emplace(fMoleculeID_OHm, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/BD/OHm")));
	else
		fMoleculeDamageProb_BD.emplace(fMoleculeID_OHm, 0.0);

	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/SSB/e_aq")))
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_e_aq, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/SSB/e_aq")));
	else
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_e_aq, 0.0);
	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/BD/e_aq")))
		fMoleculeDamageProb_BD.emplace(fMoleculeID_e_aq, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/BD/e_aq")));
	else
		fMoleculeDamageProb_BD.emplace(fMoleculeID_e_aq, 0.0);

	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/SSB/H2")))
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_H2, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/SSB/H2")));
	else
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_H2, 0.0);
	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/BD/H2")))
		fMoleculeDamageProb_BD.emplace(fMoleculeID_H2, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/BD/H2")));
	else
		fMoleculeDamageProb_BD.emplace(fMoleculeID_H2, 0.0);

	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/SSB/H3Op")))
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_H3Op, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/SSB/H3Op")));
	else
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_H3Op, 0.0);
	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/BD/H3Op")))
		fMoleculeDamageProb_BD.emplace(fMoleculeID_H3Op, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/BD/H3Op")));
	else
		fMoleculeDamageProb_BD.emplace(fMoleculeID_H3Op, 0.0);

	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/SSB/H")))
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_H, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/SSB/H")));
	else
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_H, 0.0);
	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/BD/H")))
		fMoleculeDamageProb_BD.emplace(fMoleculeID_H, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/BD/H")));
	else
		fMoleculeDamageProb_BD.emplace(fMoleculeID_H, 0.0);

	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/SSB/H2O2")))
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_H2O2, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/SSB/H2O2")));
	else
		fMoleculeDamageProb_SSB.emplace(fMoleculeID_H2O2, 0.0);
	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/BD/H2O2")))
		fMoleculeDamageProb_BD.emplace(fMoleculeID_H2O2, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/BD/H2O2")));
	else
		fMoleculeDamageProb_BD.emplace(fMoleculeID_H2O2, 0.0);

// ========================================================================================================================================
// ========================================================================================================================================
	// TODO: check if Chemistry is Extended
	// try (G4MoleculeTable::Instance()->GetConfiguration("HO2") != NULL) {
	// 	fMoleculeID_HO2 = G4MoleculeTable::Instance()->GetConfiguration("HO2")->GetMoleculeID();
	//
	// 	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/SSB/HO2")))
	// 		fMoleculeDamageProb_SSB.emplace(fMoleculeID_HO2, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/SSB/HO2")));
	// 	else
	// 		fMoleculeDamageProb_SSB.emplace(fMoleculeID_HO2, 0.0);
	//
	// 	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/BD/HO2")))
	// 		fMoleculeDamageProb_BD.emplace(fMoleculeID_HO2, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/BD/HO2")));
	// 	else
	// 		fMoleculeDamageProb_BD.emplace(fMoleculeID_HO2, 0.0);
	// }

	// if (G4MoleculeTable::Instance()->GetConfiguration("HO2m") != NULL) {
	// 	fMoleculeID_HO2m = G4MoleculeTable::Instance()->GetConfiguration("HO2m")->GetMoleculeID();
	//
	// 	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/SSB/HO2m")))
	// 		fMoleculeDamageProb_SSB.emplace(fMoleculeID_HO2m, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/SSB/HO2m")));
	// 	else
	// 		fMoleculeDamageProb_SSB.emplace(fMoleculeID_HO2m, 0.0);
	//
	// 	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/BD/HO2m")))
	// 		fMoleculeDamageProb_BD.emplace(fMoleculeID_HO2m, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/BD/HO2m")));
	// 	else
	// 		fMoleculeDamageProb_BD.emplace(fMoleculeID_HO2m, 0.0);
	// }
	//
	// if (G4MoleculeTable::Instance()->GetConfiguration("O2") != NULL) {
	// 	fMoleculeID_O2 = G4MoleculeTable::Instance()->GetConfiguration("O2")->GetMoleculeID();
	//
	// 	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/SSB/O2")))
	// 		fMoleculeDamageProb_SSB.emplace(fMoleculeID_O2, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/SSB/O2")));
	// 	else
	// 		fMoleculeDamageProb_SSB.emplace(fMoleculeID_O2, 0.0);
	//
	// 	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/BD/O2")))
	// 		fMoleculeDamageProb_BD.emplace(fMoleculeID_O2, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/BD/O2")));
	// 	else
	// 		fMoleculeDamageProb_BD.emplace(fMoleculeID_O2, 0.0);
	// }
	//
	// if (G4MoleculeTable::Instance()->GetConfiguration("O2m") != NULL) {
	// 	fMoleculeID_O2m = G4MoleculeTable::Instance()->GetConfiguration("O2m")->GetMoleculeID();
	//
	// 	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/SSB/O2m")))
	// 		fMoleculeDamageProb_SSB.emplace(fMoleculeID_O2m, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/SSB/O2m")));
	// 	else
	// 		fMoleculeDamageProb_SSB.emplace(fMoleculeID_O2m, 0.0);
	//
	// 	if (fPm->ParameterExists(GetFullParmName("DamageProbabilityOnInteraction/BD/O2m")))
	// 		fMoleculeDamageProb_BD.emplace(fMoleculeID_O2m, fPm->GetUnitlessParameter(GetFullParmName("DamageProbabilityOnInteraction/BD/O2m")));
	// 	else
	// 		fMoleculeDamageProb_BD.emplace(fMoleculeID_O2m, 0.0);
	// }

// ========================================================================================================================================
// ========================================================================================================================================

	//----------------------------------------------------------------------------------------------
	// Material of DNA residue volumes in which to score
	//----------------------------------------------------------------------------------------------
	G4String DNAMaterialName = "G4_WATER";
	if ( fPm->ParameterExists(GetFullParmName("DNAMaterialName")))
		DNAMaterialName = fPm->GetStringParameter(GetFullParmName("DNAMaterialName"));
	fDNAMaterial = GetMaterial(DNAMaterialName);

	//----------------------------------------------------------------------------------------------
	// Parameters not specific to this extension. Can't/don't need to use GetFullParmName()
	//----------------------------------------------------------------------------------------------
	fNumberOfThreads = fPm->GetIntegerParameter("Ts/NumberOfThreads");
}


//--------------------------------------------------------------------------------------------------
// Erase contents of various output files (and their corresponding header files)
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::ClearOutputFiles() {
	std::ofstream fileToClear;

	// Run summary
	fileToClear.open(fFileRunSummary+fOutFileExtension, std::ofstream::trunc);
	fileToClear.close();

	// Complex DSB
	fileToClear.open(fFileComplexDSB+fOutFileExtension, std::ofstream::trunc);
	fileToClear.close();

	// Non-DSB clusters
	fileToClear.open(fFileNonDSBCluster+fOutFileExtension, std::ofstream::trunc);
	fileToClear.close();

	// Headers
	if (fOutputHeaders) {
		fileToClear.open(fFileRunSummary+fOutHeaderExtension, std::ofstream::trunc);
		fileToClear.close();
		fileToClear.open(fFileComplexDSB+fOutHeaderExtension, std::ofstream::trunc);
		fileToClear.close();
		fileToClear.open(fFileNonDSBCluster+fOutHeaderExtension, std::ofstream::trunc);
		fileToClear.close();
	}
}


//--------------------------------------------------------------------------------------------------
// Calculate cubic volume of component attached to the scorer. Use parameter values to perform
// calculation according to the shape of the volume
//--------------------------------------------------------------------------------------------------
G4double ScoreClusteredDNADamage::CalculateComponentVolume() {
	G4double componentVolume = pow(2*fVoxelSideLength*fNumVoxelsPerSide,3);

	return componentVolume;
}

//--------------------------------------------------------------------------------------------------
// Handle conversion of dose threshold to energy threshold using the cubic volume and density of
// the geometry component. The resulting energy threshold is divided evenly amongst the worker
// threads. E.g. If total energy threshold is 1 MeV and there are 4 worker threads. The returned
// value is 0.25 MeV.
//--------------------------------------------------------------------------------------------------
G4double ScoreClusteredDNADamage::ConvertDoseThresholdToEnergy() {
	G4double energyThreshold;

	// Throw error if no dose threshold provided
	if (fDoseThreshold < 0) {
		G4cerr << "Topas is exiting due to a serious error in the scoring parameter DoseThreshold." << G4endl;
		G4cerr << "No valid dose threshold has been provided." << G4endl;
		fPm->AbortSession(1);
	}

	// Convert dose to energy thershold (Treating whole volume as same material)
	energyThreshold = GetMaterial("G4_WATER")->GetDensity() * fComponentVolume * fDoseThreshold;
	// energyThreshold = fDNAMaterial->GetDensity() * volume * fDoseThreshold;

	// Divide by number of threads
	if (fNumberOfThreads > 1){
		energyThreshold /= fNumberOfThreads;
	}

	// Output information
	// G4cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << G4endl;
	// G4cout << "density = " << fDNAMaterial->GetDensity()/(kg/m3) << G4endl;
	// G4cout << "volume = " << volume/m3 << G4endl;
	// G4cout << "fDoseThreshold = " << fDoseThreshold/gray << G4endl;
	// G4cout << "energyThreshold = " << energyThreshold/eV << G4endl;
	// G4cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << G4endl;

	return energyThreshold;
}

//--------------------------------------------------------------------------------------------------
// The ProcessHits method is called for every hit in the sensitive detector volume, i.e. when an
// interaction occurs in a sensitive volume (may or may not be energy deposit). In this case simply
// record energy deposited in a given bp index, in either the backbone or base of either  strand # 1
// or strand #2 (4 separate energy deposition maps). Multiple energy depositions occurring in the
// same volume are added together cumulatively
//
// Note that the copy number of the volume is used to determine in which volume the energy
// deposition took place. This is faster than using string comparisons. Tis method is only called
// for energy depositions in the sensitive volumes (i.e. residues) by using material filtering, as
// defined in the parameter file with "OnlyIncludeIfInMaterial" parameter.
//--------------------------------------------------------------------------------------------------
G4bool ScoreClusteredDNADamage::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	fNumProcessHitsCalls++;
	G4double edep = aStep->GetTotalEnergyDeposit(); // In eV;
	fTotalEdep += edep; // running sum of energy deposition in entire volume

	// Check if volume is histone
	G4int volCopyNo = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo() / 1000000;
	G4int trackID = aStep->GetTrack()->GetTrackID();
	if (volCopyNo == 2 && trackID < 0){ // histone and chemical track

		G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
		G4String volumeName = touchable->GetVolume()->GetName();
		G4cout << " volumeName: " << volumeName << G4endl;
		G4String processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
		G4cout << "\tprocessName: " << processName << G4endl;
		G4int volumeCopyNo = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
		G4cout << "\tvolumeCopyNo: " << volumeCopyNo << G4endl;

		aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		G4cout << "\tTrack KILLED! trackID: " << trackID << G4endl;
		return false;
	}

	// Material filtering (only proceed if in sensitive DNA volumes)
	G4Material* material = aStep->GetPreStepPoint()->GetMaterial();
	if ( material != fDNAMaterial) {
		return false;
	}

	//----------------------------------------------------------------------------------------------
	// If this hit deposits energy (in sensitive DNA volume), update the appropriate energy deposition
	// map
	//----------------------------------------------------------------------------------------------
	if (fIncludeDirectDamage && edep > 0) {
		// G4Touchable provides access to parent volumes, etc.
		G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
		// Determine unique voxel ID number using replica ID of parent volumes
		if (fBuildNucleus) {
			G4int voxIDZ = touchable->GetReplicaNumber(fParentIndexVoxelZ);
			G4int voxIDX = touchable->GetReplicaNumber(fParentIndexVoxelX);
			G4int voxIDY = touchable->GetReplicaNumber(fParentIndexVoxelY);
			fVoxelID = voxIDZ + (fNumVoxelsPerSide*voxIDX) + (fNumVoxelsPerSide*fNumVoxelsPerSide*voxIDY);
		}
		// Determine unique fiber ID number using copy ID of parent
		if (fNumFibers > 1)
			fFiberID = touchable->GetCopyNumber(fParentIndexFiber);
		//------------------------------------------------------------------------------------------
		// Determine the indices defining the volume in which energy was deposited by parsing the
		// copy number of the Physical Volume.
		//------------------------------------------------------------------------------------------
		// G4int volID = touchable->GetVolume()->GetCopyNo();
		G4int volID = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
		G4int num_strand = volID / 1000000;
		G4int num_res = (volID - (num_strand*1000000)) / 100000;
		G4int num_nucleotide = volID - (num_strand*1000000) - (num_res*100000);

		//------------------------------------------------------------------------------------------
		// Use the DNA strand ID, residue ID, and nucleotide ID to increment the energy deposited
		// in the appropriate energy deposition map. Maps are indexed as follows:
		// First index specifies the voxel
		// Second index specifies DNA fibre
		// Third index specifies the bp index
		//------------------------------------------------------------------------------------------
		if ( num_strand == 0 ){ // first strand
			if (num_res == fVolIdPhosphate || num_res == fVolIdDeoxyribose){
				fMapEdepStrand1Backbone[fVoxelID][fFiberID][num_nucleotide] += edep;
			}
			else if (num_res == fVolIdBase) {
				fMapEdepStrand1Base[fVoxelID][fFiberID][num_nucleotide] += edep;
			}
		}
		else{ // second strand
			if (num_res == fVolIdPhosphate || num_res == fVolIdDeoxyribose){
				fMapEdepStrand2Backbone[fVoxelID][fFiberID][num_nucleotide] += edep;
			}
			else if (num_res == fVolIdBase) {
				fMapEdepStrand2Base[fVoxelID][fFiberID][num_nucleotide] += edep;
			}
		}
		return true;
	}

	// Indirect damage starting here
	if (fIncludeIndirectDamage && trackID < 0) { // chemical tracks
		// G4Touchable provides access to parent volumes, etc.
		G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
		// Determine unique voxel ID number using replica ID of parent volumes
		if (fBuildNucleus) {
			G4int voxIDZ = touchable->GetReplicaNumber(fParentIndexVoxelZ);
			G4int voxIDX = touchable->GetReplicaNumber(fParentIndexVoxelX);
			G4int voxIDY = touchable->GetReplicaNumber(fParentIndexVoxelY);
			fVoxelID = voxIDZ + (fNumVoxelsPerSide*voxIDX) + (fNumVoxelsPerSide*fNumVoxelsPerSide*voxIDY);
		}
		// Determine unique fiber ID number using copy ID of parent
		if (fNumFibers > 1)
			fFiberID = touchable->GetCopyNumber(fParentIndexFiber);
		//------------------------------------------------------------------------------------------
		// Determine the indices defining the volume in which energy was deposited by parsing the
		// copy number of the Physical Volume.
		//------------------------------------------------------------------------------------------
		// G4int volID = touchable->GetVolume()->GetCopyNo();
		G4int volID = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
		G4int num_strand = volID / 1000000;
		G4int num_res = (volID - (num_strand*1000000)) / 100000;
		G4int num_nucleotide = volID - (num_strand*1000000) - (num_res*100000);

		// Get molecule info
		G4Track* aTrack = aStep->GetTrack();
		G4int moleculeID = GetMolecule(aTrack)->GetMoleculeID();

		// Determine if damage is inflicted
		G4bool damaged = IsDamageInflicted(moleculeID, num_res);

		// Record nucleotide damage.
		if (damaged) {
			G4cout << "\tDAMAGE INDUCED!" << G4endl;
			G4cout << "\teventID: " << GetEventID() << G4endl;
			G4cout << "\tthreadID: " << G4Threading::G4GetThreadId() << G4endl;

			// Do not record if backbone or base has already been "damaged" previously via indirect actions
			if ( num_strand == 0 ) { // first strand
				if (num_res == fVolIdPhosphate || num_res == fVolIdDeoxyribose) {
					if (!IsElementInVector(num_nucleotide, fIndicesSSB1_indirect)) {
						fMapIndDamageStrand1Backbone[fVoxelID][fFiberID].push_back(num_nucleotide);
						G4cout << "\tRECORDED SSB1! num_nucleotide: " << num_nucleotide << G4endl;
						G4cout << "\tNumber of recorded indices in vector: " << fMapIndDamageStrand1Backbone[fVoxelID][fFiberID].size() << G4endl;
					}
					else {
						G4cout << "\tALREADY DAMAGED SSB1! num_nucleotide: " << num_nucleotide << G4endl;
						return false;
					}
				}
				else if (num_res == fVolIdBase) {
					if (!IsElementInVector(num_nucleotide, fIndicesBD1_indirect)) {
						fMapIndDamageStrand1Base[fVoxelID][fFiberID].push_back(num_nucleotide);
						G4cout << "\tRECORDED BD1! num_nucleotide: " << num_nucleotide << G4endl;
						G4cout << "\tNumber of recorded indices in vector: " << fMapIndDamageStrand1Base[fVoxelID][fFiberID].size() << G4endl;
					}
					else {
						G4cout << "\tALREADY DAMAGED BD1! num_nucleotide: " << num_nucleotide << G4endl;
						return false;
					}
				}
			}
			else { // second strand
				if (num_res == fVolIdPhosphate || num_res == fVolIdDeoxyribose) {
					if (!IsElementInVector(num_nucleotide, fIndicesSSB2_indirect)) {
						fMapIndDamageStrand2Backbone[fVoxelID][fFiberID].push_back(num_nucleotide);
						G4cout << "\tRECORDED SSB2! num_nucleotide: " << num_nucleotide << G4endl;
						G4cout << "\tNumber of recorded indices in vector: " << fMapIndDamageStrand2Backbone[fVoxelID][fFiberID].size() << G4endl;
					}
					else {
						G4cout << "\tALREADY DAMAGED SSB2! num_nucleotide: " << num_nucleotide << G4endl;
						return false;
					}
				}
				else if (num_res == fVolIdBase) {
					if (!IsElementInVector(num_nucleotide, fIndicesBD2_indirect)) {
						fMapIndDamageStrand2Base[fVoxelID][fFiberID].push_back(num_nucleotide);
						G4cout << "\tRECORDED BD2! num_nucleotide: " << num_nucleotide << G4endl;
						G4cout << "\tNumber of recorded indices in vector: " << fMapIndDamageStrand2Base[fVoxelID][fFiberID].size() << G4endl;
					}
					else {
						G4cout << "\tALREADY DAMAGED BD2! num_nucleotide: " << num_nucleotide << G4endl;
						return false;
					}
				}
			}
		} // damaged
		// Kill track whether the interaction induced damage or not
		aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		G4cout << "\tTrack KILLED! trackID: " << trackID << G4endl;
		// G4cout << "\tProcessHits calls: " << fNumProcessHitsCalls << G4endl;
		return true;
	} // negative trackID

	return false;
}


//--------------------------------------------------------------------------------------------------
// This helper method checks whether an element is in a vector.
//--------------------------------------------------------------------------------------------------
G4bool ScoreClusteredDNADamage::IsElementInVector(G4int pElement, std::vector<G4int> pVector) {
	if (std::find(pVector.begin(), pVector.end(), pElement) != pVector.end())
		return true;
	else
		return false;
}


//--------------------------------------------------------------------------------------------------
// This helper method checks whether an indirect damage has been induced given the moleculeID
//--------------------------------------------------------------------------------------------------
G4bool ScoreClusteredDNADamage::IsDamageInflicted(G4int pMoleculeID, G4int pDNAVolumeID) {
	G4float prob_damage = 0.;

	if ( fMoleculeDamageProb_SSB.find(pMoleculeID) == fMoleculeDamageProb_SSB.end() && fMoleculeDamageProb_BD.find(pMoleculeID) == fMoleculeDamageProb_BD.end()) {
  	// pMoleculeID not found in fMoleculeDamageProb_SSB and fMoleculeDamageProb_BD
		G4cout << "\tmoleculeID NOT LISTED: " << pMoleculeID << G4endl;
		return false;
	}

	if (pDNAVolumeID == fVolIdPhosphate || pDNAVolumeID == fVolIdDeoxyribose)
		prob_damage = fMoleculeDamageProb_SSB[pMoleculeID]; // returns 0 if pMoleculeID is not a "key" in the map
	if (pDNAVolumeID == fVolIdBase)
		prob_damage = fMoleculeDamageProb_BD[pMoleculeID]; // returns 0 if pMoleculeID is not a "key" in the map

	if (prob_damage == 0.)
		return false;

	// Generate random float between 0 and 1
	G4float prob_random = (float) std::rand()/RAND_MAX;

	G4cout << "\tmoleculeID: " << pMoleculeID << G4endl;
	G4cout << "\tprob_damage: " << prob_damage << G4endl;
	G4cout << "\tprob_random: " << prob_random << G4endl;

	if (prob_random <= prob_damage)
		return true;
	else
		return false;
}


//--------------------------------------------------------------------------------------------------
// This method is called at the end of the event (primary particle & its secondaries).
//
// If doing event-by-event scoring: process the energy deposoition maps to determine DNA damage
// yields and output the results.
//
// If using a dose threshold on the run, check to see if threshold has been exceeded & abort this
// worker thread if so.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::UserHookForEndOfEvent() {
	fEventID = GetEventID();
	fThreadID = G4Threading::G4GetThreadId();

	fNumEvents++;

	// Analyze damage if doing event-by-event scoring
	if (fRecordDamagePerEvent) {
		RecordDamage();
		if (fScoreClusters) {
			OutputComplexDSBToFile();
			OutputNonDSBClusterToFile();
		}
		ResetMemberVariables(); // Necessary to reset variables before proceeding to next event
	}

	// Check if dose threshold has been met
	if (fUseDoseThreshold && fTotalEdep > fEnergyThreshold) {
		G4cout << "Aborting worker #" << G4Threading::G4GetThreadId() << " because dose threshold has been met" << G4endl;
		G4RunManager::GetRunManager()->AbortRun(true);
	}
}

//--------------------------------------------------------------------------------------------------
// This method is called at the end of the run (all primary particles). Only called by the master
// thread, not the worker threads.
//
// If scoring over the whole run: process the energy deposoition maps (combined over all worker
// threads) to determine DNA yields and output the results.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::UserHookForEndOfRun() {
	// fEventID = GetEventID();
	fThreadID = G4Threading::G4GetThreadId();

	OutputRunSummaryToFile();
	G4cout << "Run summary has been written to: " << fFileRunSummary << G4endl;

	// Analyze damage if scoring over the whole run
	if (!fRecordDamagePerEvent) {
		fEventID = fAggregateValueIndicator;
		RecordDamage();
		if (fScoreClusters) {
			OutputComplexDSBToFile();
			OutputNonDSBClusterToFile();
		}
	}

	if (fScoreClusters) {
		G4cout << "Complex DSB details have been written to: " << fFileComplexDSB << G4endl;
		G4cout << "Non-DSB cluster details have been written to: " << fFileNonDSBCluster << G4endl;
	}
}


//--------------------------------------------------------------------------------------------------
// This method outputs the details of the run to a header file and a data file.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::OutputRunSummaryToFile() {
	//----------------------------------------------------------------------------------------------
	// Header file
	//----------------------------------------------------------------------------------------------
	if (fOutputHeaders) {
		std::ofstream outHeader;
		G4String headerFileName = fFileRunSummary + fOutHeaderExtension;

		outHeader.open(headerFileName, std::ofstream::trunc);

		// Catch file I/O error
		if (!outHeader.good()) {
			G4cerr << "Topas is exiting due to a serious error in file output." << G4endl;
			G4cerr << "Output file: " << headerFileName << " cannot be opened" << G4endl;
			fPm->AbortSession(1);
		}

		outHeader << "# events" << fDelimiter;
		outHeader << "Dose (Gray)" << fDelimiter;
		outHeader << "Energy (eV)" << G4endl;
		outHeader.close();
	}

	//----------------------------------------------------------------------------------------------
	// Data file
	//----------------------------------------------------------------------------------------------
	G4String outputFileName = fFileRunSummary + fOutFileExtension;
	std::ofstream outFile(outputFileName, std::ios_base::app);

	// Catch file I/O error
	if (!outFile.good()) {
		G4cerr << "Topas is exiting due to a serious error in file output." << G4endl;
		G4cerr << "Output file: " << outputFileName << " cannot be opened" << G4endl;
		fPm->AbortSession(1);
	}

	G4double doseDep = fTotalEdep / GetMaterial("G4_WATER")->GetDensity() / fComponentVolume;

	outFile << fNumEvents << fDelimiter;
	outFile << doseDep/gray << fDelimiter;
	outFile << fTotalEdep/eV << G4endl;

	outFile.close();
}


//--------------------------------------------------------------------------------------------------
// This method outputs the details of scored Complex DSBs to a header file and data file. Each line
// in the data file contains information for a single cluster.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::OutputComplexDSBToFile() {
	//----------------------------------------------------------------------------------------------
	// Header file
	//----------------------------------------------------------------------------------------------
	if (fOutputHeaders) {
		std::ofstream outHeader;
		G4String headerFileName = fFileComplexDSB + fOutHeaderExtension;

		outHeader.open(headerFileName, std::ofstream::trunc);

		// Catch file I/O error
		if (!outHeader.good()) {
			G4cerr << "Topas is exiting due to a serious error in file output." << G4endl;
			G4cerr << "Output file: " << headerFileName << " cannot be opened" << G4endl;
			fPm->AbortSession(1);
		}

		outHeader << "Size (bp)" << fDelimiter;
		outHeader << "Total # of damages" << fDelimiter;
		outHeader << "# of SSB" << fDelimiter;
		outHeader << "# of BD" << fDelimiter;
		outHeader << "# of DSB" << G4endl;
		outHeader.close();
	}

	//----------------------------------------------------------------------------------------------
	// Data file
	//----------------------------------------------------------------------------------------------
	G4String outputFileName = fFileComplexDSB + fOutFileExtension;
	std::ofstream outFile(outputFileName, std::ios_base::app);

	// Catch file I/O error
	if (!outFile.good()) {
		G4cerr << "Topas is exiting due to a serious error in file output." << G4endl;
		G4cerr << "Output file: " << outputFileName << " cannot be opened" << G4endl;
		fPm->AbortSession(1);
	}

	// Record data
	for (G4int i = 0; i < fComplexDSBSizes.size(); i++) {
		outFile << fComplexDSBSizes[i] << fDelimiter;
		outFile << fComplexDSBNumDamage[i] << fDelimiter;
		outFile << fComplexDSBNumSSB[i] << fDelimiter;
		outFile << fComplexDSBNumBD[i] << fDelimiter;
		outFile << fComplexDSBNumDSB[i] << G4endl;
	}
	outFile.close();
}


//--------------------------------------------------------------------------------------------------
// This method outputs the details of scored Non-DSB clusters to a header file and data file. Each
// line in the data file contains information for a single cluster.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::OutputNonDSBClusterToFile() {
	//----------------------------------------------------------------------------------------------
	// Header file
	//----------------------------------------------------------------------------------------------
	if (fOutputHeaders) {
		std::ofstream outHeader;
		G4String headerFileName = fFileNonDSBCluster + fOutHeaderExtension;

		outHeader.open(headerFileName, std::ofstream::trunc);

		// Catch file I/O error
		if (!outHeader.good()) {
			G4cerr << "Topas is exiting due to a serious error in file output." << G4endl;
			G4cerr << "Output file: " << headerFileName << " cannot be opened" << G4endl;
			fPm->AbortSession(1);
		}

		outHeader << "Size (bp)" << fDelimiter;
		outHeader << "Total # of damages" << fDelimiter;
		outHeader << "# of SSB" << fDelimiter;
		outHeader << "# of BD" << G4endl;
		outHeader.close();
	}

	//----------------------------------------------------------------------------------------------
	// Data file
	//----------------------------------------------------------------------------------------------
	G4String outputFileName = fFileNonDSBCluster + fOutFileExtension;
	std::ofstream outFile(outputFileName, std::ios_base::app);

	// Catch file I/O error
	if (!outFile.good()) {
		G4cerr << "Topas is exiting due to a serious error in file output." << G4endl;
		G4cerr << "Output file: " << outputFileName << " cannot be opened" << G4endl;
		fPm->AbortSession(1);
	}

	// Record data
	for (G4int i = 0; i < fNonDSBClusterSizes.size(); i++) {
		outFile << fNonDSBClusterSizes[i] << fDelimiter;
		outFile << fNonDSBClusterNumDamage[i] << fDelimiter;
		outFile << fNonDSBClusterNumSSB[i] << fDelimiter;
		outFile << fNonDSBClusterNumBD[i] << G4endl;
	}
	outFile.close();
}


//--------------------------------------------------------------------------------------------------
// This method transfers information from worker threads to the master thread, which allows results
// be processed on a per-run basis. This method is called once per worker thread, at the end of that
// threads execution. The class members populated here are those of the master thread.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer)
{
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer); // run the parent version

	ScoreClusteredDNADamage* myWorkerScorer = dynamic_cast<ScoreClusteredDNADamage*>(workerScorer);

	// Absorb various worker thread data
	fTotalEdep += myWorkerScorer->fTotalEdep;
	fNumEvents += myWorkerScorer->fNumEvents;
	fNumProcessHitsCalls += myWorkerScorer->fNumProcessHitsCalls;

    // Absorb the energy deposition maps from this worker
    if (!fRecordDamagePerEvent) {
	    AbsorbDirDmgMapFromWorkerScorer(fMapEdepStrand1Backbone,myWorkerScorer->fMapEdepStrand1Backbone);
	    AbsorbDirDmgMapFromWorkerScorer(fMapEdepStrand2Backbone,myWorkerScorer->fMapEdepStrand2Backbone);
	    AbsorbDirDmgMapFromWorkerScorer(fMapEdepStrand1Base,myWorkerScorer->fMapEdepStrand1Base);
	    AbsorbDirDmgMapFromWorkerScorer(fMapEdepStrand2Base,myWorkerScorer->fMapEdepStrand2Base);

			AbsorbIndDmgMapFromWorkerScorer(fMapIndDamageStrand1Backbone,myWorkerScorer->fMapIndDamageStrand1Backbone);
			AbsorbIndDmgMapFromWorkerScorer(fMapIndDamageStrand2Backbone,myWorkerScorer->fMapIndDamageStrand2Backbone);
			AbsorbIndDmgMapFromWorkerScorer(fMapIndDamageStrand1Base,myWorkerScorer->fMapIndDamageStrand1Base);
			AbsorbIndDmgMapFromWorkerScorer(fMapIndDamageStrand2Base,myWorkerScorer->fMapIndDamageStrand2Base);
	}
}

//--------------------------------------------------------------------------------------------------
// This method transfers the contents of a map of energy depositions from the worker thread to the
// master thread
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::AbsorbDirDmgMapFromWorkerScorer(
	std::map<G4int,std::map<G4int,std::map<G4int, G4double>>> &masterMap,
	std::map<G4int,std::map<G4int,std::map<G4int, G4double>>> &workerMap)
{
	// Loop over all voxels in nucleus
	std::map<G4int,std::map<G4int,std::map<G4int, G4double>>>::iterator itVoxel = workerMap.begin();
	while (itVoxel != workerMap.end()) {
		G4int indexVoxel = itVoxel->first;
		std::map<G4int,std::map<G4int, G4double>> workerMapVoxel = itVoxel->second;

		// Loop over all fibers in voxel
		std::map<G4int,std::map<G4int, G4double>>::iterator itFiber = workerMapVoxel.begin();
		while (itFiber != workerMapVoxel.end()) {
			G4int indexFiber = itFiber->first;
			std::map<G4int, G4double> workerMapFiber = itFiber->second;

			// Loop over all base pairs fiber
			std::map<G4int, G4double>::iterator itBP = workerMapFiber.begin();
			while (itBP != workerMapFiber.end()) {
				G4int indexBP = itBP->first;
				G4double eDep = itBP->second;

				// Increment master thread energy map
				masterMap[indexVoxel][indexFiber][indexBP] += eDep;
				itBP++;
			}
			itFiber++;
		}
		itVoxel++;
	}
}


//--------------------------------------------------------------------------------------------------
// This method transfers the contents of a map of indirect hits indices from the worker thread to the
// master thread
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::AbsorbIndDmgMapFromWorkerScorer(
	std::map<G4int,std::map<G4int,std::vector<G4int>>> &masterMap,
	std::map<G4int,std::map<G4int,std::vector<G4int>>> &workerMap)
{
	// Loop over all voxels in nucleus
	std::map<G4int,std::map<G4int,std::vector<G4int>>>::iterator itVoxel = workerMap.begin();
	while (itVoxel != workerMap.end()) {
		G4int indexVoxel = itVoxel->first;
		std::map<G4int,std::vector<G4int>> workerMapVoxel = itVoxel->second;

		// Loop over all fibers in voxel
		std::map<G4int,std::vector<G4int>>::iterator itFiber = workerMapVoxel.begin();
		while (itFiber != workerMapVoxel.end()) {
			G4int indexFiber = itFiber->first;
			std::vector<G4int> workerMapFiber = itFiber->second;

			// Print1DVectorContents(workerMapFiber); // debugging

			// Loop over all base pairs fiber
			std::vector<G4int>::iterator itBP = workerMapFiber.begin();
			while (itBP != workerMapFiber.end()) {
				G4int indexBP = *itBP;

				// Increment master thread damage index map
				if (!IsElementInVector(indexBP, masterMap[indexVoxel][indexFiber]))
					masterMap[indexVoxel][indexFiber].push_back(indexBP);
				else
					fDoubleCountsII++;

				itBP++;
			}
			itFiber++;
		}
		itVoxel++;
	}
}


//--------------------------------------------------------------------------------------------------
// Process maps of energy depositions and record DNA damage yields to member variables. Damage is
// processed within a single DNA fiber at a time (i.e. damages in subsequent fibers are not
// processed together). Simple damages contained in aggregate damages (DSBs and clusters) are not
// included in their own counters (i.e. the 2 SSBs comprising a DSB do not count towards fTotalSSB).
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::RecordDamage() {
	// Include following line if want to create a fake, predefined energy map to validate scoring
	// CreateFakeEnergyMap();

	// Iterate over all voxels
	G4int numVoxels = pow(fNumVoxelsPerSide,3);
	for (G4int iVoxel = 0; iVoxel < numVoxels; iVoxel++) {
		fVoxelID = iVoxel;

		// Iterate over all DNA fibres
		for (G4int iFiber = 0; iFiber < fNumFibers; iFiber++) {
			fFiberID = iFiber;

			// As only a single fiber is processed at a time, copy the 1D maps for the current fiber
			// from the larger 3D maps for easier processing.
			fFiberMapEdepStrand1Backbone = fMapEdepStrand1Backbone[iVoxel][iFiber];
			fFiberMapEdepStrand2Backbone = fMapEdepStrand2Backbone[iVoxel][iFiber];
			fFiberMapEdepStrand1Base = fMapEdepStrand1Base[iVoxel][iFiber];
			fFiberMapEdepStrand2Base = fMapEdepStrand2Base[iVoxel][iFiber];

			// Determine yields of simple damages (SSB and BD) in both strands
			fIndicesSSB1_direct = RecordSimpleDamage(fThresEdepForSSB,fFiberMapEdepStrand1Backbone);
			fIndicesSSB2_direct = RecordSimpleDamage(fThresEdepForSSB,fFiberMapEdepStrand2Backbone);
			fIndicesBD1_direct = RecordSimpleDamage(fThresEdepForBD,fFiberMapEdepStrand1Base);
			fIndicesBD2_direct = RecordSimpleDamage(fThresEdepForBD,fFiberMapEdepStrand2Base);

			fIndicesSSB1_indirect = fMapIndDamageStrand1Backbone[iVoxel][iFiber];
			fIndicesSSB2_indirect = fMapIndDamageStrand2Backbone[iVoxel][iFiber];
			fIndicesBD1_indirect = fMapIndDamageStrand1Base[iVoxel][iFiber];
			fIndicesBD2_indirect = fMapIndDamageStrand2Base[iVoxel][iFiber];

			// Merge and resolve damage yields from direct and indirect damage
			if (fIncludeDirectDamage && fIncludeIndirectDamage){
				fIndicesSSB1 = MergeDamageIndices(fIndicesSSB1_direct, fIndicesSSB1_indirect);
				fIndicesSSB2 = MergeDamageIndices(fIndicesSSB2_direct, fIndicesSSB2_indirect);
				fIndicesBD1 = MergeDamageIndices(fIndicesBD1_direct, fIndicesBD1_indirect);
				fIndicesBD2 = MergeDamageIndices(fIndicesBD2_direct, fIndicesBD2_indirect);
			}
			else if (fIncludeIndirectDamage) {
				fIndicesSSB1 = fIndicesSSB1_indirect;
				fIndicesSSB2 = fIndicesSSB2_indirect;
				fIndicesBD1 = fIndicesBD1_indirect;
				fIndicesBD2 = fIndicesBD2_indirect;
			}
			else {
				fIndicesSSB1 = fIndicesSSB1_direct;
				fIndicesSSB2 = fIndicesSSB2_direct;
				fIndicesBD1 = fIndicesBD1_direct;
				fIndicesBD2 = fIndicesBD2_direct;
			}

			// debugging
			// G4cout << "\tSSB1-indirect: " << G4endl;
			// Print1DVectorContents(fIndicesSSB1_indirect);
			// G4cout << "\tSSB1-direct: " << G4endl;
			// Print1DVectorContents(fIndicesSSB1_direct);
			// G4cout << "\tSSB1-merged: " << G4endl;
			// Print1DVectorContents(fIndicesSSB1);
			// G4cout << G4endl;
			//
			// G4cout << "\tSSB2-indirect: " << G4endl;
			// Print1DVectorContents(fIndicesSSB2_indirect);
			// G4cout << "\tSSB2-direct: " << G4endl;
			// Print1DVectorContents(fIndicesSSB2_direct);
			// G4cout << "\tSSB2-merged: " << G4endl;
			// Print1DVectorContents(fIndicesSSB2);
			// G4cout << G4endl;
			//
			// G4cout << "\tBD1-indirect: " << G4endl;
			// Print1DVectorContents(fIndicesBD1_indirect);
			// G4cout << "\tBD1-direct: " << G4endl;
			// Print1DVectorContents(fIndicesBD1_direct);
			// G4cout << "\tBD1-merged: " << G4endl;
			// Print1DVectorContents(fIndicesBD1);
			// G4cout << G4endl;
			//
			// G4cout << "\tBD2-indirect: " << G4endl;
			// Print1DVectorContents(fIndicesBD2_indirect);
			// G4cout << "\tBD2-direct: " << G4endl;
			// Print1DVectorContents(fIndicesBD2_direct);
			// G4cout << "\tBD2-merged: " << G4endl;
			// Print1DVectorContents(fIndicesBD2);
			// G4cout << G4endl;

			// G4cout << "\tDouble counts DI: " << fDoubleCountsDI << G4endl;
			// G4cout << "\tDouble counts II: " << fDoubleCountsII << G4endl;
			// G4cout << "\tProcessHits calls: " << fNumProcessHitsCalls << G4endl;

			// Process SSBs in both strands to determine if there are any DSB
			fIndicesDSB = RecordDSB();

			// Determine the present total number of SSB, BD, and DSB. These may be later modified
			// by RecordClusteredDamage() if there are clusters of damage.
			fTotalSSB += fIndicesSSB1.size() + fIndicesSSB2.size();
			fTotalBD += fIndicesBD1.size() + fIndicesBD2.size();
			fTotalDSB += fIndicesDSB.size(); // Note is currently 2x number of DSB. Division by 2 is performed later

			// If recording clustered damage, combine all damages into a single, sequential vector
			// of damage that indicates the type and bp index. Then process this vector to determine
			// clustered damage yields
			if (fScoreClusters) {
				fIndicesSimple = CombineSimpleDamage();
				RecordClusteredDamage();
			}

			// If recording damage on a fiber-by-fiber basis, fill the output ntuple
			if (fRecordDamagePerFiber) {
				fTotalDSB = fTotalDSB/2;
				fNtuple->Fill(); // Move this to outside loop if aggregating over all fibres

				// Reset variables before next fibre (not aggregating over all fibres)
				ResetDamageCounterVariables();
			}
		}
	}

	// If recording damage aggregated over all fibers, fill the output ntuple
	if (!fRecordDamagePerFiber) {
		fTotalDSB = fTotalDSB/2;
		fFiberID = fAggregateValueIndicator;
		fNtuple->Fill(); // Move this to outside loop if aggregating over all fibres
	}
	// PrintDNADamageToConsole(); // debugging;
}


//--------------------------------------------------------------------------------------------------
// This method merges and resolves duplicates of the damage yields from direct and indirect damage.
//--------------------------------------------------------------------------------------------------
std::vector<G4int> ScoreClusteredDNADamage::MergeDamageIndices(std::vector<G4int> pDamageIndices1,
	std::vector<G4int> pDamageIndices2)
{
	std::vector<G4int> indicesDamage_merged;

	for (G4int element : pDamageIndices1) {
		if (!IsElementInVector(element, indicesDamage_merged))
			indicesDamage_merged.push_back(element);
		else
			fDoubleCountsDI++;
	}

	for (G4int element : pDamageIndices2) {
		if (!IsElementInVector(element, indicesDamage_merged))
			indicesDamage_merged.push_back(element);
		else
			fDoubleCountsDI++;
	}

	return indicesDamage_merged;
}


//--------------------------------------------------------------------------------------------------
// This method resets member variable values, which is necessary if processing damage on an
// event-by-event basis.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::ResetMemberVariables() {
	fFiberMapEdepStrand1Backbone.erase(fFiberMapEdepStrand1Backbone.begin(), fFiberMapEdepStrand1Backbone.end());
	fFiberMapEdepStrand2Backbone.erase(fFiberMapEdepStrand2Backbone.begin(), fFiberMapEdepStrand2Backbone.end());
	fFiberMapEdepStrand1Base.erase(fFiberMapEdepStrand1Base.begin(), fFiberMapEdepStrand1Base.end());
	fFiberMapEdepStrand2Base.erase(fFiberMapEdepStrand2Base.begin(), fFiberMapEdepStrand2Base.end());

	ResetDamageCounterVariables();

	fComplexDSBSizes.clear();
	fComplexDSBNumSSB.clear();
	fComplexDSBNumBD.clear();
	fComplexDSBNumDSB.clear();
	fComplexDSBNumDamage.clear();

	fNonDSBClusterSizes.clear();
	fNonDSBClusterNumSSB.clear();
	fNonDSBClusterNumBD.clear();
	fNonDSBClusterNumDamage.clear();

	fMapEdepStrand1Backbone.erase(fMapEdepStrand1Backbone.begin(), fMapEdepStrand1Backbone.end());
	fMapEdepStrand2Backbone.erase(fMapEdepStrand2Backbone.begin(), fMapEdepStrand2Backbone.end());
	fMapEdepStrand1Base.erase(fMapEdepStrand1Base.begin(), fMapEdepStrand1Base.end());
	fMapEdepStrand2Base.erase(fMapEdepStrand2Base.begin(), fMapEdepStrand2Base.end());

	fMapIndDamageStrand1Backbone.erase(fMapIndDamageStrand1Backbone.begin(), fMapIndDamageStrand1Backbone.end());
	fMapIndDamageStrand2Backbone.erase(fMapIndDamageStrand2Backbone.begin(), fMapIndDamageStrand2Backbone.end());
	fMapIndDamageStrand1Base.erase(fMapIndDamageStrand1Base.begin(), fMapIndDamageStrand1Base.end());
	fMapIndDamageStrand2Base.erase(fMapIndDamageStrand2Base.begin(), fMapIndDamageStrand2Base.end());
}


//--------------------------------------------------------------------------------------------------
// This method resets variables that count the yields for various types of DNA damage.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::ResetDamageCounterVariables() {
	fTotalSSB = 0;
	fTotalBD = 0;
	fTotalDSB = 0;
	fTotalComplexDSB = 0;
	fTotalNonDSBCluster = 0;

	fDoubleCountsDI = 0;
	fDoubleCountsII = 0;
}


//--------------------------------------------------------------------------------------------------
// Record bp indices of one type of simple DNA damage (SSB or BD) in a single strand to a 1D vector.
// This method deletes the content in the provided map.
//--------------------------------------------------------------------------------------------------
std::vector<G4int> ScoreClusteredDNADamage::RecordSimpleDamage(G4double ThreshEDep,
	std::map<G4int,G4double> mapEDep)
{
	std::vector<G4int> indicesDamage;

	G4int indexBP;
	G4double eDep;

	// Loop through map of energy depositions. Add indices of energy depositions over the appropriate
	// threshold to a vector, which is returned. The entries in the map are deleted as they are
	// processed
	while ( !mapEDep.empty() )
	{
		indexBP = mapEDep.begin()->first;
		eDep = mapEDep.begin()->second;

		if (eDep >= ThreshEDep) {
			indicesDamage.push_back(indexBP);
		}
		mapEDep.erase(mapEDep.begin());
	}

	return indicesDamage;
}


//--------------------------------------------------------------------------------------------------
// Record indices of DSBs in a 1D vector. Size should always be even, corresponding to two damage
// sites per DSB. The lowest bp index is recorded first, regardless of whether in strand 1 or 2.
//--------------------------------------------------------------------------------------------------
std::vector<G4int> ScoreClusteredDNADamage::RecordDSB()
{
	std::vector<G4int> indicesDSB1D;

	std::vector<G4int>::iterator site1 = fIndicesSSB1.begin();
	std::vector<G4int>::iterator site2 = fIndicesSSB2.begin();

	// Proceed until have completely processed SSBs in either strand
	while (site1 != fIndicesSSB1.end() && site2 != fIndicesSSB2.end()) {
		G4int siteDiff = *site2 - *site1; // separation in number of bp

		// Damage in site 2 is within range of site 1 to count as DSB (either before or after)
		if (abs(siteDiff) <= fThresDistForDSB){
			// Damage in site 1 is earlier or parallel to damage in site 2
			if (*site1 <= *site2) {
				indicesDSB1D.push_back(*site1);
				indicesDSB1D.push_back(*site2);
			}
			// Damage in site 2 is earler to damage in site 1
			else {
				indicesDSB1D.push_back(*site2);
				indicesDSB1D.push_back(*site1);
			}
			site1 = fIndicesSSB1.erase(site1);
			site2 = fIndicesSSB2.erase(site2);
		}
		// Damage in site 2 is earlier than site 1 and outside range to be considered DSB
		else if (siteDiff < 0) {
			site2++;
		}
		// Damage in site 1 is earlier than site 2 and outside range to be considered DSB
		else { // if siteDiff > 0
			site1++;
		}
	}

	return indicesDSB1D;
}


//--------------------------------------------------------------------------------------------------
// Process a single sequential vector of damage indices (labelled according to damage types) to
// enable recording of two types of clustered DNA damage: Complex DSB and Non-DSB Clusters.
// Definitions are equivalent, except the former contains one or more DSB. Clustering is performed
// by calculating distances (in bp) between subsequent damage sites and comparing with a maximum
// clustering distance.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::RecordClusteredDamage()
{
	// Only 1 or 0 damages, so no clustering.
	if (fIndicesSimple.size() < 2){
		return;
	}

	 // start iterating at second index (clusters contain > 1 damage)
	std::vector<std::array<G4int,2>>::iterator site = fIndicesSimple.begin()+1;

	DamageCluster cluster;

	G4bool newCluster = true; // flag to indicate a new cluster is formed
	G4bool buildingCluster = false; // flag to indicate building on an existing cluster

	// Loop over all damages arranged in order along the strand (SSBs, BDs, and DSBs)
	while (site != fIndicesSimple.end()) {
		std::array<G4int,2> sitePrev = *(site-1); // previous damage site
		std::array<G4int,2> siteCur = *site; // current damage site

		// If damage sites are close enough to form a cluster
		if ((siteCur[0]-sitePrev[0]) <= fThresDistForCluster) {
			// A new cluster is being formed, so the "previous" damage must be added to the start
			// of the cluster
			if (newCluster) {
				AddDamageToCluster(cluster,sitePrev[0],sitePrev[1],newCluster);
				newCluster = false;
				buildingCluster = true;
			}

			// Add the current damage to the cluster
			AddDamageToCluster(cluster,siteCur[0],siteCur[1],newCluster);
		}
		// Damage sites are too far away to form a cluster
		else {
			// Previous site was part of a cluster, but now cluster has ended so record it
			if (buildingCluster) {
				buildingCluster = false;
				newCluster = true;
				RecordCluster(cluster);
			}
		}
		site++;
	}
	// Handle case if was building a cluster when reached end of list
	if (buildingCluster) {
		buildingCluster = false;
		RecordCluster(cluster);
	}
}


//--------------------------------------------------------------------------------------------------
// Add a new DNA damage site to a cluster.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::AddDamageToCluster(DamageCluster& cluster, G4int damageSite,
												G4int damageType, G4bool newCluster) {
	// If this cluster is a new cluster, set the new damage as the start site. Otherwise, set the
	// new damage as the end site.
	if (newCluster) {
		cluster.start = damageSite;
	}
	else {
		cluster.end = damageSite;
	}

	// Increment the correct damage counter for this cluster, and correspondingly decrement the
	// correct individual damage counter
	if (damageType == fIdSSB) {
		cluster.numSSB++;
		fTotalSSB--;
	}
	else if (damageType == fIdBD) {
		cluster.numBD++;
		fTotalBD--;
	}
	else if (damageType == fIdDSB) {
		cluster.numDSB++;
		fTotalDSB--;
	}
	// Throw an error if an unrecognized damage type is added
	else {
		G4cout << "An error has arisen while scoring clustered DNA damage." << G4endl;
		G4cout << "The following integer damage label is unrecognized:" << G4endl;
		G4cout << damageType << G4endl;
		exit(0);
	}
}


//--------------------------------------------------------------------------------------------------
// Add the details of a cluster to the appropriate member variables (distinguishing a Complex DSB
// from a Non-DSB Cluster). Update counts of appropriate type of cluster. Reset the cluster
// variable.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::RecordCluster(DamageCluster& cluster) {
	// Handle Simple DSB (i.e. not a cluster, so record nothing & reset)
	if (cluster.numDSB == 2 && cluster.numSSB == 0 && cluster.numBD == 0) {
		fTotalDSB += 2;
		cluster = DamageCluster();
	}
	// Handle Complex DSB
	else if (cluster.numDSB > 0) {
		cluster.numDSB = cluster.numDSB/2;

		fComplexDSBSizes.push_back(cluster.end - cluster.start + 1);
		fComplexDSBNumSSB.push_back(cluster.numSSB);
		fComplexDSBNumBD.push_back(cluster.numBD);
		fComplexDSBNumDSB.push_back(cluster.numDSB);
		fComplexDSBNumDamage.push_back(cluster.numSSB + cluster.numBD + cluster.numDSB);

		fTotalComplexDSB++;

		cluster = DamageCluster();
	}
	// Handle Non-DSB cluster
	else {
		fNonDSBClusterSizes.push_back(cluster.end - cluster.start + 1);
		fNonDSBClusterNumSSB.push_back(cluster.numSSB);
		fNonDSBClusterNumBD.push_back(cluster.numBD);
		fNonDSBClusterNumDamage.push_back(cluster.numSSB + cluster.numBD);

		fTotalNonDSBCluster++;

		cluster = DamageCluster();
	}
}


//--------------------------------------------------------------------------------------------------
// Combine class member vectors containing various types of damages into a single, ordered, vector
// of all damges in a DNA fibre (both strands). Each vector element is a 2-element array containing
// (i) the bp index of the damage and (ii) an integer indicating the type of damage.
//--------------------------------------------------------------------------------------------------
std::vector<std::array<G4int,2>> ScoreClusteredDNADamage::CombineSimpleDamage() {
	// G4cout << "COMBINING SIMPLE DAMAGE" << G4endl;
	std::vector<std::array<G4int,2>> indicesSimple;

	// Fake data for testing
	// fIndicesSSB1 = {3,4,5,7,9};
	// fIndicesBD1 = {1,4,7};
	// fIndicesSSB2 = {1,12};
	// fIndicesBD2 = {6,7,10,12};
	// fIndicesDSB = {3,5,11,14};

	// SSBs in strand 1
	std::vector<G4int>::iterator itSSB1 = fIndicesSSB1.begin();
	while (itSSB1 != fIndicesSSB1.end()) {
		indicesSimple.push_back({*itSSB1,fIdSSB});
		itSSB1++;
	}

	// BDs in strand 1
	std::vector<G4int>::iterator itBD1 = fIndicesBD1.begin();
	std::vector<std::array<G4int,2>>::iterator itSimple = indicesSimple.begin();
	while (itBD1 != fIndicesBD1.end() && itSimple != indicesSimple.end()) {
		if (*itBD1 < (*itSimple)[0]) {
			itSimple = indicesSimple.insert(itSimple,{*itBD1,fIdBD});
			itBD1++;
		}
		else {
			itSimple++;
		}
	}
	while (itBD1 != fIndicesBD1.end()) {
		indicesSimple.push_back({*itBD1,fIdBD});
		itBD1++;
	}

	// SSBs in strand 2
	std::vector<G4int>::iterator itSSB2 = fIndicesSSB2.begin();
	itSimple = indicesSimple.begin();
	while (itSSB2 != fIndicesSSB2.end() && itSimple != indicesSimple.end()) {
		if (*itSSB2 < (*itSimple)[0]) {
			itSimple = indicesSimple.insert(itSimple,{*itSSB2,fIdSSB});
			itSSB2++;
		}
		else {
			itSimple++;
		}
	}
	while (itSSB2 != fIndicesSSB2.end()) {
		indicesSimple.push_back({*itSSB2,fIdSSB});
		itSSB2++;
	}

	// BDs in strand 2
	std::vector<G4int>::iterator itBD2 = fIndicesBD2.begin();
	itSimple = indicesSimple.begin();
	while (itBD2 != fIndicesBD2.end() && itSimple != indicesSimple.end()) {
		if (*itBD2 < (*itSimple)[0]) {
			itSimple = indicesSimple.insert(itSimple,{*itBD2,fIdBD});
			itBD2++;
		}
		else {
			itSimple++;
		}
	}
	while (itBD2 != fIndicesBD2.end()) {
		indicesSimple.push_back({*itBD2,fIdBD});
		itBD2++;
	}

	// DSBs
	std::vector<G4int>::iterator itDSB = fIndicesDSB.begin();
	itSimple = indicesSimple.begin();
	while (itDSB != fIndicesDSB.end() && itSimple != indicesSimple.end()) {
		if (*itDSB < (*itSimple)[0]) {
			itSimple = indicesSimple.insert(itSimple,{*itDSB,fIdDSB});
			itDSB++;
		}
		else {
			itSimple++;
		}
	}
	while (itDSB != fIndicesDSB.end()) {
		indicesSimple.push_back({*itDSB,fIdDSB});
		itDSB++;
	}

	// If want to output full list of simple damages:
	// G4cout << "------------------------------------------" << G4endl;
	// for (int i = 0; i < indicesSimple.size(); i++) {
	// 	G4cout << "site = " << indicesSimple[i][0] << ", type = " << indicesSimple[i][1] << G4endl;
	// }
	// G4cout << "------------------------------------------" << G4endl;

	return indicesSimple;
}


//--------------------------------------------------------------------------------------------------
// Calculate the order of magnitude (base 10) of a positive integer value.
//--------------------------------------------------------------------------------------------------
G4int ScoreClusteredDNADamage::CalculateIntegerMagnitude(G4int value) {
	G4int orderOfMagnitude = 0;
	G4double base = 10.0;
	G4double valueD = value*1.0;

	while (valueD >= base) {
		valueD = valueD/base;
		orderOfMagnitude++;
	}

	return static_cast<int>(pow(10,orderOfMagnitude));
}


//--------------------------------------------------------------------------------------------------
// Print DNA damage inforamtion for the current event (as per fEventID) to console
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::PrintDNADamageToConsole() {
	G4cout << "##########################################################" << G4endl;
	G4cout << "Event #" << fEventID << G4endl;
	G4cout << "# of SSBs = " << fTotalSSB << G4endl;
	G4cout << "# of BDs = " << fTotalBD << G4endl;
	G4cout << "# of simple DSBs = " << fTotalDSB << G4endl;
	G4cout << "# of complex DSBs = " << fTotalComplexDSB << G4endl;
	for (G4int i=0; i < fTotalComplexDSB; i++) {
		G4cout << "-----------------------------" << G4endl;
		G4cout << "\tSize = " << fComplexDSBSizes[i] << G4endl;
		G4cout << "\tNum SSB = " << fComplexDSBNumSSB[i] << G4endl;
		G4cout << "\tNum BD = " << fComplexDSBNumBD[i] << G4endl;
		G4cout << "\tNum DSB = " << fComplexDSBNumDSB[i] << G4endl;
		G4cout << "\tNum damage = " << fComplexDSBNumDamage[i] << G4endl;
	}
	G4cout << "# of non-DSB clusters = " << fTotalNonDSBCluster << G4endl;
	for (G4int i=0; i < fTotalNonDSBCluster; i++) {
		G4cout << "-----------------------------" << G4endl;
		G4cout << "\tSize = " << fNonDSBClusterSizes[i] << G4endl;
		G4cout << "\tNum SSB = " << fNonDSBClusterNumSSB[i] << G4endl;
		G4cout << "\tNum BD = " << fNonDSBClusterNumBD[i] << G4endl;
		G4cout << "\tNum damage = " << fNonDSBClusterNumDamage[i] << G4endl;
	}
	G4cout << "##########################################################" << G4endl;
}


//--------------------------------------------------------------------------------------------------
// Create a fake map of energy depositions to validate scoring & clustering behaviour.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::CreateFakeEnergyMap() {
	// Test cases to validate damage clustering algorithm is functioning properly
	G4double eng = 20./eV;
	fThresDistForCluster = 4;
	fThresDistForDSB = 4;

	// Test case 1
	std::map<G4int,G4double> back1 = {{1,eng},{5,eng},{10,eng}};
	std::map<G4int,G4double> base1 = {{10,eng},{12,eng}};
	std::map<G4int,G4double> base2 = {{11,eng},{17,eng}};
	std::map<G4int,G4double> back2 = {{1,eng},{3,eng}};

	// Test case 2
	// std::map<G4int,G4double> back1 = {{6,eng}};
	// std::map<G4int,G4double> base1 = {{15,eng},{20,eng}};
	// std::map<G4int,G4double> base2 = {{15,eng}};
	// std::map<G4int,G4double> back2 = {{2,eng},{10,eng}};


	fMapEdepStrand1Backbone[0][0] = back1;
	fMapEdepStrand1Base[0][0] = base1;
	fMapEdepStrand2Base[0][0] = base2;
	fMapEdepStrand2Backbone[0][0] = back2;
	// fFiberMapEdepStrand1Backbone[0] = back1;
	// fFiberMapEdepStrand1Base[0] = base1;
	// fFiberMapEdepStrand2Base[0] = base2;
	// fFiberMapEdepStrand2Backbone[0] = back2;
}

void ScoreClusteredDNADamage::Print1DVectorContents(std::vector<G4int> pVector) {
	G4cout << "\t";
	for (G4int element : pVector) {
		G4cout << element << ' ';
	}
	G4cout << G4endl;
}
