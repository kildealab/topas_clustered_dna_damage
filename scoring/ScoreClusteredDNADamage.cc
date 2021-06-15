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
	DamageCluster() : numSSB(0), numSSB_direct(0), numSSB_indirect(0),
										numBD(0), numBD_direct(0), numBD_indirect(0),
										numDSB(0), numDSB_hybrid(0), numDSB_direct(0), numDSB_indirect(0),
										start(0), end(0), size(0) {}

	G4int numSSB; // # of SSB in cluster
	G4int numSSB_direct;
	G4int numSSB_indirect;

	G4int numBD; // # of base damages in cluster
	G4int numBD_direct;
	G4int numBD_indirect;

	G4int numDSB; // # of DSB in cluster
	G4int numDSB_hybrid;
	G4int numDSB_direct;
	G4int numDSB_indirect;

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
	fTotalSSB_direct = 0;
	fTotalSSB_indirect = 0;

	fTotalBD = 0;
	fTotalBD_direct = 0;
	fTotalBD_indirect = 0;

	fTotalDSB = 0;
	fTotalDSB_direct = 0;
	fTotalDSB_indirect = 0;
	fTotalDSB_hybrid = 0;

	fTotalComplexDSB = 0;
	fTotalComplexDSB_direct = 0;
	fTotalComplexDSB_indirect = 0;
	fTotalComplexDSB_hybrid = 0;

	fTotalNonDSBCluster = 0;
	fTotalNonDSBCluster_direct = 0;
	fTotalNonDSBCluster_indirect = 0;
	fTotalNonDSBCluster_hybrid = 0;

	fDoubleCountsDD = 0;
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
	fNtuple->RegisterColumnI(&fTotalSSB, "Total single strand breaks"); // Number of SSB caused by this primary particle
	if (fIncludeDirectDamage)
		fNtuple->RegisterColumnI(&fTotalSSB_direct, "SSBs direct");
	if (fIncludeIndirectDamage)
		fNtuple->RegisterColumnI(&fTotalSSB_indirect, "SSBs indirect");

	fNtuple->RegisterColumnI(&fTotalDSB, "Total double strand breaks"); // Number of simple DSB caused by this primary particle
	if (fIncludeDirectDamage)
		fNtuple->RegisterColumnI(&fTotalDSB_direct, "DSBs direct");
	if (fIncludeIndirectDamage)
		fNtuple->RegisterColumnI(&fTotalDSB_indirect, "DSBs indirect");
	if (fIncludeDirectDamage && fIncludeIndirectDamage)
		fNtuple->RegisterColumnI(&fTotalDSB_hybrid, "DSBs hybrid");

	fNtuple->RegisterColumnI(&fTotalBD, "Total base damages"); // Number of BD caused by this primary particle
	if (fIncludeDirectDamage)
		fNtuple->RegisterColumnI(&fTotalBD_direct, "BDs direct");
	if (fIncludeIndirectDamage)
		fNtuple->RegisterColumnI(&fTotalBD_indirect, "BDs indirect");

	if (fScoreClusters) {
		fNtuple->RegisterColumnI(&fTotalComplexDSB, "Complex DSBs"); // Number of Complex DSB caused by this primary particle
		if (fIncludeDirectDamage)
			fNtuple->RegisterColumnI(&fTotalComplexDSB_direct, "Complex DSBs direct");
		if (fIncludeIndirectDamage)
			fNtuple->RegisterColumnI(&fTotalComplexDSB_indirect, "Complex DSBs indirect");
		if (fIncludeDirectDamage && fIncludeIndirectDamage)
			fNtuple->RegisterColumnI(&fTotalComplexDSB_hybrid, "Complex DSBs hybrid");

		fNtuple->RegisterColumnI(&fTotalNonDSBCluster, "Non-DSB clusters"); // Number of Non-DSB Clusters caused by this primary particle
		if (fIncludeDirectDamage)
			fNtuple->RegisterColumnI(&fTotalNonDSBCluster_direct, "Non-DSB clusters direct");
		if (fIncludeIndirectDamage)
			fNtuple->RegisterColumnI(&fTotalNonDSBCluster_indirect, "Non-DSB clusters indirect");
		if (fIncludeDirectDamage && fIncludeIndirectDamage)
			fNtuple->RegisterColumnI(&fTotalNonDSBCluster_hybrid, "Non-DSB clusters hybrid");
	}

	if (fIncludeDirectDamage)
		fNtuple->RegisterColumnI(&fDoubleCountsDD, "Double counts direct-direct");
	if (fIncludeIndirectDamage)
		fNtuple->RegisterColumnI(&fDoubleCountsII, "Double counts indirect-indirect");
	if (fIncludeDirectDamage && fIncludeIndirectDamage)
		fNtuple->RegisterColumnI(&fDoubleCountsDI, "Double counts direct-indirect");
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
	// Chemistry-related parameters
	//----------------------------------------------------------------------------------------------
	G4String* modules = fPm->GetStringVector(GetFullParmName("Modules"));
	G4int numModules = fPm->GetVectorLength(GetFullParmName("Modules"));
	fHasChemistryModule = false;
	for (int i = 0; i < numModules; i++) {
		modules[i].toLower();
		if (modules[i].contains("chemistry")) {
			fHasChemistryModule = true;
			G4cout << " Chemistry module found: " << modules[i] << G4endl;
		}
	}

	if (fHasChemistryModule) {
		G4ConfigurationIterator mol_iterator = G4MoleculeTable::Instance()->GetConfigurationIterator();
		// G4MoleculeDefinitionIterator mol_iterator =	G4MoleculeTable::Instance()->GetDefintionIterator();
		while ((mol_iterator)()) {
			G4String mol_name = mol_iterator.value()->GetUserID();
			G4int mol_ID =  mol_iterator.value()->GetMoleculeID();
			G4String paramNameSSB = "DamageProbabilityOnInteraction/SSB/" + mol_name;
			G4String paramNameBD = "DamageProbabilityOnInteraction/BD/" + mol_name;

			G4cout << " Registering molecule: " << mol_name << "\tID: " << mol_ID << G4endl;

			if (mol_name == "OH") {
				if (fPm->ParameterExists(GetFullParmName(paramNameSSB)))
					fMoleculeDamageProb_SSB.emplace(mol_ID, fPm->GetUnitlessParameter(GetFullParmName(paramNameSSB)));
				else
					fMoleculeDamageProb_SSB.emplace(mol_ID, 0.4);
			}
			else {
				if (fPm->ParameterExists(GetFullParmName(paramNameSSB)))
					fMoleculeDamageProb_SSB.emplace(mol_ID, fPm->GetUnitlessParameter(GetFullParmName(paramNameSSB)));
				else
					fMoleculeDamageProb_SSB.emplace(mol_ID, 0.0);
			}

			if (fPm->ParameterExists(GetFullParmName(paramNameBD)))
				fMoleculeDamageProb_BD.emplace(mol_ID, fPm->GetUnitlessParameter(GetFullParmName(paramNameBD)));
			else
				fMoleculeDamageProb_BD.emplace(mol_ID, 0.0);
		}

		if (fPm->ParameterExists(GetFullParmName("SpeciesToKillByDNAVolumes"))) {
			G4String* speciesToKillbyDNAVolumes_names = fPm->GetStringVector(GetFullParmName("SpeciesToKillByDNAVolumes"));
			G4int vectorLength = fPm->GetVectorLength(GetFullParmName("SpeciesToKillByDNAVolumes"));
			for (int i = 0; i < vectorLength; i++){
				G4String mol_name = speciesToKillbyDNAVolumes_names[i];
				fSpeciesToKillByDNAVolumes.push_back(G4MoleculeTable::Instance()->GetConfiguration(mol_name)->GetMoleculeID());
				G4cout << " To be killed by DNA volume: " << mol_name << G4endl;
			}
		}
		else {
			fSpeciesToKillByDNAVolumes.push_back(G4MoleculeTable::Instance()->GetConfiguration("OH")->GetMoleculeID());
		}

		if (fPm->ParameterExists(GetFullParmName("SpeciesToKillByHistones"))) {
			G4String* speciestoKillbyHistones_names = fPm->GetStringVector(GetFullParmName("SpeciesToKillByHistones"));
			G4int vectorLength = fPm->GetVectorLength(GetFullParmName("SpeciesToKillByHistones"));
			for (int i = 0; i < vectorLength; i++){
				G4String mol_name = speciestoKillbyHistones_names[i];
				fspeciesToKillByHistones.push_back(G4MoleculeTable::Instance()->GetConfiguration(mol_name)->GetMoleculeID());
				G4cout << " To be killed by histone volume: " << mol_name << G4endl;
			}
		}
		else {
			fspeciesToKillByHistones.push_back(G4MoleculeTable::Instance()->GetConfiguration("OH")->GetMoleculeID());
			fspeciesToKillByHistones.push_back(G4MoleculeTable::Instance()->GetConfiguration("e_aq")->GetMoleculeID());
			fspeciesToKillByHistones.push_back(G4MoleculeTable::Instance()->GetConfiguration("H")->GetMoleculeID());
		}

		if (fPm->ParameterExists(GetFullParmName("HistonesAsScavenger")))
			fHistonesAsScavenger = fPm->GetBooleanParameter(GetFullParmName("HistonesAsScavenger"));
		else
			fHistonesAsScavenger = true;
	}

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
		fIncludeIndirectDamage = false;
	if (fIncludeIndirectDamage && !fHasChemistryModule) {
		G4cerr << "Error: Indirect damage scoring cannot be enabled if there is no chemistry module." << G4endl;
		exit(0);
	}

	G4cout << "StepStatuses: " << G4endl;
	G4cout << "\tfWorldBoundary: " << fWorldBoundary  << G4endl;
	G4cout << "\tfGeomBoundary: " << fGeomBoundary   << G4endl;
	G4cout << "\tfAtRestDoItProc: " << fAtRestDoItProc   << G4endl;
	G4cout << "\tfAlongStepDoItProc: " << fAlongStepDoItProc   << G4endl;
	G4cout << "\tfPostStepDoItProc: " << fPostStepDoItProc   << G4endl;
	G4cout << "\tfUserDefinedLimit: " << fUserDefinedLimit   << G4endl;
	G4cout << "\tfExclusivelyForcedProc: " << fExclusivelyForcedProc   << G4endl;
	G4cout << "\tfUndefined: " << fUndefined   << G4endl;

	//----------------------------------------------------------------------------------------------
	// Material of DNA residue volumes in which to score
	//----------------------------------------------------------------------------------------------
	G4String DNAMaterialName = "G4_WATER";
	if ( fPm->ParameterExists(GetFullParmName("DNAMaterialName")))
		DNAMaterialName = fPm->GetStringParameter(GetFullParmName("DNAMaterialName"));
	fDNAMaterial = GetMaterial(DNAMaterialName);

	G4String HistoneMaterialName = "G4_WATER";
	if ( fPm->ParameterExists(GetFullParmName("HistoneMaterialName")))
		HistoneMaterialName = fPm->GetStringParameter(GetFullParmName("HistoneMaterialName"));
	fHistoneMaterial = GetMaterial(HistoneMaterialName);

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
// deposition took place. This is faster than using string comparisons. This method is only called
// for energy depositions in the sensitive volumes (i.e. residues) by using material filtering, as
// defined in the parameter file with "OnlyIncludeIfInMaterial" parameter.
//--------------------------------------------------------------------------------------------------
G4bool ScoreClusteredDNADamage::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	fNumProcessHitsCalls++;
	G4double edep = aStep->GetTotalEnergyDeposit(); // In eV;
	fTotalEdep += edep; // running sum of energy deposition in entire volume

	// Material filtering (only proceed if in sensitive DNA volumes)
	G4Material* material = aStep->GetPreStepPoint()->GetMaterial();
	G4bool isDNAMaterial = (material == fDNAMaterial);
	G4bool isHistoneMaterial = (material == fHistoneMaterial);

	if (!isDNAMaterial && !isHistoneMaterial) {
		// G4cout << "\tNot DNA nor histone material: " << aStep->GetPreStepPoint()->GetMaterial()->GetName() << G4endl;
		return false;
	}

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
	if (fNumFibers > 1) {
		fFiberID = touchable->GetCopyNumber(fParentIndexFiber);
	}
	// Determine the indices defining the volume in which hit occured by parsing the copy number of the Physical Volume.
	G4int volID = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
	G4int strandID = volID / 1000000;
	G4int residueID = (volID - (strandID*1000000)) / 100000;
	G4int bpID = volID - (strandID*1000000) - (residueID*100000);

	// Determines whether track is physical or chemical
	G4int trackID = aStep->GetTrack()->GetTrackID();

	//----------------------------------------------------------------------------------------------
	// If this hit deposits energy (in sensitive DNA volume), update the appropriate energy deposition
	// map
	//----------------------------------------------------------------------------------------------
	if (fIncludeDirectDamage && edep > 0 && trackID >= 0 && isDNAMaterial) { // energy deposition should be from physical tracks
		// if (aStep->IsFirstStepInVolume()) {
		// 	G4cout << "First step in volume by physical track: " << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << G4endl;
		// }
		//------------------------------------------------------------------------------------------
		// Use the DNA strand ID, residue ID, and nucleotide ID to increment the energy deposited
		// in the appropriate energy deposition map. Maps are indexed as follows:
		// First index specifies the voxel
		// Second index specifies DNA fibre
		// Third index specifies the bp index
		//------------------------------------------------------------------------------------------
		if ( strandID == 0 ){ // first strand
			if (residueID == fVolIdPhosphate || residueID == fVolIdDeoxyribose){
				fMapEdepStrand1Backbone[fVoxelID][fFiberID][bpID] += edep;
			}
			else if (residueID == fVolIdBase) {
				fMapEdepStrand1Base[fVoxelID][fFiberID][bpID] += edep;
			}
		}
		else if ( strandID == 1 ) { // second strand
			if (residueID == fVolIdPhosphate || residueID == fVolIdDeoxyribose){
				fMapEdepStrand2Backbone[fVoxelID][fFiberID][bpID] += edep;
			}
			else if (residueID == fVolIdBase) {
				fMapEdepStrand2Base[fVoxelID][fFiberID][bpID] += edep;
			}
		}
		else {
			G4cerr << "Error: (While scoring direct damage) The following strandID is unrecognized: " << strandID << G4endl;
			exit(0);
		}
		return true;
	}

	// Indirect damage
	if (fIncludeIndirectDamage && trackID < 0) { // chemical tracks
		// Get molecule info
		G4int moleculeID = GetMolecule(aStep->GetTrack())->GetMoleculeID();
		G4String moleculeName = GetMolecule(aStep->GetTrack())->GetName();
		G4String postStepProcessName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
		G4String volumeName = touchable->GetVolume()->GetName();

		if (aStep->IsFirstStepInVolume()) {
			G4cout << "First step in volume by chemical track: " << moleculeName << G4endl; // never called :(
		}

		// Kill species generated inside DNA volumes and histones
		G4bool justGotGenerated	= (aStep->GetTrack()->GetCurrentStepNumber() == 1);
		if (justGotGenerated) {
			G4cout << "Killing " << moleculeName << " (" << trackID << ") upon generation in " << volumeName << G4endl;
			PrintStepInfo(aStep); //debugging
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return false;
		}

		// Determine if molecule just entered a volume
		G4bool justEnteredVolume	= (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary);

		// Determine if damage is inflicted
		G4bool isDamaged = IsDamageInflicted(moleculeID, residueID);
		if (isDNAMaterial && justEnteredVolume && isDamaged) {
			// Check which damage map to update
			if ( strandID == 0 ) { // first strand
				if (residueID == fVolIdPhosphate || residueID == fVolIdDeoxyribose){ // backbone damage
					fIndices = &fMapIndDamageStrand1Backbone[fVoxelID][fFiberID];
				}
				else if (residueID == fVolIdBase){ // base damage
					fIndices = &fMapIndDamageStrand1Base[fVoxelID][fFiberID];
				}
			}
			else { // second strand
				if (residueID == fVolIdPhosphate || residueID == fVolIdDeoxyribose){ // backbone damage
					fIndices = &fMapIndDamageStrand2Backbone[fVoxelID][fFiberID];
				}
				else if (residueID == fVolIdBase){ // base damage
					fIndices = &fMapIndDamageStrand2Base[fVoxelID][fFiberID];
				}
			}

			// Check if backbone or base has already been damaged previously via indirect action
			if (IsElementInVector(bpID, *fIndices)){
				G4cout << "\tBase pair already damaged. voxel: " << fVoxelID << " fiber: " << fFiberID << " bp: " << bpID << G4endl;
				return false;
			}

			// Record damaged nucleotide
			fIndices->push_back(bpID);
			G4cout << "\tDAMAGE INDUCED! voxel: " << fVoxelID << " fiber: " << fFiberID << " bp: " << bpID << G4endl;
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return true;
		}

		// Kill certain species interacting with DNA volumes
		G4bool isSpeciesToKillByDNA = IsElementInVector(moleculeID, fSpeciesToKillByDNAVolumes);
		if (isDNAMaterial && justEnteredVolume && isSpeciesToKillByDNA) {
			G4cout << "No damage. Track killed. moleculeName: " << moleculeName << " (" << trackID << ")" << G4endl;
			PrintStepInfo(aStep); //debugging
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return false;
		}

		// Kill certain species diffusing in histone volumes
		if (fHistonesAsScavenger && isHistoneMaterial) {
			G4bool isSpeciesToKillByHistone = IsElementInVector(moleculeID, fspeciesToKillByHistones);
			if (isSpeciesToKillByHistone && !justEnteredVolume) {
				G4cout << "Track killed because diffusing in histone. moleculeName:" << moleculeName << " (" << trackID << ")" << G4endl;
				PrintStepInfo(aStep); //debugging
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				return false;
			}
		}
	} // negative trackID

	return false;
}


//--------------------------------------------------------------------------------------------------
// This helper method checks whether an element is in a vector.
//--------------------------------------------------------------------------------------------------
G4bool ScoreClusteredDNADamage::IsElementInVector(G4int pElement, std::vector<G4int> &pVector) {
	if (std::find(pVector.begin(), pVector.end(), pElement) != pVector.end())
		return true;
	else
		return false;
}

void ScoreClusteredDNADamage::RemoveElementFromVector(G4int pElement, std::vector<G4int> &pVector) {
	if (!IsElementInVector(pElement, pVector)) {
		G4cout << "Element (" << pElement << ") not in vector" << G4endl;
		Print1DVectorContents(pVector);
		return;
	}
	pVector.erase(std::remove(pVector.begin(), pVector.end(), pElement), pVector.end());
}


//--------------------------------------------------------------------------------------------------
// This helper method checks whether an indirect damage has been induced given the moleculeID
//--------------------------------------------------------------------------------------------------
G4bool ScoreClusteredDNADamage::IsDamageInflicted(G4int pMoleculeID, G4int pDNAVolumeID) {
	G4float prob_damage = 0.;

	if ( fMoleculeDamageProb_SSB.find(pMoleculeID) == fMoleculeDamageProb_SSB.end() && fMoleculeDamageProb_BD.find(pMoleculeID) == fMoleculeDamageProb_BD.end()) {
  	// pMoleculeID not found in fMoleculeDamageProb_SSB and fMoleculeDamageProb_BD
		G4cerr << "\tmoleculeID NOT LISTED: " << pMoleculeID << G4endl;
		return false;
	}

	if (pDNAVolumeID == fVolIdPhosphate || pDNAVolumeID == fVolIdDeoxyribose)
		prob_damage = fMoleculeDamageProb_SSB[pMoleculeID]; // returns 0 if pMoleculeID is not a "key" in the map
	else if (pDNAVolumeID == fVolIdBase)
		prob_damage = fMoleculeDamageProb_BD[pMoleculeID]; // returns 0 if pMoleculeID is not a "key" in the map

	if (prob_damage == 0.)
		return false;

	// Generate random float between 0 and 1
	G4float prob_random = G4UniformRand();

	if (prob_random <= prob_damage)
		return true;
	else
		return false;
}

void ScoreClusteredDNADamage::PrintStepInfo(G4Step* aStep) {
	G4cout << "\tpreStep volume: " << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
	G4cout << "\tpostStep volume: " << aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
	G4cout << "\tpreStep status: " << aStep->GetPreStepPoint()->GetStepStatus() << G4endl;
	G4cout << "\tpostStep status: " << aStep->GetPostStepPoint()->GetStepStatus() << G4endl;
	if (aStep->GetPreStepPoint()->GetStepStatus() != fUndefined)
		G4cout << "\tPreStep process name: " << aStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
	else
		G4cout << "\tPreStep process name: Undefined" << G4endl;
	G4cout << "\tPostStep processname: " << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
	G4cout << "\tCurrentStepNumber of track: " << aStep->GetTrack()->GetCurrentStepNumber() << G4endl;
	G4cout << "\tParent trackID: " << aStep->GetTrack()->GetParentID() << G4endl;
	if (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName())
		G4cout << "\tPreStep volume is NOT the same as PostStep volume." << G4endl;
	else
		G4cout << "\tPreStep volume is the same as PostStep volume." << G4endl;
	if (aStep->GetTrack()->GetCurrentStepNumber() == 1)
		G4cout << "\tCurrentStepNumber of track is 1." << G4endl;
	if (aStep->GetTrack()->GetParentID() >= 0)
		G4cout << "\tParent track is physical." << G4endl;
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
		outHeader << "# of direct SSB" << fDelimiter;
		outHeader << "# of indirect SSB" << fDelimiter;
		outHeader << "# of BD" << fDelimiter;
		outHeader << "# of direct BD" << fDelimiter;
		outHeader << "# of indirect BD" << fDelimiter;
		outHeader << "# of DSB" << fDelimiter;
		outHeader << "# of hybrid DSB" << fDelimiter;
		outHeader << "# of direct DSB" << fDelimiter;
		outHeader << "# of indirect DSB" << G4endl;
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
		outFile << fComplexDSBNumSSB_direct[i] << fDelimiter;
		outFile << fComplexDSBNumSSB_indirect[i] << fDelimiter;
		outFile << fComplexDSBNumBD[i] << fDelimiter;
		outFile << fComplexDSBNumBD_direct[i] << fDelimiter;
		outFile << fComplexDSBNumBD_indirect[i] << fDelimiter;
		outFile << fComplexDSBNumDSB[i] << fDelimiter;
		outFile << fComplexDSBNumDSB_hybrid[i] << fDelimiter;
		outFile << fComplexDSBNumDSB_direct[i] << fDelimiter;
		outFile << fComplexDSBNumDSB_indirect[i] << G4endl;
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
		outHeader << "# of direct SSB" << fDelimiter;
		outHeader << "# of indirect SSB" << fDelimiter;
		outHeader << "# of BD" << fDelimiter;
		outHeader << "# of direct BD" << fDelimiter;
		outHeader << "# of indirect BD" << G4endl;
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
		outFile << fNonDSBClusterNumSSB_direct[i] << fDelimiter;
		outFile << fNonDSBClusterNumSSB_indirect[i] << fDelimiter;
		outFile << fNonDSBClusterNumBD[i] << fDelimiter;
		outFile << fNonDSBClusterNumBD_direct[i] << fDelimiter;
		outFile << fNonDSBClusterNumBD_indirect[i] << G4endl;
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

			// Process SSBs in both strands to determine if there are any DSB
			// fIndicesDSB = RecordDSB();
			G4int totalFiberSSB_direct = 0, totalFiberSSB_indirect = 0, totalFiberBD_direct = 0, totalFiberBD_indirect = 0;
			G4int totalFiberDSB_direct = 0, totalFiberDSB_indirect = 0, totalFiberDSB_hybrid = 0;

			// G4bool fIncludeHybridDamage = true;
			if (fIncludeDirectDamage && fIncludeIndirectDamage) { // && fIncludeHybridDamage) {
				fIndicesDSB_hybrid = RecordDSB(fIdHybrid);
				totalFiberDSB_hybrid = fIndicesDSB_hybrid.size();
				fTotalDSB_hybrid += totalFiberDSB_hybrid;
				fTotalDSB += totalFiberDSB_hybrid;
			}

			if (fIncludeDirectDamage) {
				fIndicesDSB_direct = RecordDSB(fIdDirect);
				totalFiberDSB_direct = fIndicesDSB_direct.size();
				fTotalDSB_direct += totalFiberDSB_direct;
				fTotalDSB += totalFiberDSB_direct;

				totalFiberSSB_direct = fIndicesSSB1_direct.size() + fIndicesSSB2_direct.size();
				fTotalSSB_direct += totalFiberSSB_direct;
				fTotalSSB += totalFiberSSB_direct;

				totalFiberBD_direct = fIndicesBD1_direct.size() + fIndicesBD2_direct.size();
				fTotalBD_direct += totalFiberBD_direct;
				fTotalBD += totalFiberBD_direct;
			}

			if (fIncludeIndirectDamage) {
				fIndicesDSB_indirect = RecordDSB(fIdIndirect);
				totalFiberDSB_indirect = fIndicesDSB_indirect.size();
				fTotalDSB_indirect += totalFiberDSB_indirect;
				fTotalDSB += totalFiberDSB_indirect;

				totalFiberSSB_indirect = fIndicesSSB1_indirect.size() + fIndicesSSB2_indirect.size();
				fTotalSSB_indirect += totalFiberSSB_indirect;
				fTotalSSB += totalFiberSSB_indirect;

				totalFiberBD_indirect = fIndicesBD1_indirect.size() + fIndicesBD2_indirect.size();
				fTotalBD_indirect += totalFiberBD_indirect;
				fTotalBD += totalFiberBD_indirect;
			}

			// If recording clustered damage, combine all damages into a single, sequential vector
			// of damage that indicates the type and bp index. Then process this vector to determine
			// clustered damage yields
			if (fScoreClusters) {
				fIndicesSimple = CombineSimpleDamage();
				RecordClusteredDamage();
			}

			// if (totalFiberSSB_direct + totalFiberSSB_indirect > 0) {
			// 	G4cout << "\tfTotalSSB: " << fTotalSSB << " voxel: " << iVoxel << " fiber: " << iFiber << G4endl;
			// 	G4cout << "\tfTotalSSB_direct: " << fTotalSSB_direct << " voxel: " << iVoxel << " fiber: " << iFiber << G4endl;
			// }
			// if (totalFiberBD_direct + totalFiberBD_indirect > 0)
			// 	G4cout << "\tfTotalBD: " << fTotalBD << " voxel: " << iVoxel << " fiber: " << iFiber << G4endl;
			// if (totalFiberDSB_hybrid + totalFiberDSB_direct + totalFiberDSB_indirect > 0)
			// 	G4cout << "\tfTotalDSB: " << fTotalDSB/2 << " voxel: " << iVoxel << " fiber: " << iFiber << G4endl;

			// If recording damage on a fiber-by-fiber basis, fill the output ntuple
			if (fRecordDamagePerFiber) {
				fTotalDSB = fTotalDSB/2;
				fTotalDSB_hybrid = fTotalDSB_hybrid/2;
				fTotalDSB_direct = fTotalDSB_direct/2;
				fTotalDSB_indirect = fTotalDSB_indirect/2;
				fNtuple->Fill(); // Move this to outside loop if aggregating over all fibres

				// Reset variables before next fibre (not aggregating over all fibres)
				ResetDamageCounterVariables();
			}
		}
	}

	// If recording damage aggregated over all fibers, fill the output ntuple
	if (!fRecordDamagePerFiber) {
		fTotalDSB = fTotalDSB/2;
		fTotalDSB_hybrid = fTotalDSB_hybrid/2;
		fTotalDSB_direct = fTotalDSB_direct/2;
		fTotalDSB_indirect = fTotalDSB_indirect/2;
		fFiberID = fAggregateValueIndicator;
		fNtuple->Fill(); // Move this to outside loop if aggregating over all fibres
	}
	// PrintDNADamageToConsole(); // debugging;
}


//--------------------------------------------------------------------------------------------------
// This method merges and resolves duplicates of the damage yields from direct and indirect damage.
//--------------------------------------------------------------------------------------------------
std::vector<G4int> ScoreClusteredDNADamage::MergeDamageIndices(std::vector<G4int> &pDamageIndices_direct,
	std::vector<G4int> &pDamageIndices_indirect)
{
	std::vector<G4int> indicesDamage_merged;
	std::vector<G4int> damageIndices_direct_copy = pDamageIndices_direct;
	std::vector<G4int> damageIndices_indirect_copy = pDamageIndices_indirect;

	for (G4int element : damageIndices_direct_copy) {
		if (!IsElementInVector(element, indicesDamage_merged))
			indicesDamage_merged.push_back(element);
		else {
			fDoubleCountsDD++;
			G4cout << "\tDirect damage site already recorded: " << element << G4endl;
			RemoveElementFromVector(element, pDamageIndices_direct);
		}
	}

	for (G4int element : damageIndices_indirect_copy) {
		if (!IsElementInVector(element, indicesDamage_merged))
			indicesDamage_merged.push_back(element);
		else {
			fDoubleCountsDI++;
			G4cout << "\tIndirect damage site already recorded: " << element << G4endl;
			RemoveElementFromVector(element, pDamageIndices_indirect);
		}
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
	fTotalSSB_direct = 0;
	fTotalSSB_indirect = 0;

	fTotalBD = 0;
	fTotalBD_direct = 0;
	fTotalBD_indirect = 0;

	fTotalDSB = 0;
	fTotalDSB_hybrid = 0;
	fTotalDSB_direct = 0;
	fTotalDSB_indirect = 0;

	fTotalComplexDSB = 0;
	fTotalComplexDSB_direct = 0;
	fTotalComplexDSB_indirect = 0;
	fTotalComplexDSB_hybrid = 0;

	fTotalNonDSBCluster = 0;
	fTotalNonDSBCluster_direct = 0;
	fTotalNonDSBCluster_indirect = 0;
	fTotalNonDSBCluster_hybrid = 0;

	fDoubleCountsDD = 0;
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
std::vector<G4int> ScoreClusteredDNADamage::RecordDSB(G4int pDamageCause)
{
	if (pDamageCause == fIdDirect) { // direct
		fIndicesSSB1 = &fIndicesSSB1_direct;
		fIndicesSSB2 = &fIndicesSSB2_direct;
	}
	else if (pDamageCause == fIdIndirect) { // indirect
		fIndicesSSB1 = &fIndicesSSB1_indirect;
		fIndicesSSB2 = &fIndicesSSB2_indirect;
	}
	else if (pDamageCause == fIdHybrid) { // hybrid
		fIndicesSSB1_merged = MergeDamageIndices(fIndicesSSB1_direct, fIndicesSSB1_indirect);
		fIndicesSSB2_merged = MergeDamageIndices(fIndicesSSB2_direct, fIndicesSSB2_indirect);
		fIndicesSSB1 = &fIndicesSSB1_merged;
		fIndicesSSB2 = &fIndicesSSB2_merged;
	}
	else {
		G4cerr << "Error: (While scoring DSB) The following integer damage cause label is unrecognized: " << pDamageCause << G4endl;
		exit(0);
	}

	std::vector<G4int> indicesDSB1D;
	if (fIndicesSSB1->size() == 0 || fIndicesSSB2->size() == 0)
		return indicesDSB1D;

	if (fIndicesSSB1->size() > 1 && !std::is_sorted(fIndicesSSB1->begin(), fIndicesSSB1->end()))
		std::sort(fIndicesSSB1->begin(), fIndicesSSB1->end());
	if (fIndicesSSB2->size() > 1 && !std::is_sorted(fIndicesSSB2->begin(), fIndicesSSB2->end()))
		std::sort(fIndicesSSB2->begin(), fIndicesSSB2->end());

	std::vector<G4int>::iterator site1 = fIndicesSSB1->begin();
	std::vector<G4int>::iterator site2 = fIndicesSSB2->begin();

	// Proceed until have completely processed SSBs in either strand
	while (site1 != fIndicesSSB1->end() && site2 != fIndicesSSB2->end()) {
		G4int siteDiff = *site2 - *site1; // separation in number of bp

		G4bool isDSB = abs(siteDiff) <= fThresDistForDSB;
		G4bool isDSBrecorded = isDSB;

		G4bool isHybrid = false;
		G4bool isSite1Direct = IsElementInVector(*site1, fIndicesSSB1_direct);
		G4bool isSite1Indirect = IsElementInVector(*site1, fIndicesSSB1_indirect);
		G4bool isSite2Direct = IsElementInVector(*site2, fIndicesSSB2_direct);
		G4bool isSite2Indirect = IsElementInVector(*site2, fIndicesSSB2_indirect);
		if (pDamageCause == fIdHybrid) { // hybrid
			if (isSite1Direct && isSite1Indirect) {
				G4cerr << "Error: (While scoring DSB) CONTRADICTION! site1: " << *site1 << " is both direct and indirect!" << G4endl;
				exit(0);
			}
			if (isSite2Direct && isSite2Indirect) {
				G4cerr << "Error: (While scoring DSB) CONTRADICTION! site2: " << *site2 << " is both direct and indirect!" << G4endl;
				exit(0);
			}
			isHybrid = (isSite1Direct && isSite2Indirect) || (isSite1Indirect && isSite2Direct);
			isDSBrecorded = isHybrid;
		}

		// Damage in site 2 is within range of site 1 to count as DSB (either before or after)
		if (isDSB && isDSBrecorded) {
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

			if (isHybrid) {
				G4cout << "\tRecorded for hybrid damage: " << *site1 << ", " << *site2 << G4endl;
				if (isSite1Direct && isSite2Indirect) {
					RemoveElementFromVector(*site1, fIndicesSSB1_direct);
					RemoveElementFromVector(*site2, fIndicesSSB2_indirect);
				}
				else if (isSite1Indirect && isSite2Direct) {
					RemoveElementFromVector(*site1, fIndicesSSB1_indirect);
					RemoveElementFromVector(*site2, fIndicesSSB2_direct);
				}
				else {
					G4cerr << "Error: (While scoring DSB) False count of hybrid damage!" << G4endl;
					exit(0);
				}
			}
			site1 = fIndicesSSB1->erase(site1);
			site2 = fIndicesSSB2->erase(site2);
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
	std::vector<std::array<G4int,3>>::iterator site = fIndicesSimple.begin()+1;

	DamageCluster cluster;

	G4bool newCluster = true; // flag to indicate a new cluster is formed
	G4bool buildingCluster = false; // flag to indicate building on an existing cluster

	// Loop over all damages arranged in order along the strand (SSBs, BDs, and DSBs)
	while (site != fIndicesSimple.end()) {
		std::array<G4int,3> sitePrev = *(site-1); // previous damage site
		std::array<G4int,3> siteCur = *site; // current damage site

		// If damage sites are close enough to form a cluster
		if ((siteCur[0]-sitePrev[0]) <= fThresDistForCluster) {
			// A new cluster is being formed, so the "previous" damage must be added to the start
			// of the cluster
			if (newCluster) {
				AddDamageToCluster(cluster,sitePrev[0],sitePrev[1],sitePrev[2],newCluster);
				newCluster = false;
				buildingCluster = true;
			}

			// Add the current damage to the cluster
			AddDamageToCluster(cluster,siteCur[0],siteCur[1],siteCur[2],newCluster);
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
												G4int damageType, G4int damageCause, G4bool newCluster) {
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
	if (damageType == fIdSSB) { // SSBs
		cluster.numSSB++;
		fTotalSSB--;
		if (damageCause == fIdDirect){
			cluster.numSSB_direct++;
			fTotalSSB_direct--;
		}
		else if (damageCause == fIdIndirect){
			cluster.numSSB_indirect++;
			fTotalSSB_indirect--;
		}
		else {
			G4cout << "Error: (While scoring clustered DMA damage) The following integer damage cause label is unrecognized: " << damageCause << G4endl;
			exit(0);
		}
	}
	else if (damageType == fIdBD) { // BDs
		cluster.numBD++;
		fTotalBD--;
		if (damageCause == fIdDirect){
			cluster.numBD_direct++;
			fTotalBD_direct--;
		}
		else if (damageCause == fIdIndirect){
			cluster.numBD_indirect++;
			fTotalBD_indirect--;
		}
		else {
			G4cout << "Error: (While scoring clustered DMA damage) The following integer damage cause label is unrecognized: " << damageCause << G4endl;
			exit(0);
		}
	}
	else if (damageType == fIdDSB) { // DSBs
		cluster.numDSB++;
		fTotalDSB--;
		if (damageCause == fIdDirect){
			cluster.numDSB_direct++;
			fTotalDSB_direct--;
		}
		else if (damageCause == fIdIndirect){
			cluster.numDSB_indirect++;
			fTotalDSB_indirect--;
		}
		else if (damageCause == fIdHybrid){
			cluster.numDSB_hybrid++;
			fTotalDSB_hybrid--;
		}
		else {
			G4cout << "Error: (While scoring clustered DMA damage) The following integer damage cause label is unrecognized: " << damageCause << G4endl;
			exit(0);
		}
	}
	// Throw an error if an unrecognized damage type is added
	else {
		G4cout << "An error has arisen while scoring clustered DNA damage." << G4endl;
		G4cout << "The following integer damage type label is unrecognized:" << G4endl;
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
	G4bool hasDirectSSB = cluster.numSSB_direct > 0;
	G4bool hasIndirectSSB = cluster.numSSB_indirect > 0;

	G4bool hasDirectBD = cluster.numBD_direct > 0;
	G4bool hasIndirectBD = cluster.numBD_indirect > 0;

	G4bool hasDirectDSB = cluster.numDSB_direct > 0;
	G4bool hasIndirectDSB = cluster.numDSB_indirect > 0;
	G4bool hasHybridDSB = cluster.numDSB_hybrid > 0;

	// Handle Simple DSB (i.e. not a cluster, so record nothing & reset)
	if (cluster.numDSB == 2 && cluster.numSSB == 0 && cluster.numBD == 0) {
		fTotalDSB += 2;

		if (hasDirectDSB && !hasIndirectDSB && !hasHybridDSB)
			fTotalDSB_direct += 2;
		else if (!hasDirectDSB && hasIndirectDSB && !hasHybridDSB)
			fTotalDSB_indirect += 2;
		else
			fTotalDSB_hybrid += 2;

		cluster = DamageCluster();
	}
	// Handle Complex DSB
	else if (cluster.numDSB > 0) {
		cluster.numDSB = cluster.numDSB/2;

		fComplexDSBSizes.push_back(cluster.end - cluster.start + 1);

		fComplexDSBNumSSB.push_back(cluster.numSSB);
		fComplexDSBNumSSB_direct.push_back(cluster.numSSB_direct);
		fComplexDSBNumSSB_indirect.push_back(cluster.numSSB_indirect);

		fComplexDSBNumBD.push_back(cluster.numBD);
		fComplexDSBNumBD_direct.push_back(cluster.numBD_direct);
		fComplexDSBNumBD_indirect.push_back(cluster.numBD_indirect);

		fComplexDSBNumDSB.push_back(cluster.numDSB);
		fComplexDSBNumDSB_direct.push_back(cluster.numDSB_direct);
		fComplexDSBNumDSB_indirect.push_back(cluster.numDSB_indirect);
		fComplexDSBNumDSB_hybrid.push_back(cluster.numDSB_hybrid);

		fComplexDSBNumDamage.push_back(cluster.numSSB + cluster.numBD + cluster.numDSB);

		fTotalComplexDSB++;

		if ( hasDirectDSB && !hasIndirectSSB && !hasIndirectBD && !hasIndirectDSB && !hasHybridDSB)
			fTotalComplexDSB_direct++;
		else if (!hasDirectSSB && !hasDirectBD && !hasDirectDSB && hasIndirectDSB && !hasHybridDSB)
			fTotalComplexDSB_indirect++;
		else
			fTotalComplexDSB_hybrid++;

		cluster = DamageCluster();
	}
	// Handle Non-DSB cluster
	else {
		fNonDSBClusterSizes.push_back(cluster.end - cluster.start + 1);

		fNonDSBClusterNumSSB.push_back(cluster.numSSB);
		fNonDSBClusterNumSSB_direct.push_back(cluster.numSSB_direct);
		fNonDSBClusterNumSSB_indirect.push_back(cluster.numSSB_indirect);

		fNonDSBClusterNumBD.push_back(cluster.numBD);
		fNonDSBClusterNumBD_direct.push_back(cluster.numBD_direct);
		fNonDSBClusterNumBD_indirect.push_back(cluster.numBD_indirect);

		fNonDSBClusterNumDamage.push_back(cluster.numSSB + cluster.numBD);

		fTotalNonDSBCluster++;

		if ( (hasDirectSSB || hasDirectBD) && !hasIndirectSSB && !hasIndirectBD)
			fTotalNonDSBCluster_direct++;
		else if (!hasDirectSSB && !hasDirectBD && (hasIndirectSSB || hasIndirectBD) )
			fTotalNonDSBCluster_indirect++;
		else
			fTotalNonDSBCluster_hybrid++;

		cluster = DamageCluster();
	}
}


//--------------------------------------------------------------------------------------------------
// Combine class member vectors containing various types of damages into a single, ordered, vector
// of all damges in a DNA fibre (both strands). Each vector element is a 2-element array containing
// (i) the bp index of the damage and (ii) an integer indicating the type of damage (BD or SSB).
//--------------------------------------------------------------------------------------------------
std::vector<std::array<G4int,3>> ScoreClusteredDNADamage::CombineSimpleDamage() {
	// G4cout << "COMBINING SIMPLE DAMAGE" << G4endl;
	std::vector<std::array<G4int,3>> indicesSimple;

	// Fake data for testing
	// fIndicesSSB1 = {3,4,5,7,9};
	// fIndicesBD1 = {1,4,7};
	// fIndicesSSB2 = {1,12};
	// fIndicesBD2 = {6,7,10,12};
	// fIndicesDSB = {3,5,11,14};

	// SSBs in strand 1 (direct)
	std::vector<G4int>::iterator iter = fIndicesSSB1_direct.begin();
	while (iter != fIndicesSSB1_direct.end()) {
		indicesSimple.push_back({*iter,fIdSSB,fIdDirect});
		iter++;
	}

	// SSBs in strand 1 (indirect)
	iter = fIndicesSSB1_indirect.begin();
	while (iter != fIndicesSSB1_indirect.end()) {
		indicesSimple.push_back({*iter,fIdSSB,fIdIndirect});
		iter++;
	}

	// BDs in strand 1 (direct)
	iter = fIndicesBD1_direct.begin();
	while (iter != fIndicesBD1_direct.end()) {
		indicesSimple.push_back({*iter,fIdBD,fIdDirect});
		iter++;
	}

	// BDs in strand 1 (indirect)
	iter = fIndicesBD1_indirect.begin();
	while (iter != fIndicesBD1_indirect.end()) {
		indicesSimple.push_back({*iter,fIdBD,fIdIndirect});
		iter++;
	}

	// SSBs in strand 2 (direct)
	iter = fIndicesSSB2_direct.begin();
	while (iter != fIndicesSSB2_direct.end()) {
		indicesSimple.push_back({*iter,fIdSSB,fIdDirect});
		iter++;
	}

	// SSBs in strand 2 (indirect)
	iter = fIndicesSSB2_indirect.begin();
	while (iter != fIndicesSSB2_indirect.end()) {
		indicesSimple.push_back({*iter,fIdSSB,fIdIndirect});
		iter++;
	}

	// BDs in strand 2 (direct)
	iter = fIndicesBD2_direct.begin();
	while (iter != fIndicesBD2_direct.end()) {
		indicesSimple.push_back({*iter,fIdBD,fIdDirect});
		iter++;
	}

	// BDs in strand 2 (indirect)
	iter = fIndicesBD2_indirect.begin();
	while (iter != fIndicesBD2_indirect.end()) {
		indicesSimple.push_back({*iter,fIdBD,fIdIndirect});
		iter++;
	}

	// DSBs (hybrid)
	iter = fIndicesDSB_hybrid.begin();
	while (iter != fIndicesDSB_hybrid.end()) {
		indicesSimple.push_back({*iter,fIdDSB,fIdHybrid});
		iter++;
	}

	// DSBs (direct)
	iter = fIndicesDSB_direct.begin();
	while (iter != fIndicesDSB_direct.end()) {
		indicesSimple.push_back({*iter,fIdDSB,fIdDirect});
		iter++;
	}

	// DSBs (indirect)
	iter = fIndicesDSB_indirect.begin();
	while (iter != fIndicesDSB_indirect.end()) {
		indicesSimple.push_back({*iter,fIdDSB,fIdIndirect});
		iter++;
	}

	if (indicesSimple.size() > 1)
		std::sort(indicesSimple.begin(), indicesSimple.end()); // sorts by first element in each vector item by default (i.e site index)

	// If want to output full list of simple damages:
	// if (indicesSimple.size() > 0) {
	// 	G4cout << "------------------------------------------" << G4endl;
	// 	for (int i = 0; i < indicesSimple.size(); i++) {
	// 		G4cout << "site = " << indicesSimple[i][0] << ", type = " << indicesSimple[i][1] << ", cause = " << indicesSimple[i][2] << G4endl;
	// 	}
	// 	G4cout << "------------------------------------------" << G4endl;
	// }

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

	G4cout << "Total # of SSBs = " << fTotalSSB << G4endl;
	G4cout << "# of direct SSBs = " << fTotalSSB_direct << G4endl;
	G4cout << "# of indirect SSBs = " << fTotalSSB_indirect << G4endl;

	G4cout << "Total # of BDs = " << fTotalBD << G4endl;
	G4cout << "# of direct BDs = " << fTotalBD_direct << G4endl;
	G4cout << "# of indirect BDs = " << fTotalBD_indirect << G4endl;

	G4cout << "Total # of simple DSBs = " << fTotalDSB << G4endl;
	G4cout << "# of hybrid DSBs = " << fTotalDSB_hybrid << G4endl;
	G4cout << "# of direct DSBs = " << fTotalDSB_direct << G4endl;
	G4cout << "# of indirect DSBs = " << fTotalDSB_indirect << G4endl;

	G4cout << "Total # of complex DSBs = " << fTotalComplexDSB << G4endl;
	G4cout << "# of hybrid complex DSBs = " << fTotalComplexDSB_hybrid << G4endl;
	G4cout << "# of direct complex DSBs = " << fTotalComplexDSB_direct << G4endl;
	G4cout << "# of indirect complex DSBs = " << fTotalComplexDSB_indirect << G4endl;
	for (G4int i=0; i < fTotalComplexDSB; i++) {
		G4cout << "-----------------------------" << G4endl;
		G4cout << "\tSize = " << fComplexDSBSizes[i] << G4endl;
		G4cout << "\tNum SSB = " << fComplexDSBNumSSB[i] << G4endl;
		G4cout << "\t Num direct SSB = " << fComplexDSBNumSSB_direct[i] << G4endl;
		G4cout << "\t Num indirect SSB = " << fComplexDSBNumSSB_indirect[i] << G4endl;
		G4cout << "\tNum BD = " << fComplexDSBNumBD[i] << G4endl;
		G4cout << "\t Num direct BD = " << fComplexDSBNumBD_direct[i] << G4endl;
		G4cout << "\t Num indirect BD = " << fComplexDSBNumBD_indirect[i] << G4endl;
		G4cout << "\tNum DSB = " << fComplexDSBNumDSB[i] << G4endl;
		G4cout << "\t Num hybrid DSB = " << fComplexDSBNumDSB_hybrid[i] << G4endl;
		G4cout << "\t Num direct DSB = " << fComplexDSBNumDSB_direct[i] << G4endl;
		G4cout << "\t Num indirect DSB = " << fComplexDSBNumDSB_indirect[i] << G4endl;
		G4cout << "\tNum damage = " << fComplexDSBNumDamage[i] << G4endl;
	}

	G4cout << "Total # of non-DSB clusters = " << fTotalNonDSBCluster << G4endl;
	G4cout << "# of hybrid non-DSB clusters = " << fTotalNonDSBCluster_hybrid << G4endl;
	G4cout << "# of direct non-DSB clusters = " << fTotalNonDSBCluster_direct << G4endl;
	G4cout << "# of indirect non-DSB clusters = " << fTotalNonDSBCluster_indirect << G4endl;
	for (G4int i=0; i < fTotalNonDSBCluster; i++) {
		G4cout << "-----------------------------" << G4endl;
		G4cout << "\tSize = " << fNonDSBClusterSizes[i] << G4endl;
		G4cout << "\tNum SSB = " << fNonDSBClusterNumSSB[i] << G4endl;
		G4cout << "\t Num direct SSB = " << fNonDSBClusterNumSSB_direct[i] << G4endl;
		G4cout << "\t Num indirect SSB = " << fNonDSBClusterNumSSB_indirect[i] << G4endl;
		G4cout << "\tNum BD = " << fNonDSBClusterNumBD[i] << G4endl;
		G4cout << "\t Num direct BD = " << fNonDSBClusterNumBD_direct[i] << G4endl;
		G4cout << "\t Num indirect BD = " << fNonDSBClusterNumBD_indirect[i] << G4endl;
		G4cout << "\tNum damage = " << fNonDSBClusterNumDamage[i] << G4endl;
	}

	G4cout << "Total # of direct-direct double counts = " << fDoubleCountsDD << G4endl;
	G4cout << "Total # of direct-indirect double counts = " << fDoubleCountsDI << G4endl;
	G4cout << "Total # of indirect-indirect double counts = " << fDoubleCountsII << G4endl;
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
