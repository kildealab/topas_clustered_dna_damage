// Scorer for ScoreSDD
//
//**************************************************************************************************
// Author: Nicolas Desjardins, Logan Montgomery, James Manalad
// Based on code written by:
//      - H Zhu et al. (2020). DOI:10.1667/RR15531.1
//      - J Schuemann et al. (2019). DOI:10.1667/RR15226.1
//		- A McNamara et al. (2018). DOI:10.1088/1361-6560/aad8eb
//      - E Delage et al. (2015). DOI:10.1016/j.cpc.2015.02.026
//
// This class is a custom Topas ntuple scorer that records DNA damage induced in a DNA
// model in the SDD format (VoxelizedNuclearDNA).
//**************************************************************************************************

#include "ScoreSDD.hh"
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

#include <iostream>
#include <stdio.h>




//--------------------------------------------------------------------------------------------------
// Constructor. Initialize member variables using a variety of methods. Specify which data is output
// to the main data output file.
//--------------------------------------------------------------------------------------------------
ScoreSDD::ScoreSDD(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
											   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
:TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)

{
	//----------------------------------------------------------------------------------------------
	// Initialize member variables
	//----------------------------------------------------------------------------------------------
	ResolveParams(); // initialize some member variables using Topas parameter file

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

	G4cout << "Conversion Factor: " << GetMaterial("G4_WATER")->GetDensity() * fComponentVolume * gray << G4endl;

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

}

//--------------------------------------------------------------------------------------------------
// Destructor
//--------------------------------------------------------------------------------------------------
ScoreSDD::~ScoreSDD() {
}


//--------------------------------------------------------------------------------------------------
// Read in parameters from the Topas parameter file & save values in some member variables. Default
// values are provided in case the user does not specify a value.
//--------------------------------------------------------------------------------------------------
void ScoreSDD::ResolveParams() {
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
	// Specify whether to report damage on a per-fibre basis (vs. grouping together). Note that
	// damage in separate fibers are never considered together when clustering damage.
	//----------------------------------------------------------------------------------------------
	if ( fPm->ParameterExists(GetFullParmName("RecordDamagePerFiber")))
		fRecordDamagePerFiber = fPm->GetBooleanParameter(GetFullParmName("RecordDamagePerFiber"));
	else
		fRecordDamagePerFiber = false;


	//----------------------------------------------------------------------------------------------
	//Name of output file for all events and name out SDD ouput file
	//----------------------------------------------------------------------------------------------

	if ( fPm->ParameterExists(GetFullParmName("FileAllEvents")))
		fFileAllEvents = fPm->GetStringParameter(GetFullParmName("FileAllEvents"));
	else
		fFileAllEvents = "FileAllEvents.txt";


	if ( fPm->ParameterExists(GetFullParmName("SDDOutputFile")))
		fSDDOutputFile = fPm->GetStringParameter(GetFullParmName("SDDOutputFile"));
	else
		fSDDOutputFile = "SDDOutputFile.sdd";


	if ( fPm->ParameterExists(GetFullParmName("RealDoseFile")))
		fRealDoseFile = fPm->GetStringParameter(GetFullParmName("RealDoseFile"));
	else
		fRealDoseFile = "RealDoseFile.txt";


	if ( fPm->ParameterExists(GetFullParmName("MapFile")))
		fMapFile = fPm->GetStringParameter(GetFullParmName("MapFile"));
	else
		fMapFile= "mapping.csv";


	if ( fPm->ParameterExists(GetFullParmName("MeanEnergy")))
		MeanEnergy = fPm->GetStringParameter(GetFullParmName("MeanEnergy"));
	else
		MeanEnergy = "Non specified";


	if ( fPm->ParameterExists(GetFullParmName("PrimaryParticle")))
		PrimaryParticle = fPm->GetStringParameter(GetFullParmName("PrimaryParticle"));
	else
		PrimaryParticle = "proton";



	if (PrimaryParticle == "proton")
		IncidentParticles = "2012";

	else if (PrimaryParticle  == "e-")
		IncidentParticles = "11";

	else if (PrimaryParticle == "alpha")
		IncidentParticles = "100204";

	else 
		IncidentParticles = "22";





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

	// Check for chemistry module
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

		// Damage probabilities of radiolytic species to induce SSB and BD
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
		// Species killed by DNA volumes
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
		// Species killed by histone volumes
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
	// Toggles to score direct or indirect damage
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


	//This let the user wether of not to include damage blocks that do not contain a DSB
	if ( fPm->ParameterExists(GetFullParmName("DSBOnly")))
		DsbOnly = fPm->GetBooleanParameter(GetFullParmName("DSBOnly"));
	else
		DsbOnly = false;

	//----------------------------------------------------------------------------------------------
	// Material of DNA residue and histone volumes in which to score
	//----------------------------------------------------------------------------------------------
	G4String DNAMaterialName = "G4_WATER_DNA";
	if ( fPm->ParameterExists(GetFullParmName("DNAMaterialName")))
		DNAMaterialName = fPm->GetStringParameter(GetFullParmName("DNAMaterialName"));
	fDNAMaterial = GetMaterial(DNAMaterialName);

	G4String HistoneMaterialName = "G4_WATER_HISTONE";
	if ( fPm->ParameterExists(GetFullParmName("HistoneMaterialName")))
		HistoneMaterialName = fPm->GetStringParameter(GetFullParmName("HistoneMaterialName"));
	fHistoneMaterial = GetMaterial(HistoneMaterialName);

	//----------------------------------------------------------------------------------------------
	// Parameters not specific to this extension. Can't/don't need to use GetFullParmName()
	//----------------------------------------------------------------------------------------------
	fNumberOfThreads = fPm->GetIntegerParameter("Ts/NumberOfThreads");

	GetChromoMap();
}  


//--------------------------------------------------------------------------------------------------
// Erase contents of various output files 
//--------------------------------------------------------------------------------------------------
void ScoreSDD::ClearOutputFiles() {
	std::ofstream fileToClear;

	fileToClear.open(fFileAllEvents, std::ofstream::trunc);
 	fileToClear.close();

 	fileToClear.open(fSDDOutputFile, std::ofstream::trunc);
 	fileToClear.close();

 	fileToClear.open(fRealDoseFile, std::ofstream::trunc);
 	fileToClear.close();
}


//--------------------------------------------------------------------------------------------------
// Calculate cubic volume of component attached to the scorer. Use parameter values to perform
// calculation according to the shape of the volume
//--------------------------------------------------------------------------------------------------
G4double ScoreSDD::CalculateComponentVolume() {
	G4double componentVolume = pow(2*fVoxelSideLength*fNumVoxelsPerSide,3);

	return componentVolume;
}

//--------------------------------------------------------------------------------------------------
// Handle conversion of dose threshold to energy threshold using the cubic volume and density of
// the geometry component. The resulting energy threshold is divided evenly amongst the worker
// threads. E.g. If total energy threshold is 1 MeV and there are 4 worker threads. The returned
// value is 0.25 MeV.
//--------------------------------------------------------------------------------------------------
G4double ScoreSDD::ConvertDoseThresholdToEnergy() {
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
G4bool ScoreSDD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	fNumProcessHitsCalls++; // for debugging purposes
	G4double edep = aStep->GetTotalEnergyDeposit(); // In eV;
	fTotalEdep += edep; // running sum of energy deposition in entire volume
	// Material filtering of pre-step and post-step (only proceed if a sensitive volume is involved)
	G4Material* materialPreStep = aStep->GetPreStepPoint()->GetMaterial();
	G4bool isPreStepDNAMaterial = (materialPreStep == fDNAMaterial);
	G4bool isPreStepHistoneMaterial = (materialPreStep == fHistoneMaterial);
	G4Material* materialPostStep = aStep->GetPostStepPoint()->GetMaterial();
	G4bool isPostStepDNAMaterial = (materialPostStep == fDNAMaterial);
	if (!isPreStepDNAMaterial && !isPostStepDNAMaterial && !isPreStepHistoneMaterial) {
		return false;
	}

	// G4Touchable provides access to parent volumes, etc.
	G4TouchableHistory* touchable;
	if (isPreStepDNAMaterial || isPreStepHistoneMaterial) {
		touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
	}
	else if (isPostStepDNAMaterial) {
		touchable = (G4TouchableHistory*)(aStep->GetPostStepPoint()->GetTouchable());
	}
	else {
		G4cerr << "Error: Material check failed. PreStep: " << materialPreStep << " PostStep: " << materialPostStep << G4endl;
		exit(0);
	}

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
	G4int volID = touchable->GetVolume()->GetCopyNo();
	G4int strandID = volID / 1000000;
	G4int residueID = (volID - (strandID*1000000)) / 100000;
	G4int bpID = volID - (strandID*1000000) - (residueID*100000);

	// Determines whether track is physical or chemical
	G4int trackID = aStep->GetTrack()->GetTrackID();

	//Stores the position (x,y,z) of the pre-step 
	G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
	G4ThreeVector pos2 = aStep->GetPostStepPoint()->GetPosition();
	//----------------

	//----------------------------------------------------------------------------------------------
	// If this hit deposits energy (in sensitive DNA volume), update the appropriate energy deposition
	// map
	//----------------------------------------------------------------------------------------------
	if (fIncludeDirectDamage && edep > 0 && trackID >= 0 && isPreStepDNAMaterial) { // energy deposition should be from physical tracks
		//------------------------------------------------------------------------------------------
		// Use the DNA strand ID, residue ID, and nucleotide ID to increment the energy deposited
		// in the appropriate energy deposition map. Maps are indexed as follows:
		// First index specifies the voxel
		// Second index specifies DNA fibre
		// Third index specifies the bp index
		//------------------------------------------------------------------------------------------
		if ( strandID == 0 ){ // first strand
			if (residueID == fVolIdPhosphate || residueID == fVolIdDeoxyribose) {
				fMapEdepStrand1Backbone[fVoxelID][fFiberID][bpID].first += edep;
				fMapEdepStrand1Backbone[fVoxelID][fFiberID][bpID].second = pos;
	

			}
			else if (residueID == fVolIdBase) {
				fMapEdepStrand1Base[fVoxelID][fFiberID][bpID].first += edep;
				fMapEdepStrand1Base[fVoxelID][fFiberID][bpID].second = pos;
				
		
			}
		}
		else if ( strandID == 1 ) { // second strand
			if (residueID == fVolIdPhosphate || residueID == fVolIdDeoxyribose) {
				fMapEdepStrand2Backbone[fVoxelID][fFiberID][bpID].first += edep;
				fMapEdepStrand2Backbone[fVoxelID][fFiberID][bpID].second = pos;
			}
			else if (residueID == fVolIdBase) {
				fMapEdepStrand2Base[fVoxelID][fFiberID][bpID].first += edep;
				fMapEdepStrand2Base[fVoxelID][fFiberID][bpID].second = pos;
			}
		}
		else {
			G4cerr << "Error: (While scoring direct damage) The following strandID is unrecognized: " << strandID << G4endl;
			exit(0);
		}
		return true;
	}

	// Indirect damage
	else if (fIncludeIndirectDamage && trackID < 0) { // chemical tracks
		// Get molecule info
		G4int moleculeID = GetMolecule(aStep->GetTrack())->GetMoleculeID();

		// Kill species generated inside DNA volumes and histones by not letting them exit
		G4Material* materialTrackVertex = aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetMaterial();
		G4bool isPreStepInTrackVertexVolume = (materialPreStep == materialTrackVertex);
		G4bool isPostStepInNewVolume	= (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);
		if ( isPostStepInNewVolume && isPreStepInTrackVertexVolume && (isPreStepDNAMaterial || isPreStepHistoneMaterial) ) {
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return false;
		}

		// Determine if damage is inflicted
		G4bool isDamaged = IsDamageInflicted(moleculeID, residueID);
		if (isDamaged && isPostStepDNAMaterial && isPostStepInNewVolume && !isPreStepDNAMaterial && !isPreStepHistoneMaterial) {
			// Check which damage map to update
			if ( strandID == 0 ) { // first strand
				if (residueID == fVolIdPhosphate || residueID == fVolIdDeoxyribose) { // backbone damage
					fMapIndDamageStrand1Backbone[fVoxelID][fFiberID][bpID] = pos;
					//------------------
				}
				else if (residueID == fVolIdBase){ // base damage
					fMapIndDamageStrand1Base[fVoxelID][fFiberID][bpID] = pos;
				}
			}
			else if ( strandID == 1 ) { // second strand
				if (residueID == fVolIdPhosphate || residueID == fVolIdDeoxyribose) { // backbone damage
					fMapIndDamageStrand2Backbone[fVoxelID][fFiberID][bpID] = pos;
				}
				else if (residueID == fVolIdBase){ // base damage
					fMapIndDamageStrand2Base[fVoxelID][fFiberID][bpID] = pos;
				}
			}
			else {
				G4cerr << "Error: (While scoring indirect damage) The following strandID is unrecognized: " << strandID << G4endl;
				exit(0);
			}

		}

		// Kill certain species interacting with DNA volumes
		G4bool isSpeciesToKillByDNA = IsElementInVector(moleculeID, fSpeciesToKillByDNAVolumes);
		if (isPostStepDNAMaterial && isPostStepInNewVolume && isSpeciesToKillByDNA && !isPreStepDNAMaterial && !isPreStepHistoneMaterial) {
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return false;
		}

		// Kill certain species diffusing in histone volumes
		if (fHistonesAsScavenger && isPreStepHistoneMaterial) {
			G4bool isSpeciesToKillByHistone = IsElementInVector(moleculeID, fspeciesToKillByHistones);
			G4bool isPreStepInNewVolume	= (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary);
			if (isSpeciesToKillByHistone && isPreStepInNewVolume) {
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
G4bool ScoreSDD::IsElementInVector(G4int pElement, std::vector<G4int> &pVector) {
	if (std::find(pVector.begin(), pVector.end(), pElement) != pVector.end())
		return true;
	else
		return false;
}


//--------------------------------------------------------------------------------------------------
// This helper method removes a specified element from a vector.
//--------------------------------------------------------------------------------------------------
void ScoreSDD::RemoveElementFromVector(G4int pElement, std::vector<G4int> &pVector) {
	if (!IsElementInVector(pElement, pVector)) {
		G4cout << "Element (" << pElement << ") not in vector" << G4endl;
		Print1DVectorContents(pVector);
		return;
	}
	pVector.erase(std::remove(pVector.begin(), pVector.end(), pElement), pVector.end());
}


//--------------------------------------------------------------------------------------------------
// This helper method checks whether an indirect damage has been induced given the moleculeID.
//--------------------------------------------------------------------------------------------------
G4bool ScoreSDD::IsDamageInflicted(G4int pMoleculeID, G4int pDNAVolumeID) {
	G4float prob_damage = 0.;

	if ( fMoleculeDamageProb_SSB.find(pMoleculeID) == fMoleculeDamageProb_SSB.end() && fMoleculeDamageProb_BD.find(pMoleculeID) == fMoleculeDamageProb_BD.end()) {
  	// pMoleculeID not found in fMoleculeDamageProb_SSB and fMoleculeDamageProb_BD
		G4cerr << "\tmoleculeID NOT LISTED: " << pMoleculeID << G4endl;
		exit(0);
	}

	// Check if SSB or BD
	if (pDNAVolumeID == fVolIdPhosphate || pDNAVolumeID == fVolIdDeoxyribose)
		prob_damage = fMoleculeDamageProb_SSB[pMoleculeID];
	else if (pDNAVolumeID == fVolIdBase)
		prob_damage = fMoleculeDamageProb_BD[pMoleculeID];

	if (prob_damage == 0.)
		return false;

	// Generate random float between 0 and 1
	G4float prob_random = G4UniformRand();

	if (prob_random <= prob_damage)
		return true;
	else
		return false;
}


//--------------------------------------------------------------------------------------------------
// This helper method prints useful information about a given G4step.
//--------------------------------------------------------------------------------------------------
void ScoreSDD::PrintStepInfo(G4Step* aStep) { // for debugging
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
	G4cout << "\tTrack length: " << aStep->GetTrack()->GetTrackLength() << G4endl;
	G4cout << "\tTrack status: " << aStep->GetTrack()->GetTrackStatus() << G4endl;
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
void ScoreSDD::UserHookForEndOfEvent() {
	fEventID = GetEventID();
	fThreadID = G4Threading::G4GetThreadId();

	fNumEvents++;

	G4String outputFileName = fFileAllEvents;
	std::ofstream outFile(outputFileName, std::ios_base::app);
	outFile.close();

	NewEvent = true;
	RecordDamageSDD();
	ResetMemberVariables(); // Necessary to reset variables before proceeding to next event
	

	// Check if dose threshold has been met
	if (fUseDoseThreshold && fTotalEdep > fEnergyThreshold) {
		G4cout << "Aborting worker #" << G4Threading::G4GetThreadId() << " because dose threshold has been met" << G4endl;
		std::ofstream outFile(fRealDoseFile, std::ios_base::app);
			outFile << fTotalEdep / (GetMaterial("G4_WATER")->GetDensity() * fComponentVolume * gray) << G4endl;
			outFile.close();
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
void ScoreSDD::UserHookForEndOfRun() {
	// fEventID = GetEventID();
	fThreadID = G4Threading::G4GetThreadId();

	OutputSDDHeader();
	PrintSDDdata();
}


void ScoreSDD::OutputEvent(G4int ModeOfAction, G4int PositionOnBp, G4ThreeVector position, G4int voxelID, G4int index){

	//----------------------------------------------------------------------------------------------
	// Data file
	//----------------------------------------------------------------------------------------------
	G4String outputFileName = fFileAllEvents;
	std::ofstream outFile(outputFileName, std::ios_base::app);

	// Catch file I/O error
	if (!outFile.good()) {
		G4cerr << "Topas is exiting due to a serious error in file output." << G4endl;
		G4cerr << "Output file: " << outputFileName << " cannot be opened" << G4endl;
		fPm->AbortSession(1);
	}

	G4int Chromosome = ChromoMapping[fVoxelID][0];
	G4int arm = ChromoMapping[fVoxelID][2];
	G4int BpPerFiber = 18000;
	G4int BpPerVoxel = 20*BpPerFiber;

	G4int BpInChromo = ChromoMapping[fVoxelID][1]*BpPerVoxel + fFiberID*BpPerFiber + index;

	// Record data
	outFile << ModeOfAction  << fDelimiter << position[0]/um << fDelimiter << position[1]/um  << fDelimiter << position[2]/um 
	<< fDelimiter << Chromosome << fDelimiter << BpInChromo << fDelimiter << PositionOnBp << G4endl;
	
	outFile.close();
}


//--------------------------------------------------------------------------------------------------
// This method transfers information from worker threads to the master thread.
//--------------------------------------------------------------------------------------------------
void ScoreSDD::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer)
{
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer); // run the parent version

	ScoreSDD* myWorkerScorer = dynamic_cast<ScoreSDD*>(workerScorer);

	// Absorb various worker thread data
	fTotalEdep += myWorkerScorer->fTotalEdep;
	fNumEvents += myWorkerScorer->fNumEvents;
	fNumProcessHitsCalls += myWorkerScorer->fNumProcessHitsCalls;

	AbsorbDamageBlock(block_pairs, myWorkerScorer->block_pairs);
}



//--------------------------------------------------------------------------------------------------
// This method merges and resolves duplicates of the damage yields from direct and indirect damage.
//--------------------------------------------------------------------------------------------------
std::vector<G4int> ScoreSDD::MergeDamageIndices(std::vector<G4int> &pDamageIndices_direct,
	std::vector<G4int> &pDamageIndices_indirect)
{
	std::vector<G4int> indicesDamage_merged;
	std::vector<G4int> damageIndices_direct_copy = pDamageIndices_direct;
	std::vector<G4int> damageIndices_indirect_copy = pDamageIndices_indirect;

	for (G4int element : damageIndices_direct_copy) {
		if (!IsElementInVector(element, indicesDamage_merged))
			indicesDamage_merged.push_back(element);
		else {
			G4cout << "\tDirect damage site already recorded: " << element << G4endl;
			RemoveElementFromVector(element, pDamageIndices_direct);
		}
	}

	for (G4int element : damageIndices_indirect_copy) {
		if (!IsElementInVector(element, indicesDamage_merged))
			indicesDamage_merged.push_back(element);
		else {
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
void ScoreSDD::ResetMemberVariables() {
	fFiberMapEdepStrand1Backbone.erase(fFiberMapEdepStrand1Backbone.begin(), fFiberMapEdepStrand1Backbone.end());
	fFiberMapEdepStrand2Backbone.erase(fFiberMapEdepStrand2Backbone.begin(), fFiberMapEdepStrand2Backbone.end());
	fFiberMapEdepStrand1Base.erase(fFiberMapEdepStrand1Base.begin(), fFiberMapEdepStrand1Base.end());
	fFiberMapEdepStrand2Base.erase(fFiberMapEdepStrand2Base.begin(), fFiberMapEdepStrand2Base.end());

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
// Record bp indices of one type of simple DNA damage (SSB or BD) in a single strand to a 1D vector.
// This method deletes the content in the provided map.
//--------------------------------------------------------------------------------------------------


std::vector<G4int> ScoreSDD::RecordSimpleDamage(G4double ThreshEDep,
	std::map<G4int, std::pair<G4double, G4ThreeVector>> mapEDep, G4int VoxelID, G4int PositionOnBp)
{
	std::vector<G4int> indicesDamage;

	G4int indexBP;
	G4double eDep;
	G4ThreeVector position;



	// Loop through map of energy depositions. Add indices of energy depositions over the appropriate
	// threshold to a vector, which is returned. The entries in the map are deleted as they are
	// processed
	while ( !mapEDep.empty())
	{
		indexBP = mapEDep.begin()->first;
		eDep = mapEDep.begin()->second.first;
		position = mapEDep.begin()->second.second;

		if (eDep >= ThreshEDep) {
			indicesDamage.push_back(indexBP);
			OutputEvent(0, PositionOnBp, position, VoxelID, indexBP);
		}
		mapEDep.erase(mapEDep.begin());
	
	}

	return indicesDamage;
}

std::vector<G4int> ScoreSDD::RecordBD(std::map<G4int, G4ThreeVector> IndMap , G4int VoxelID, G4int PositionOnBp){

	std::vector<G4int> indicesDamage;

	G4int indexBP;
	G4ThreeVector position;


	while ( !IndMap.empty())
		{
			indexBP = IndMap.begin()->first;
			position = IndMap.begin()->second;
			indicesDamage.push_back(indexBP);
			OutputEvent(1, PositionOnBp, position, VoxelID, indexBP);
			IndMap.erase(IndMap.begin());
		}

		return indicesDamage;
}


//--------------------------------------------------------------------------------------------------
// Record indices of DSBs in a 1D vector. Size should always be even, corresponding to two damage
// sites per DSB. The lowest bp index is recorded first, regardless of whether in strand 1 or 2.
//--------------------------------------------------------------------------------------------------
std::vector<G4int> ScoreSDD::RecordDSB(G4int pDamageCause)
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

	// sort SSBs according to index
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
		// Check damage type of SSB on strand 1 and on strand 2
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
			// Check if DSB is hybrid
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
			// Damage in site 2 is earlier to damage in site 1
			else {
				indicesDSB1D.push_back(*site2);
				indicesDSB1D.push_back(*site1);
			}
			// Remove already-counted damage sites from list of uncounted damage indices
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
// Calculate the order of magnitude (base 10) of a positive integer value.
//--------------------------------------------------------------------------------------------------
G4int ScoreSDD::CalculateIntegerMagnitude(G4int value) {
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
// Helper function to print the contents of a 1D vector.
//--------------------------------------------------------------------------------------------------
void ScoreSDD::Print1DVectorContents(std::vector<G4int> pVector) {
	G4cout << "\t";
	for (G4int element : pVector) {
		G4cout << element << ' ';
	}
	G4cout << G4endl;
}

void ScoreSDD::OutputSDDHeader()
{

	std::ofstream outFile;
	outFile.open(fSDDOutputFile);
	outFile << "SDD Version, SDD" <<  "v1.0"<< ";\n"; //done
	outFile << "Software, TOPAS-nBio;\n"; //done
	outFile << "Author, " << "NICE Team"<< ";\n"; //done
	outFile << "Simulation Details, " << "Neutron Irradiation of a cubic cell nucleus model" << ";\n"; //done
	outFile << "Source, " << "Secondary particles of photoneutrons from a LINAC" << ";\n"; //done
	outFile << "Source type, " << "2" << ";\n"; //done
	outFile << "Incident particles, " << IncidentParticles << ";\n";
	outFile << "Mean particle energy, " << MeanEnergy << ";\n";
	outFile << "Energy distribution, " <<"URL: https://github.com/McGillMedPhys/topas_clustered_dna_damage/tree/master/spectra" << ";\n"; //done
	outFile << "Particle fraction, 1.0;\n"; //done
	outFile << "Dose or fluence,  " << fDoseThreshold/gray << ";\n";
	outFile << "Irradiation target, " << "Full human cubic shape nuclear DNA model containing ~6.3 Gbp constructed using voxels"  << ";\n";//done
	outFile << "Volumes, 1,7.816,7.816,7.816,0,0,0,1,3.901,3.901,3.901,0,0,0;\n"; //done
	outFile << "Chromosome sizes, 46, 263.2,258.1,212.4,203.4,192.2,181.8,168.8,155.5,149.4,144.0,142.9,140.8,121.3,113.0,106.9,94.7,83.5,81.0,68.0,66.6,"
	<< "50.0,52.6,263.2,258.1,212.4,203.4,192.2,181.8,168.8,155.5,149.4,144.0,142.9,140.8,121.3,113.0,106.9,94.7,83.5,81.0,68.0,66.6,50.0,52.6,61.6,165.2"  << ";\n";
	outFile << "DNA Density, " << "13.3" << ";\n"; //done
	outFile << "Cell Cycle Phase, " << "0,1"  << ";\n"; //done
	outFile << "DNA Structure, " << "0,1"  << ";\n"; //done
	outFile << "In vitro / in vivo, " << "0"  << ";\n"; //done
	outFile << "Proliferation status, " << "1"  << ";\n"; //done
	outFile << "Microenvironment, "  << ";\n"; //done
	outFile << "Damage definition, " << "1,0,10,10,17.5" << ";\n"; //done
	outFile << "Time " << "0" << ";\n"; //done
	outFile << "Damage and primary count,  "<< block_pairs.size() << ", " << fNumEvents << ";\n";
	outFile << "Data entries, " << "1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0"<< ";\n";//done
	outFile << "Additional information, " << ";\n"; //done
	outFile << "***EndOfHeader***;\n";
	outFile.close();
}



void ScoreSDD::RecordDamageSDD() {
	// Include following line if want to create a fake, predefined energy map to validate scoring
	//CreateFakeEnergyMap();

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
			fIndicesSSB1_direct = RecordSimpleDamage(fThresEdepForSSB,fFiberMapEdepStrand1Backbone,iVoxel, 1);
			fIndicesSSB2_direct = RecordSimpleDamage(fThresEdepForSSB,fFiberMapEdepStrand2Backbone,iVoxel,4);
			fIndicesBD1_direct = RecordSimpleDamage(fThresEdepForBD,fFiberMapEdepStrand1Base,iVoxel,2);
			fIndicesBD2_direct = RecordSimpleDamage(fThresEdepForBD,fFiberMapEdepStrand2Base,iVoxel,3);


			fIndicesSSB1_indirect = RecordBD(fMapIndDamageStrand1Backbone[iVoxel][iFiber],iVoxel,1);
			fIndicesSSB2_indirect = RecordBD(fMapIndDamageStrand2Backbone[iVoxel][iFiber],iVoxel,4);
			fIndicesBD1_indirect = RecordBD(fMapIndDamageStrand1Base[iVoxel][iFiber],iVoxel,2);
			fIndicesBD2_indirect = RecordBD(fMapIndDamageStrand2Base[iVoxel][iFiber],iVoxel,3);

 			ComputeDSB();

 			if (DsbOnly == false){
	 			fAllSimpleDamage = AllSimpleDamageSDD();
	 			SimpleClusterSDD();
	 		}

 		}
	}
}


void ScoreSDD::ComputeDSB()
{


	std::vector<std::array<G4int,3>> dsbBlock;

	fIndicesSSB1_merged = MergeDamageIndices(fIndicesSSB1_direct, fIndicesSSB1_indirect);
	fIndicesSSB2_merged = MergeDamageIndices(fIndicesSSB2_direct, fIndicesSSB2_indirect);
	fIndicesSSB1 = &fIndicesSSB1_merged;
	fIndicesSSB2 = &fIndicesSSB2_merged;

	


	
	if (fIndicesSSB1->size() == 0 || fIndicesSSB2->size() == 0)
		return;

	// sort SSBs according to index
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

		
		G4bool isSite1Direct = IsElementInVector(*site1, fIndicesSSB1_direct);
		G4bool isSite1Indirect = IsElementInVector(*site1, fIndicesSSB1_indirect);
		G4bool isSite2Direct = IsElementInVector(*site2, fIndicesSSB2_direct);
		G4bool isSite2Indirect = IsElementInVector(*site2, fIndicesSSB2_indirect);

	
		if (isSite1Direct && isSite1Indirect) {
			G4cerr << "Error: (While scoring DSB) CONTRADICTION! site1: " << *site1 << " is both direct and indirect!" << G4endl;
			exit(0);
		}
		if (isSite2Direct && isSite2Indirect) {
			G4cerr << "Error: (While scoring DSB) CONTRADICTION! site2: " << *site2 << " is both direct and indirect!" << G4endl;
			exit(0);
			}



		// Damage in site 2 is within range of site 1 to count as DSB (either before or after)
		if (isDSB && isDSBrecorded) {

				if (isSite1Direct) {
					dsbBlock.push_back({*site1,fIdBckBone1,fIdDirect});
					RemoveElementFromVector(*site1, fIndicesSSB1_direct);
				}
				else {
					dsbBlock.push_back({*site1,fIdBckBone1,fIdIndirect});
					RemoveElementFromVector(*site1, fIndicesSSB1_indirect);
				}


				if (isSite2Direct) {
					dsbBlock.push_back({*site2,fIdBckBone2,fIdDirect});
					RemoveElementFromVector(*site2, fIndicesSSB2_direct);
				}
				else {
					dsbBlock.push_back({*site2,fIdBckBone2,fIdIndirect});
					RemoveElementFromVector(*site2, fIndicesSSB2_indirect);
				}

			min_site = std::min(*site1, *site2);
			max_site = std::max(*site1, *site2);


			if (*site1 < *site2){
    		auto nextsite = std::next(site1);
    		while (nextsite != fIndicesSSB1->end() && *nextsite <= *site2){
      		  if (IsElementInVector(*nextsite, fIndicesSSB1_direct)){
        		    dsbBlock.push_back({*nextsite, fIdBckBone1, fIdDirect});
         		   RemoveElementFromVector(*nextsite, fIndicesSSB1_direct);
        		}
        		else if (IsElementInVector(*nextsite, fIndicesSSB1_indirect)){ 
           			 dsbBlock.push_back({*nextsite, fIdBckBone1, fIdIndirect});
           			 RemoveElementFromVector(*nextsite, fIndicesSSB1_indirect);
        		}
        		else{
            		G4cout << "!!! Issue, the damage isn't contained in any list !!! " << G4endl;
        		}
        		nextsite = fIndicesSSB1->erase(nextsite);
   			 }
			}

			if (*site2 < *site1){
				auto nextsite = std::next(site2);
				while (nextsite != fIndicesSSB2->end() && *nextsite <= *site1){
					if (IsElementInVector(*nextsite, fIndicesSSB2_direct)){
						dsbBlock.push_back({*nextsite, fIdBckBone2, fIdDirect});
						RemoveElementFromVector(*nextsite, fIndicesSSB2_direct);
					}
					else if (IsElementInVector(*nextsite, fIndicesSSB2_indirect)){ 
						dsbBlock.push_back({*nextsite, fIdBckBone2, fIdIndirect});
						RemoveElementFromVector(*nextsite, fIndicesSSB2_indirect);
					}
					else{
						G4cout << "!!! Issue, the damage isn't contained in any list !!! " << G4endl;
					}
					nextsite = fIndicesSSB2->erase(nextsite);
				}
			}

			site1 = fIndicesSSB1->erase(site1);
			site2 = fIndicesSSB2->erase(site2);



			InBlock_direct_bd1 = FindBaseDamage(fIndicesBD1_direct, min_site, max_site);
			InBlock_indirect_bd1 = FindBaseDamage(fIndicesBD1_indirect, min_site, max_site);
			InBlock_direct_bd2 = FindBaseDamage(fIndicesBD2_direct, min_site, max_site);
			InBlock_indirect_bd2 = FindBaseDamage(fIndicesBD2_indirect, min_site, max_site);

			for (auto bp : InBlock_direct_bd1) {
        		dsbBlock.push_back({bp, fIdBase1, fIdDirect});
        		RemoveElementFromVector(bp, fIndicesBD1_direct);
    		}

    		for (auto bp : InBlock_indirect_bd1) {
        		dsbBlock.push_back({bp, fIdBase1, fIdIndirect});
        		RemoveElementFromVector(bp, fIndicesBD1_indirect);
    		}

    		for (auto bp : InBlock_direct_bd2) {
        		dsbBlock.push_back({bp, fIdBase2, fIdDirect});
        		RemoveElementFromVector(bp, fIndicesBD2_direct);
    		}

    		for (auto bp : InBlock_indirect_bd2) {
        		dsbBlock.push_back({bp, fIdBase2, fIdIndirect});
        		RemoveElementFromVector(bp, fIndicesBD2_indirect);
    		}


    		RecordBlock(fEventID,dsbBlock,1);
			dsbBlock.clear();

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
}

	


std::vector<G4int> ScoreSDD::FindBaseDamage(std::vector<G4int> indices, G4int MinBp, G4int MaxBp) {

    std::vector<G4int> BaseDamage;

    for (auto BpID : indices) {
        if (BpID >= MinBp && BpID <= MaxBp) {
            BaseDamage.push_back(BpID);
       }

    }
    return BaseDamage;
}

std::vector<std::array<G4int,3>> ScoreSDD::AllSimpleDamageSDD() {
	// G4cout << "COMBINING SIMPLE DAMAGE" << G4endl;
	std::vector<std::array<G4int,3>> AllSimpleDamage;

	// Fake data for testing
	// fIndicesSSB1 = {3,4,5,7,9};
	// fIndicesBD1 = {1,4,7};
	// fIndicesSSB2 = {1,12};
	// fIndicesBD2 = {6,7,10,12};
	// fIndicesDSB = {3,5,11,14};

	// SSBs in strand 1 (direct)
	std::vector<G4int>::iterator iter = fIndicesSSB1_direct.begin();
	while (iter != fIndicesSSB1_direct.end()) {
		AllSimpleDamage.push_back({*iter,fIdBckBone1,fIdDirect});
		iter++;
	}

	// SSBs in strand 1 (indirect)
	iter = fIndicesSSB1_indirect.begin();
	while (iter != fIndicesSSB1_indirect.end()) {
		AllSimpleDamage.push_back({*iter,fIdBckBone1,fIdIndirect});
		iter++;
	}

	// BDs in strand 1 (direct)
	iter = fIndicesBD1_direct.begin();
	while (iter != fIndicesBD1_direct.end()) {
		AllSimpleDamage.push_back({*iter,fIdBase1,fIdDirect});
		iter++;
	}

	//Direct-Indirect duplicate not dealt with for bases yet
	// BDs in strand 1 (indirect)
	iter = fIndicesBD1_indirect.begin();
	while (iter != fIndicesBD1_indirect.end()) {
		if(!IsElementInVector(*iter, fIndicesBD1_direct)){
			AllSimpleDamage.push_back({*iter,fIdBase1,fIdIndirect});
		}
		iter++;
	}

	// SSBs in strand 2 (direct)
	iter = fIndicesSSB2_direct.begin();
	while (iter != fIndicesSSB2_direct.end()) {
		AllSimpleDamage.push_back({*iter,fIdBckBone2,fIdDirect});
		iter++;
	}

	// SSBs in strand 2 (indirect)
	iter = fIndicesSSB2_indirect.begin();
	while (iter != fIndicesSSB2_indirect.end()) {
		AllSimpleDamage.push_back({*iter,fIdBckBone2,fIdIndirect});
		iter++;
	}

	// BDs in strand 2 (direct)
	iter = fIndicesBD2_direct.begin();
	while (iter != fIndicesBD2_direct.end()) {
		AllSimpleDamage.push_back({*iter,fIdBase2,fIdDirect});
		iter++;
	}

	// BDs in strand 2 (indirect)
	iter = fIndicesBD2_indirect.begin();
	while (iter != fIndicesBD2_indirect.end()) {
		if(!IsElementInVector(*iter, fIndicesBD2_direct)){
			AllSimpleDamage.push_back({*iter,fIdBase2,fIdIndirect});
		}
		iter++;
	}

	if (AllSimpleDamage.size() > 1)
		std::sort(AllSimpleDamage.begin(), AllSimpleDamage.end()); // sorts by first element in each vector item by default (i.e site index)

	return AllSimpleDamage;
}

void ScoreSDD::SimpleClusterSDD()
{
	// Only 1 or 0 damages, so no clustering.
	if (fAllSimpleDamage.size() == 0){
		return;
	}

	if (fAllSimpleDamage.size() == 1){


		//std::vector<std::array<G4int,3>>::iterator damage_pointer = fAllSimpleDamage.begin();
		//std::array<G4int,3> damage = *damage_pointer;

		std::array<G4int, 3> damage = *(fAllSimpleDamage.begin());

		RecordBlock(fEventID,fAllSimpleDamage,0);
		return;
	}

	std::vector<std::array<G4int,3>>::iterator site = fAllSimpleDamage.begin();

	std::vector<std::array<G4int,3>> CurrentCluster = {*site};

	while (site != fAllSimpleDamage.end()) {
		std::array<G4int,3> NextSite = *(site+1); 

		if ((NextSite[0] - CurrentCluster[0][0]) <= 10) {
			CurrentCluster.push_back(NextSite);
		}

		else {
			RecordBlock(fEventID,CurrentCluster,0);
			CurrentCluster = {NextSite};
		}
		site++;

		if ((site + 1) == fAllSimpleDamage.end()) {
			RecordBlock(fEventID,CurrentCluster,0);
			site++;  
		}
	}
}


void ScoreSDD::RecordBlock(G4int EventID, std::vector<std::array<G4int,3>> block, G4int DSBcount){

	//Field2
	G4String field2;
	if (block.size() == 1){	
		std::array<G4int,3> only_entry = *(block.begin());
		G4ThreeVector coord_only_entry = GetXYZ(only_entry[0], only_entry[1], only_entry[2]);
		field2 = std::to_string(coord_only_entry[0]/um) + ", " + std::to_string(coord_only_entry[1]/um) + ", " + std::to_string(coord_only_entry[2]/um);
	}

	else{

		std::sort(block.begin(), block.end());

		std::array<G4int,3> first_entry = block[0];
		std::array<G4int,3> last_entry = block[block.size()-1];

		G4ThreeVector coord_first_entry = GetXYZ(first_entry[0], first_entry[1], first_entry[2]);
		G4ThreeVector coord_last_entry = GetXYZ(last_entry[0], last_entry[1], last_entry[2]);

		G4ThreeVector box_middle = ((coord_first_entry + coord_last_entry)/2)/um;

		G4double x_max = std::max(coord_first_entry[0], coord_last_entry[0])/um;
		G4double x_min = std::min(coord_first_entry[0], coord_last_entry[0])/um;

		G4double y_max = std::max(coord_first_entry[1], coord_last_entry[1])/um;
		G4double y_min = std::min(coord_first_entry[1], coord_last_entry[1])/um;

		G4double z_max = std::max(coord_first_entry[2], coord_last_entry[2])/um;
		G4double z_min = std::min(coord_first_entry[2], coord_last_entry[2])/um;

		G4String middle_string = std::to_string(box_middle[0]) + ", " + std::to_string(box_middle[1]) + ", " + std::to_string(box_middle[2]);
		G4String max_string = std::to_string(x_max) + ", " + std::to_string(y_max) + ", " + std::to_string(z_max);
		G4String min_string = std::to_string(x_min) + ", " + std::to_string(y_min) + ", " + std::to_string(z_min);

		field2 = middle_string + " / " + max_string + " / " + min_string;
	}
	//Field 3 
	G4int Chromosome = ChromoMapping[fVoxelID][0];

	G4int arm = ChromoMapping[fVoxelID][2];

	G4String field3 = "0, " + std::to_string(Chromosome) + ", 1, " + std::to_string(arm); 

	//Field 4  
	G4int BpPerFiber = 18000;
	G4int BpPerVoxel = 20*BpPerFiber;

	G4int PositionInBp = ChromoMapping[fVoxelID][1]*BpPerVoxel + fFiberID*BpPerFiber + block[0][0];
	G4String field4 = std::to_string(PositionInBp);

	//Field5
	G4int BlockMode = 0;
	G4int DirectCount = 0;
	G4int IndirectCount = 0;

	//Field6
	G4int StrandCount = 0;
	G4int BaseCount = 0;


	std::vector<std::array<G4int,3>> FullBreak;
	G4int FirstBP = block[0][0];

	G4String field7 = std::to_string((block[0][0]-FirstBP+1)) +", "+ std::to_string(block[0][1]) +", "+ std::to_string(block[0][2]+1);

	//Field 6&7
	for (int i = 0; i < block.size(); i++) {

		//Field 5
		if (block[i][2] == 0){
			DirectCount++;
			}

		else{
			IndirectCount++;
			}
		

		//Field6
		if (block[i][1] == 1 || block[i][1] == 4){
			StrandCount++;
			}

		else{
			BaseCount++;
			}

		if (i > 0){
			//Field 7
			field7 += " / " + std::to_string((block[i][0]-FirstBP+1)) +", "+ std::to_string(block[i][1]) +", "+ std::to_string(block[i][2]+1);
			}
		}


	if (IndirectCount > 0 && DirectCount == 0){
		BlockMode += 1;
	}

	if (IndirectCount > 0 && DirectCount > 0){
		BlockMode += 2;
	}


	G4String field5 = std::to_string(BlockMode) + ", " + std::to_string(DirectCount) + ", " + std::to_string(IndirectCount);

	
	G4String field6 = std::to_string(BaseCount) + ", " + std::to_string(StrandCount) + ", " + std::to_string(DSBcount);


	G4String block_string = field2+"; "+field3+"; "+field4+"; "+field5+"; "+field6+"; "+field7+";";

	block_pairs.push_back(std::make_pair(EventID, block_string));

	/*std::ofstream outFile;
		outFile.open(fSDDOutputFile, std::ios::app);
		outFile << EventID << ";" << field2 << "\n";
		outFile.close();*/


}


G4ThreeVector ScoreSDD::GetXYZ(G4int BpPosition, G4int PositionOnBp, G4int Mode){

	G4ThreeVector Coord;

	if (PositionOnBp == 1 && Mode == 0 ){
		Coord = fMapEdepStrand1Backbone[fVoxelID][fFiberID][BpPosition].second;
		}

	if (PositionOnBp == 4 && Mode == 0 ){
		Coord = fMapEdepStrand2Backbone[fVoxelID][fFiberID][BpPosition].second;
		}

	if (PositionOnBp == 2 && Mode == 0 ){
		Coord = fMapEdepStrand1Base[fVoxelID][fFiberID][BpPosition].second;
		}

	if (PositionOnBp == 3 && Mode == 0 ){
		Coord = fMapEdepStrand2Base[fVoxelID][fFiberID][BpPosition].second;
		}

	if (PositionOnBp == 1 && Mode == 1 ){
		Coord = fMapIndDamageStrand1Backbone[fVoxelID][fFiberID][BpPosition];
		}

	if (PositionOnBp == 4 && Mode == 1 ){
		Coord = fMapIndDamageStrand2Backbone[fVoxelID][fFiberID][BpPosition];
		}

	if (PositionOnBp == 2 && Mode == 1 ){
		Coord = fMapIndDamageStrand1Base[fVoxelID][fFiberID][BpPosition];
		}

	if (PositionOnBp == 3 && Mode == 1 ){
		Coord = fMapIndDamageStrand2Base[fVoxelID][fFiberID][BpPosition];
		}

	return Coord;
}


void ScoreSDD::PrintSDDdata(){

	std::sort(block_pairs.begin(), block_pairs.end());

	if(block_pairs.size() == 0){
		return;
	}

	else {

		G4int EventID = block_pairs[0].first;
	    std::string block_string = block_pairs[0].second;
	    G4int NewEvent = 1;

	    std::ofstream outFile;
		outFile.open(fSDDOutputFile, std::ios::app);
		outFile << NewEvent << ", " << EventID << "; " << block_string << "\n";
		outFile.close();

	
		if (block_pairs.size() > 1){
			for (G4int i = 1; i < block_pairs.size(); i++) {
							
			EventID = block_pairs[i].first;
		    block_string = block_pairs[i].second;

		    if (EventID  == block_pairs[i-1].first){
		    	NewEvent = 0;
		    	}

		    else{
		    	NewEvent = 1;
		    }

		    std::ofstream outFile;
			outFile.open(fSDDOutputFile, std::ios::app);
			outFile << NewEvent << ", " << EventID << "; " << block_string << "\n";
			outFile.close();

			}

		}
	}

}


void ScoreSDD::AbsorbDamageBlock(
	std::vector<std::pair<int, std::string>> &masterMap,
	std::vector<std::pair<int, std::string>> &workerMap)
{
	for (G4int i = 0; i < workerMap.size(); i++){
		masterMap.push_back(workerMap[i]);

	}

}

G4int ScoreSDD::GetChromoMap(){

    std::ifstream inFile;

    inFile.open(fMapFile);

    if(inFile.fail()){
        std::cout << "ERROR: Could not read  " << fMapFile << std::endl;
        return 1;
    }

    G4int VoxelID, ChromoID, BpCumul, Arm;
    std::vector<G4int> values;

for (G4int i = 0; i < 17576; i++) {
    inFile >> VoxelID >> ChromoID >> BpCumul >> Arm;
    //std::cout << x << ",  " << y << ",  " << z << std::endl;
    values.clear();
    values.push_back(ChromoID);
    values.push_back(BpCumul);
    values.push_back(Arm);
    ChromoMapping[VoxelID] = values;
}

	inFile.close();
//testing
 //std::cout << "key: 27 ChromoID: " << ChromoMapping[27][0] << " BpCumul: " << ChromoMapping[27][1] << std::endl;

    return 0;
}





