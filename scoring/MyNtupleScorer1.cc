// Scorer for MyNtupleScorer1
//
// ********************************************************************
// *                                                                  *
// *                                                                  *
// * This file was obtained from Topas MC Inc under the license       *
// * agreement set forth at http://www.topasmc.org/registration       *
// * Any use of this file constitutes full acceptance of              *
// * this TOPAS MC license agreement.                                 *
// *                                                                  *
// ********************************************************************
//

#include "MyNtupleScorer1.hh"

#include "G4PSDirectionFlag.hh"
#include "G4VProcess.hh"

MyNtupleScorer1::MyNtupleScorer1(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
						  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
						 : TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	// SetSurfaceScorer();

	// fNtuple->RegisterColumnD(&fEnergy, "Energy", "MeV");
	// fNtuple->RegisterColumnF(&fWeight, "Weight", "");
	// fNtuple->RegisterColumnI(&fParticleType, "Particle Type (in PDG Format)");
	// fNtuple->RegisterColumnB(&fIsNewHistory, "Flag to tell if this is the First Scored Particle from this History (1 means true)");

	// fComponent->GetMotherLogical()->GetInstanceID();

	G4String fiberMaterialName = "G4_WATER";
	if ( fPm->ParameterExists(GetFullParmName("FiberMaterialName")))
		fiberMaterialName = fPm->GetStringParameter(GetFullParmName("FiberMaterialName"));

	fFiberMaterial = GetMaterial(fiberMaterialName);

	fNtuple->RegisterColumnI(&fEventID, "History ID");
	fNtuple->RegisterColumnI(&fTrackID, "Track ID");
	fNtuple->RegisterColumnD(&fEnergy, "Energy Deposited", "eV");
	fNtuple->RegisterColumnS(&fVolName, "Physical Volume");
	fNtuple->RegisterColumnS(&fProcess, "Process");
}


MyNtupleScorer1::~MyNtupleScorer1() {;}


G4bool MyNtupleScorer1::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}

	ResolveSolid(aStep); // Get the solid in which the current step resides (if using area, volume, or material)


	fEnergy = aStep->GetTotalEnergyDeposit();
	if ( fEnergy > 0. ) {
		G4Material* material = aStep->GetPostStepPoint()->GetMaterial(); // Note: get erroneus interactions in the nucleus if use PreStepPoint

		if ( material == fFiberMaterial) {
			fEventID = GetEventID();	
			fTrackID = aStep->GetTrack()->GetTrackID();	
			fVolName = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
			fProcess = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

			fNtuple->Fill();
			return true;
		}
		else {
			return false;
		}
	}
	return false;	


	// if (IsSelectedSurface(aStep)) { 	// Is there a hit on the surface of the solid
	// 	G4StepPoint* theStepPoint=0;
	// 	G4int direction = GetDirection(); // entering / exiting volume
	// 	if (direction == fFlux_In) theStepPoint = aStep->GetPreStepPoint();
	// 	else if (direction == fFlux_Out) theStepPoint = aStep->GetPostStepPoint();
	// 	else return false;

	// 	fEnergy       = theStepPoint->GetKineticEnergy();
	// 	fWeight       = theStepPoint->GetWeight();
	// 	fParticleType = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

	// 	// Check if this is a new history
	// 	fRunID   = GetRunID();
	// 	fEventID = GetEventID();
	// 	if (fEventID != fPrevEventID || fRunID != fPrevRunID) {
	// 		fIsNewHistory = true;
	// 		fPrevEventID = fEventID;
	// 		fPrevRunID = fRunID;
	// 	} else {
	// 		fIsNewHistory = false;
	// 	}

	// 	fNtuple->Fill();
	// 	return true;
	// }
	// return false;
}
