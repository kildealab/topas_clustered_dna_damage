// Scorer for clustereddnadamage
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//

#include "ScoreClusteredDNADamage.hh"
#include "TsTrackInformation.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <map>

//--------------------------------------------------------------------------------------------------
// Constructor Must define each data column & type to be recorded. In our case, record data across
// four columns:
// (1) Primary particle number (i.e. event number)
// (2) DNA fibre number (only relevant if multiple fibres used)
// (3) Total number of SSBs (SBs contributing to a DSB are not counted)
// (4) Total number of DSBs
//--------------------------------------------------------------------------------------------------
ScoreClusteredDNADamage::ScoreClusteredDNADamage(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
											   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
:TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)

{
	SetUnit("");
	
	// Default parameters
	fThresDistForDSB = 10;
	fThresEdepForSSB = 17.5 * eV;
	fThresEdepForBD = 17.5 * eV;

	fNumEdeps1 = 0;
	fTotalEdep1 = 0.;
	fNumEdeps2 = 0;
	fTotalEdep2 = 0.;
	fNumEdepsBD1 = 0;
	fTotalEdepBD1 = 0.;
	fNumEdepsBD2 = 0;
	fTotalEdepBD2 = 0.;

	
	fNbOfAlgo = 1; // Used with variance reduction
	
	if ( fPm->ParameterExists(GetFullParmName("BasePairDistanceForDefiningDSB")) )
		fThresDistForDSB = fPm->GetIntegerParameter(GetFullParmName("BasePairDistanceForDefiningDSB"));
	if ( fPm->ParameterExists(GetFullParmName("EnergyThresholdForHavingSSB")) )
		fThresEdepForSSB = fPm->GetDoubleParameter(GetFullParmName("EnergyThresholdForHavingSSB"), "Energy");
	if ( fPm->ParameterExists(GetFullParmName("EnergyThresholdForHavingBD")) )
		fThresEdepForSSB = fPm->GetDoubleParameter(GetFullParmName("EnergyThresholdToInduceBD"), "Energy");
	
	// Unsure what this is
	fBasePairDepth = 0;
	if ( fPm->ParameterExists(GetFullParmName("BasePairPositionAtGeometricHierarchy")))
		fBasePairDepth = fPm->GetIntegerParameter(GetFullParmName("BasePairPositionAtGeometricHierarchy"));
	
	// Material of DNA residue volumes in which to score
	G4String DNAMaterialName = "G4_WATER";
	if ( fPm->ParameterExists(GetFullParmName("DNAMaterialName")))
		DNAMaterialName = fPm->GetStringParameter(GetFullParmName("DNAMaterialName"));
	fDNAMaterial = GetMaterial(DNAMaterialName);

	// Determine order of magnitude of number of base pairs. Set fParser parameters accordingly
	// for use in ProcessHits()
	fNucleoNum = 0;
	fBpNum = 0;
	if ( fPm->ParameterExists(GetFullParmName("DnaNumNucleosomePerFiber")) )
		fNucleoNum = fPm->GetIntegerParameter(GetFullParmName("DnaNumNucleosomePerFiber"));
	if ( fPm->ParameterExists(GetFullParmName("DnaNumBpPerNucleosome")) )
		fBpNum = fPm->GetIntegerParameter(GetFullParmName("DnaNumBpPerNucleosome"));
	fNumBpMagnitude = CalculateIntegerMagnitude(fNucleoNum*fBpNum);
	fParserResidue = fNumBpMagnitude*10;
	fParserStrand = fNumBpMagnitude*100;
	
	// This is for variance reduction --> Not used with CharltonDNA
	if ( fPm->ParameterExists(GetFullParmName("NumberOfSplit")) )
		fNbOfAlgo = fPm->GetIntegerParameter(GetFullParmName("NumberOfSplit"));
	
	//fThresEdepForSSB /= eV;
	
	fNtuple->RegisterColumnI(&fEventID, "Event number"); // ID of primary particle / event / history
	fNtuple->RegisterColumnI(&fDNAParent, "DNA parent geometry"); // Unsure. Not initiated in header either
	fNtuple->RegisterColumnI(&fSSB,       "Single strand breaks"); // Number of SSB caused by this primary particle
	fNtuple->RegisterColumnI(&fDSB,       "Double strand breaks"); // Number of DSB caused by this primary particle
	fNtuple->RegisterColumnI(&fBD,       "Base damages"); // Number of DSB caused by this primary particle
	
	// Unsure what this does
	// From the docs: disable automatic creation & filling of output, leaving this work entirely to
	// your scorer.
	// Commenting/uncommenting seems to have no impact on CharltonDNA simulation.
	SuppressStandardOutputHandling();
}


//--------------------------------------------------------------------------------------------------
// Destructor
//--------------------------------------------------------------------------------------------------
ScoreClusteredDNADamage::~ScoreClusteredDNADamage() {
}

//--------------------------------------------------------------------------------------------------
// ProcessHits method. This method is called for every hit in the sensitive detector volume, i.e.
// when an interaction occurs in a sensitive volume (may or may not be energy deposit). In this case
// simply record energy deposited in a given bp index, in 1 of 2 strands. Multiple energy
// depositions occurring in the same volume are added together and recorded as one.
//
// Note about material filtering
// 
// From TsVScorer.hh:
// Code in this method must be written as efficiently as possible. Do not directly access parameters
// from here. Instead, access and cache them in the method UpdateForParameterChange described above.
// Do not do string comparisons if you can help it. They tend to be slow
//--------------------------------------------------------------------------------------------------
G4bool ScoreClusteredDNADamage::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	// Something to do with "gated scoring" as per TsVScorer.hh. Unused in CharltonDNA example
	if (!fIsActive) {
		G4cout << "Gated scoring!" << G4endl;
		fSkippedWhileInactive++;
		return false;
	}
	
	G4double edep = aStep->GetTotalEnergyDeposit(); ///eV;

	//----------------------------------------------------------------------------------------------
	// If this hit deposits energy: 
	//----------------------------------------------------------------------------------------------
	if (edep > 0) {
		G4StepPoint* preStep = aStep->GetPreStepPoint();

		// Get the indices defining the volume in which energy was deposited
		// Can try making these member variables --> check performance
		G4int num_strand = -1;
		G4int num_res = -1;
		G4int num_nucleotide = -1;

		G4int volID = preStep->GetPhysicalVolume()->GetCopyNo();

		// num_strand = volID / fParserStrand;
		// num_res = (volID - (num_strand*fParserStrand)) / fParserResidue;
		// num_nucleotide = volID - (num_strand*fParserStrand) - (num_res*fParserResidue);
		num_strand = volID / 1000000;
		num_res = (volID - (num_strand*1000000)) / 100000;
		num_nucleotide = volID - (num_strand*1000000) - (num_res*100000);


		// if (num_strand == 0) {
		// 	fNumEdeps1++;
		// 	fTotalEdep1 += edep;
		// }
		// else if (num_strand == 1) {
		// 	fNumEdeps2++;
		// 	fTotalEdep2 += edep;
		// }
		if (num_strand == 0) {
			if (num_res == 0 || num_res == 1) {
				fNumEdeps1++;
				fTotalEdep1 += edep;
			}
			else if (num_res == 2) {
				fNumEdepsBD1++;
				fTotalEdepBD1 += edep;
			}
		}
		else if (num_strand == 1) {
			if (num_res == 0 || num_res == 1) {
				fNumEdeps2++;
				fTotalEdep2 += edep;
			}
			else if (num_res == 2) {
				fNumEdepsBD2++;
				fTotalEdepBD2 += edep;
			}
		}

		if ( num_strand == 0 || num_strand == 1 ) {
			G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());
			// G4int numStrand  = touchable->GetVolume(fBasePairDepth)->GetCopyNo(); // This is the bp index (e.g. 184)
			G4int parentDepth = touchable->GetVolume(fBasePairDepth+1)->GetCopyNo(); // Always seems to be zero. May change with multiple DNA fibres -> i.e. index of DNA fibre?

			//--------------------------------------------------------------------------------------
			// Record energy deposited in fGenVEdepStrand (3-nested map) for either strand 1 or 2,
			// adding to existing deposits in that nucleotide (if any)
			// First index specifies the DNA fibre
			// Second index specifies someting to do with variance reduction / track splitting
			// Third index specifies the bp index
			//--------------------------------------------------------------------------------------
			// Parameters needed to handle variance reduction track splitting, if used or not
			G4int index = 1;
			if (fNbOfAlgo > 1) {
				TsTrackInformation* trackInformation = (TsTrackInformation*)aStep->GetTrack()->GetUserInformation();
				index = trackInformation->GetSplitTrackID();
			}

			// If doing variance reduction
			if (index > 2) {
				if ( num_strand == 0 ){
					fGenVEdepStrand1[parentDepth][index-3][num_nucleotide] += edep;
					// G4cout << "0 - " << num_nucleotide << G4endl;
				}
				else{
					fGenVEdepStrand2[parentDepth][index-3][num_nucleotide] += edep;
					// G4cout << "1 - " << num_nucleotide << G4endl;
				}
			} 
			// Not doing variance reduction
			else {
				if ( num_strand == 0 ){
					for ( int i = 0; i < fNbOfAlgo; i++ ){
						fGenVEdepStrand1[parentDepth][i][num_nucleotide] += edep;
						// G4cout << "0 - " << num_nucleotide << G4endl;
					}
				}
				else{
					for ( int i = 0; i < fNbOfAlgo; i++ ){
						fGenVEdepStrand2[parentDepth][i][num_nucleotide] += edep;
						// G4cout << "1 - " << num_nucleotide << G4endl;
					}
				}
			}
			
			return true;
		}

	}
	return false;
}

//--------------------------------------------------------------------------------------------------
// This method is called at the end of the event (primary particle & its secondaries). Determine the 
// number of SSB and DSB caused by the current event & record in ntuple for output.
//
// **Might be a good idea to test how long it takes to compute # of strand breaks
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::UserHookForEndOfEvent() {
	fEventID = GetEventID();

	// if (fEventID == 499) {
	// 	G4cout << "#################################################################################################################################" << G4endl;
	// 	// if (fNumEdeps1 > 0 ){
	// 	// 	G4cout << "Number of edeps in strand 1: " << fNumEdeps1 << G4endl;
	// 	// 	G4cout << "Total energy deposited: " << fTotalEdep1/eV << G4endl;
	// 	// }
	// 	// if (fNumEdeps2 > 0){
	// 	// 	G4cout << "Number of edeps in strand 2: " << fNumEdeps2 << G4endl;
	// 	// 	G4cout << "Total energy deposited: " << fTotalEdep2/eV << G4endl;
	// 	// }
	// 	G4cout << "Number of edeps in strand 1 backbone: " << fNumEdeps1 << G4endl;
	// 	G4cout << "Total energy deposited in strand 1 backbone: " << fTotalEdep1/eV << G4endl;
	// 	G4cout << "Number of edeps in strand 1 bases: " << fNumEdepsBD1 << G4endl;
	// 	G4cout << "Total energy deposited in strand 1 bases: " << fTotalEdepBD1/eV << G4endl;
	// 	G4cout << "Number of edeps in strand 2: " << fNumEdeps2 << G4endl;
	// 	G4cout << "Total energy deposited: " << fTotalEdep2/eV << G4endl;
	// 	G4cout << "Number of edeps in strand 2 bases: " << fNumEdepsBD2 << G4endl;
	// 	G4cout << "Total energy deposited in strand 2 bases: " << fTotalEdepBD2/eV << G4endl;
	// }


	// Include if want to print out map of energy depositions to validate algorithm is functioning
	// properly.
	// std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand1_COPY = fGenVEdepStrand1;
	// std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand2_COPY = fGenVEdepStrand2;

	// Iterate over all DNA fibres 
	for ( auto& energyAtStrands : fGenVEdepStrand1 ) {
		G4int parentDepth = energyAtStrands.first;

		// Unsure why this necessary. Perhaps is faster to do this rather than repeatedly indexing
		// the extra layer of mappng associated with fGenVEdepStrand1/2
		fVEdepStrand1 = fGenVEdepStrand1[parentDepth];
		fVEdepStrand2 = fGenVEdepStrand2[parentDepth];
		fDNAParent = parentDepth;
		
		// Iterate over track splitting (only 1 iteration if no splitting) and compute the # of SSB
		// and DSB induced by the current event. Fill the nTuple with these data.
		for ( int i = 0; i < fNbOfAlgo; i++ ) {
			G4int sb[2] = {0, 0}; // 2-element array to store the # of SSB and DSB
			ComputeStrandBreaks(sb, i);

			// Record SSB & DSB to appropriate variables, required to fill ntuple output
			fSSB = sb[0];
			fDSB = sb[1];
			if ( fSSB > 0 || fDSB > 0 )
				fNtuple->Fill();
		}

		// Include if want to print out map of energy depositions to validate algorithm is
		// functioning properly.
		// if (fSSB > 0 && fDSB > 0) {
		// 	G4cout << "#################################################################################################################################" << G4endl;
		// 	G4cout << "Event #" << fEventID << G4endl;
		// 	for (G4int i = 0; i < 18000; i++) {
		// 		G4cout << fGenVEdepStrand1_COPY[0][0][i]/eV << "         " << fGenVEdepStrand2_COPY[0][0][i]/eV << G4endl;
		// 	}
		// 	G4cout << "#################################################################################################################################" << G4endl;
		// }

		// Clear member variables for next fibre
		fVEdepStrand1.erase(fVEdepStrand1.begin(), fVEdepStrand1.end());
		fVEdepStrand2.erase(fVEdepStrand2.begin(), fVEdepStrand2.end());
	}
	// Clear member variables for next event
	fGenVEdepStrand1.erase(fGenVEdepStrand1.begin(), fGenVEdepStrand1.end());
	fGenVEdepStrand2.erase(fGenVEdepStrand2.begin(), fGenVEdepStrand2.end());
}


//--------------------------------------------------------------------------------------------------
// Process map of energy depositions (total energy deposited per nucleotide) spanning two strands of
// a single DNA fibre, to calculate the  number of SSB and DSB present. The number of SSBs in each
// strand are recorded separately but are summed at the end and recorded in sb. The number of DSB is
// also recorded in sb. 
//
// This class was taken from Geant4/examples/extended/medical/dna/pdb4dna
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::ComputeStrandBreaks(G4int* sb, G4int cluster)
{
	// sb quantities
	G4int ssb1=0;
	G4int ssb2=0;
	G4int dsb=0;
	
	// nucleotide id and energy deposit for each strand
	G4int nucl1;
	G4int nucl2;
	G4double edep1;
	G4double edep2;
	
	//----------------------------------------------------------------------------------------------
	// Read through energy depositions in the first DNA strand one at a time. For each, determine
	// if is SSB1, then examine the second DNA strand. Increment through depositions in the second
	// strand, recording SSB2 where applicable, and stopping when the nucleotide index in second
	// strand is either (A) within DSB-distance-threshold of SSB1 or (B) greater than that of SSB1.
	// Once either of those criteria is met, determine if the current SSB1 & SSB2 are close enough
	// to constitute a DSB. If so, record as such. If not, in the case of (B), SSB2 is not recorded
	// and is left in the queue to be processed with next SSB1 as potential DSB. Then exit iterating
	// through second strand and return to first strand.
	//
	// Note that iterators through maps proceed in order of nucleotide index, which is required
	// for the logic of this algorithm.
	//----------------------------------------------------------------------------------------------
	while ( !fVEdepStrand1[cluster].empty() )
	{
		// Extract nucleotide index & energy deposit, then delete the entry from the map
		nucl1 = fVEdepStrand1[cluster].begin()->first;
		edep1 = fVEdepStrand1[cluster].begin()->second;
		fVEdepStrand1[cluster].erase( fVEdepStrand1[cluster].begin() );

		// Temporarily record as SSB in strand1. May remove later if DSB
		if ( edep1 >= fThresEdepForSSB )
		{
			ssb1++;
		}
		
		// Look at strand2, record all SSB2 up until same index as SSB1 and potential DSB within
		// plus-minus threshold of SSB1
		if ( !fVEdepStrand2[cluster].empty() )
		{
			do
			{
				nucl2 = fVEdepStrand2[cluster].begin()->first;
				edep2 = fVEdepStrand2[cluster].begin()->second;

				// SSB in strand2
				if ( edep2 >= fThresEdepForSSB )
				{
					ssb2++;
				}
				fVEdepStrand2[cluster].erase( fVEdepStrand2[cluster].begin() );
			} while ( ((nucl1-nucl2)>fThresDistForDSB) && (!fVEdepStrand2[cluster].empty()) );
			
			// This is the scenario described above, wherein index of SSB2 exceeds SSB1 by more than
			// DSB threshold, so return SSB2 to queue for processing with next SSB1
			if ( nucl2-nucl1 > fThresDistForDSB )
			{
				fVEdepStrand2[cluster][nucl2]=edep2;
				if ( edep2 >= fThresEdepForSSB )
				{
					ssb2--;
				}
			}
			
			// Handle case of DSB. Whether SSB1 or SSB2 is "further ahead" in the chain doesn't matter
			if ( std::abs(nucl2-nucl1) <= fThresDistForDSB )
			{
				if ( ( edep2 >= fThresEdepForSSB ) &&
					( edep1 >= fThresEdepForSSB ) )
				{
					ssb1--;
					ssb2--;
					dsb++;
				}
			}
		}
	}
	
	//----------------------------------------------------------------------------------------------
	// Potentially some energy deposits (and thus SSBs) remain at end of the above iteration. Handle
	// these here, for either strand.
	//----------------------------------------------------------------------------------------------
	while ( !fVEdepStrand1[cluster].empty() )
	{
		nucl1 = fVEdepStrand1[cluster].begin()->first;
		edep1 = fVEdepStrand1[cluster].begin()->second;
		if ( edep1 >= fThresEdepForSSB )
		{
			ssb1++;
		}
		fVEdepStrand1[cluster].erase( fVEdepStrand1[cluster].begin() );
	}
	
	while ( !fVEdepStrand2[cluster].empty() )
	{
		nucl2 = fVEdepStrand2[cluster].begin()->first;
		edep2 = fVEdepStrand2[cluster].begin()->second;
		if ( edep2 >= fThresEdepForSSB )
		{
			ssb2++;
		}
		fVEdepStrand2[cluster].erase( fVEdepStrand2[cluster].begin() );
	}
	sb[0]=ssb1+ssb2;
	sb[1]=dsb;
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



