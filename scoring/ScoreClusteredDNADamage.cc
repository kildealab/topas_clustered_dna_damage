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
#include <algorithm>

#include <map>

//--------------------------------------------------------------------------------------------------
// Struct used to hold parameters of interest for a given cluter of DNA damage.
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
	fThresDistForCluster = 40;

	fTotalSSB = 0;
	fTotalBD = 0;
	fTotalDSB = 0;
	fTotalComplexDSB = 0;
	fTotalNonDSBCluster = 0;

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
	
	// This is for variance reduction
	if ( fPm->ParameterExists(GetFullParmName("NumberOfSplit")) )
		fNbOfAlgo = fPm->GetIntegerParameter(GetFullParmName("NumberOfSplit"));
	
	//fThresEdepForSSB /= eV;
	
	fNtuple->RegisterColumnI(&fEventID, "Event number"); // ID of primary particle / event / history
	// fNtuple->RegisterColumnI(&fDNAParent, "DNA parent geometry"); // Unsure. Not initiated in header either
	fNtuple->RegisterColumnI(&fTotalSSB, "Single strand breaks"); // Number of SSB caused by this primary particle
	fNtuple->RegisterColumnI(&fTotalDSB, "Double strand breaks"); // Number of simple DSB caused by this primary particle
	fNtuple->RegisterColumnI(&fTotalBD, "Base damages"); // Number of BD caused by this primary particle
	fNtuple->RegisterColumnI(&fTotalComplexDSB, "Complex DSBs"); // Number of Complex DSB caused by this primary particle
	fNtuple->RegisterColumnI(&fTotalNonDSBCluster, "Non-DSB clusters"); // Number of Non-DSB Clusters caused by this primary particle
	
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
// simply record energy deposited in a given bp index, in either the backbone or base of either 
// strand # 1 or strand #2 (4 separate energy depositio maps). Multiple energy depositions occurring
// in the same volume are added together and recorded as one.
//
// Note that the copy number of the volume is used to determine in which volume the energy
// deposition took place. This is faster than using string comparisons. Note this method is only
// called for energy depositions in the sensitive volumes (i.e. residues) by using material
// filtering, as defined in the parameter file with "OnlyIncludeIfInMaterial" parameter.
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
			// NEED TO UPDATE THIS TO MATCH NEW CODE FOR WHEN NOT DOING VARIANCE REDUCTION
			if (index > 2) {
				if ( num_strand == 0 ){
					fGenVEdepStrand1Backbone[parentDepth][index-3][num_nucleotide] += edep;
				}
				else{
					fGenVEdepStrand2Backbone[parentDepth][index-3][num_nucleotide] += edep;
				}
			} 
			// Not doing variance reduction
			else {
				if ( num_strand == 0 ){
					for ( int i = 0; i < fNbOfAlgo; i++ ){
						if (num_res == 0 || num_res == 1){
							fGenVEdepStrand1Backbone[parentDepth][i][num_nucleotide] += edep;
						}
						else if (num_res == 2) {
							fGenVEdepStrand1Base[parentDepth][i][num_nucleotide] += edep;
						}
					}
				}
				else{
					for ( int i = 0; i < fNbOfAlgo; i++ ){
						if (num_res == 0 || num_res == 1){
							fGenVEdepStrand2Backbone[parentDepth][i][num_nucleotide] += edep;
						}
						else if (num_res == 2) {
							fGenVEdepStrand2Base[parentDepth][i][num_nucleotide] += edep;
						}
					}
				}
			}
			
			return true;
		}

	}
	return false;
}


//--------------------------------------------------------------------------------------------------
// This method is called at the end of the event (primary particle & its secondaries). Process the
// energy deposoition maps populated by the ProcessHits() method to determine the number of various
// types of DNA damage. Record these results to output.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::UserHookForEndOfEvent() {
	fEventID = GetEventID();

	// Include following line if want to create a fake, predefined energy map to validate scoring
	// CreateFakeEnergyMap();

	// Include if want to print out map of energy depositions to validate algorithm is functioning
	// properly.
	// std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand1Backbone_COPY = fGenVEdepStrand1Backbone;
	// std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand2Backbone_COPY = fGenVEdepStrand2Backbone;
	// std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand1Base_COPY = fGenVEdepStrand1Base;
	// std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand2Base_COPY = fGenVEdepStrand2Base;

	// Iterate over all DNA fibres 
	for ( auto& energyAtStrands : fGenVEdepStrand1Backbone ) {

		G4int parentDepth = energyAtStrands.first;

		// Unsure why this necessary. Perhaps is faster to do this rather than repeatedly indexing
		// the extra layer of mappng associated with fGenVEdepStrand1Backbone/2
		fVEdepStrand1Backbone = fGenVEdepStrand1Backbone[parentDepth];
		fVEdepStrand2Backbone = fGenVEdepStrand2Backbone[parentDepth];
		fVEdepStrand1Base = fGenVEdepStrand1Base[parentDepth];
		fVEdepStrand2Base = fGenVEdepStrand2Base[parentDepth];
		fDNAParent = parentDepth;

		// Iterate over track splitting (only 1 iteration if no splitting), compute the # of 
		// various types of DNA damage, and fill the ntuple scorer accordingly
		for ( int i = 0; i < fNbOfAlgo; i++ ) {
			fIndicesSSB1 = RecordSimpleDamage(fThresEdepForSSB,fVEdepStrand1Backbone);
			fIndicesSSB2 = RecordSimpleDamage(fThresEdepForSSB,fVEdepStrand2Backbone);
			fIndicesBD1 = RecordSimpleDamage(fThresEdepForBD,fVEdepStrand1Base);
			fIndicesBD2 = RecordSimpleDamage(fThresEdepForBD,fVEdepStrand2Base);

			// fIndicesDSB = RecordDSB();
			fIndicesDSB1D = RecordDSB1D();

			fTotalSSB += fIndicesSSB1.size() + fIndicesSSB2.size();
			fTotalBD += fIndicesBD1.size() + fIndicesBD2.size();
			// fTotalDSB += fIndicesDSB.size();
			// fTotalDSB += fIndicesDSB1D.size()/2;
			fTotalDSB += fIndicesDSB1D.size(); // Note is currently 2x number of DSB. Division by 2 is performed after clustering

			fIndicesSimple = CombineSimpleDamage();
			RecordClusteredDamage();

			if (fTotalSSB > 0 || fTotalBD > 0 || fTotalDSB > 0 || fTotalComplexDSB > 0 || fTotalNonDSBCluster > 0) {
				// If want to output results to command line
				// PrintDNADamageToConsole();
				fNtuple->Fill();
			}
		}

		// Include if want to print out map of energy depositions to validate algorithm is
		// functioning properly.
		// if (fSSB > 0 && fDSB > 0) {
		// if (fEventID == 57) {
		// 	G4cout << "#################################################################################################################################" << G4endl;
		// 	G4cout << "Event #" << fEventID << G4endl;
		// 	int colwidth = 15;
		// 	for (G4int i = 0; i < 18000; i++) {
		// 		G4cout << std::left << "i = " << std::setw(colwidth) << i << " | " << std::setw(colwidth) << fGenVEdepStrand1Backbone_COPY[0][0][i]/eV << std::setw(colwidth) << fGenVEdepStrand1Base_COPY[0][0][i]/eV << std::setw(colwidth) << fGenVEdepStrand2Base_COPY[0][0][i]/eV << std::setw(colwidth) << fGenVEdepStrand2Backbone_COPY[0][0][i]/eV << G4endl;
		// 		// if (fGenVEdepStrand1Backbone_COPY[0][0][i]/eV > 0 || fGenVEdepStrand1Base_COPY[0][0][i]/eV > 0 || fGenVEdepStrand2Base_COPY[0][0][i]/eV > 0 || fGenVEdepStrand2Backbone_COPY[0][0][i]/eV > 0) {
		// 		// 	G4cout << std::left << "i = " << i << " | " << std::setw(colwidth) << fGenVEdepStrand1Backbone_COPY[0][0][i]/eV << std::setw(colwidth) << fGenVEdepStrand1Base_COPY[0][0][i]/eV << std::setw(colwidth) << fGenVEdepStrand2Base_COPY[0][0][i]/eV << std::setw(colwidth) << fGenVEdepStrand2Backbone_COPY[0][0][i]/eV << G4endl;
		// 		// }
		// 	}
		// 	G4cout << "#################################################################################################################################" << G4endl;
		// }

		// Clear member variables for next fibre
		fVEdepStrand1Backbone.erase(fVEdepStrand1Backbone.begin(), fVEdepStrand1Backbone.end());
		fVEdepStrand2Backbone.erase(fVEdepStrand2Backbone.begin(), fVEdepStrand2Backbone.end());
		fVEdepStrand1Base.erase(fVEdepStrand1Base.begin(), fVEdepStrand1Base.end());
		fVEdepStrand2Base.erase(fVEdepStrand2Base.begin(), fVEdepStrand2Base.end());

		fTotalSSB = 0;
		fTotalBD = 0;
		fTotalDSB = 0;
		fTotalComplexDSB = 0;
		fTotalNonDSBCluster = 0;

		fComplexDSBSizes.clear();
		fComplexDSBNumSSB.clear();
		fComplexDSBNumBD.clear();
		fComplexDSBNumDSB.clear();
		fComplexDSBNumDamage.clear();

		fNonDSBClusterSizes.clear();
		fNonDSBClusterNumSSB.clear();
		fNonDSBClusterNumBD.clear();
		fNonDSBClusterNumDamage.clear();
	}
	// Clear member variables for next event
	fGenVEdepStrand1Backbone.erase(fGenVEdepStrand1Backbone.begin(), fGenVEdepStrand1Backbone.end());
	fGenVEdepStrand2Backbone.erase(fGenVEdepStrand2Backbone.begin(), fGenVEdepStrand2Backbone.end());
	fGenVEdepStrand1Base.erase(fGenVEdepStrand1Base.begin(), fGenVEdepStrand1Base.end());
	fGenVEdepStrand2Base.erase(fGenVEdepStrand2Base.begin(), fGenVEdepStrand2Base.end());
}


//--------------------------------------------------------------------------------------------------
// Record bp indices of one type of simple DNA damage (SSB or BD) in a single strand to a 1D vector.
// This method deletes the content in the provided map.
//--------------------------------------------------------------------------------------------------
std::vector<G4int> ScoreClusteredDNADamage::RecordSimpleDamage(G4double ThreshEDep,
	std::map<G4int,std::map<G4int,G4double>> mapEDep)
{
	std::vector<G4int> indicesDamage;

	G4int indexBP;
	G4double eDep;

	// Loop through map of energy depositions. Add indices ofenergy depositions over the provided
	// threshold to a vector, which is returned. The entries in the map are deleted as they are
	// processed
	while ( !mapEDep[0].empty() )
	{
		indexBP = mapEDep[0].begin()->first;
		eDep = mapEDep[0].begin()->second;

		if (eDep >= ThreshEDep) {
			indicesDamage.push_back(indexBP);
		}
		mapEDep[0].erase(mapEDep[0].begin());
	}

	return indicesDamage;
}


//--------------------------------------------------------------------------------------------------
// Record indices of DSBs in a 1D vector. Size should always be even, corresponding to two damage
// sites per DSB. The lowest bp index is record first, regardless of whether in strand 1 or 2.
//--------------------------------------------------------------------------------------------------
std::vector<G4int> ScoreClusteredDNADamage::RecordDSB1D()
{
	std::vector<G4int> indicesDSB1D;

	// Fake data for testing
	// fIndicesSSB1 = {0,3,5,9};
	// fIndicesSSB2 = {2,8};
	// fThresDistForDSB = 1;

	std::vector<G4int>::iterator site1 = fIndicesSSB1.begin();
	std::vector<G4int>::iterator site2 = fIndicesSSB2.begin();

	while (site1 != fIndicesSSB1.end() && site2 != fIndicesSSB2.end()) {
		G4int siteDiff = *site2 - *site1;

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
// Process a single member vector of damage indices (labelled according to damage types) to enable
// recording of two types of clustered DNA damage: Complex DSB and Non-DSB Clusters. Definitions are
// equivalent, except the former contains one or more DSB. Clustering is performed by calculating
// distances (in bp) between subsequent damage sites and comparing with a maximum clustering
// distance.
//--------------------------------------------------------------------------------------------------
void ScoreClusteredDNADamage::RecordClusteredDamage()
{
	// Fake data for testing
	// fIndicesSimple = {
	// 	{1,1},
	// 	{2,0},
	// 	{5,1},
	// 	{5,2},
	// 	{6,2},
	// 	{6,1},
	// 	{7,0},
	// 	{11,1},
	// 	{12,0},
	// 	{13,1},
	// 	{15,2},
	// 	{17,2},
	// 	{18,0},
	// 	{21,2},
	// 	{22,2},
	// 	{23,1},
	// 	{26,0},
	// 	{29,1},
	// 	{31,1},
	// 	{35,0}
	// };
	// fThresDistForCluster = 2;

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
		std::array<G4int,2> sitePrev = *(site-1);
		std::array<G4int,2> siteCur = *site;

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
			// Previous site was part of a cluster, but now cluster has ended
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

	// Finally, divide number of DSB by two
	fTotalDSB = fTotalDSB/2;
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
// Add a new damage site to a particular cluster
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
	// Throw and error if an unrecognized damage type is added
	else {
		G4cout << "An error has arisen while scoring clustered DNA damage." << G4endl;
		G4cout << "The following integer damage label is unrecognized:" << G4endl;
		G4cout << damageType << G4endl;
		exit(0);
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
	// fIndicesDSB1D = {3,5,11,14};

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
	std::vector<G4int>::iterator itDSB = fIndicesDSB1D.begin();
	itSimple = indicesSimple.begin();
	while (itDSB != fIndicesDSB1D.end() && itSimple != indicesSimple.end()) {
		if (*itDSB < (*itSimple)[0]) {
			itSimple = indicesSimple.insert(itSimple,{*itDSB,fIdDSB});
			itDSB++;
		}
		else {
			itSimple++;
		}
	}
	while (itDSB != fIndicesDSB1D.end()) {
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


	fGenVEdepStrand1Backbone[0][0] = back1;
	fGenVEdepStrand1Base[0][0] = base1;
	fGenVEdepStrand2Base[0][0] = base2;
	fGenVEdepStrand2Backbone[0][0] = back2;
	// fVEdepStrand1Backbone[0] = back1;
	// fVEdepStrand1Base[0] = base1;
	// fVEdepStrand2Base[0] = base2;
	// fVEdepStrand2Backbone[0] = back2;
}


// //--------------------------------------------------------------------------------------------------
// // Record indices of DSBs in a vector, whererin each vector element is a 2D array containing the
// // both indices of a given DSB (i.e. on strand 1 and 2).
// //--------------------------------------------------------------------------------------------------
// // std::vector<std::array<G4int,2>> ScoreClusteredDNADamage::RecordDSB(G4int,std::vector<G4int>,
// // 	std::vector<G4int>)
// std::vector<std::array<G4int,2>> ScoreClusteredDNADamage::RecordDSB()
// {
// 	std::vector<std::array<G4int,2>> indicesDSB;

// 	// fIndicesSSB1 = {0,3,5,9};
// 	// fIndicesSSB2 = {2,8};
// 	// fThresDistForDSB = 1;

// 	std::vector<G4int>::iterator site1 = fIndicesSSB1.begin();
// 	std::vector<G4int>::iterator site2 = fIndicesSSB2.begin();

// 	while (site1 != fIndicesSSB1.end() && site2 != fIndicesSSB2.end()) {
// 		G4int siteDiff = *site2 - *site1;

// 		// Damage in site 2 is within range of site 1 to count as DSB (either before or after)
// 		if (abs(siteDiff) <= fThresDistForDSB){
// 			// Damage in site 1 is earlier or parallel to damage in site 2
// 			if (*site1 <= *site2) {
// 				indicesDSB.push_back({*site1,*site2});
// 			}
// 			// Damage in site 2 is earlier to damage in site 1
// 			else {
// 				indicesDSB.push_back({*site2,*site1});
// 			}
// 			site1 = fIndicesSSB1.erase(site1);
// 			site2 = fIndicesSSB2.erase(site2);
// 		}
// 		// Damage in site 2 is earlier than site 1 and outside range to be considered DSB
// 		else if (siteDiff < 0) {
// 			site2++;
// 		}
// 		// Damage in site 2 is later than site 1 and outside range to be considered DSB.
// 		else { // if siteDiff > 0
// 			site1++;
// 		}
// 	}

// 	return indicesDSB;
// }
