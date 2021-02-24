//**************************************************************************************************
// This class is essentially a wrapper for the GeoCalculation and GeoVolume classes. Those classes
// do the actual work for creating and placing DNA volumes. This class allows public methods from
// those classes to be accessed under a shared interface.
//**************************************************************************************************

#ifndef ScoreClusteredDNADamage_hh
#define ScoreClusteredDNADamage_hh

#include "TsVNtupleScorer.hh"

#include <map>


class G4Material;

class ScoreClusteredDNADamage : public TsVNtupleScorer
{
public:
    ScoreClusteredDNADamage(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
    virtual ~ScoreClusteredDNADamage();
    
    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    
    void UserHookForEndOfEvent();
    
private:
    void ComputeStrandBreaks(G4int*, G4int);

    G4int CalculateIntegerMagnitude(G4int);
    
private:
    // G4Material* fStrand1Material;
    // G4Material* fStrand2Material;
    G4Material* fDNAMaterial;

    G4int fNucleoNum;
    G4int fBpNum;
    G4int fNumBpMagnitude;
    G4int fParserStrand;
    G4int fParserResidue;

    G4int fThresDistForDSB;
    G4double fThresEdepForSSB;
    G4double fThresEdepForBD;

    G4int fNumEdeps1;
    G4double fTotalEdep1;
    G4int fNumEdeps2;
    G4double fTotalEdep2;
    G4int fNumEdepsBD1;
    G4double fTotalEdepBD1;
    G4int fNumEdepsBD2;
    G4double fTotalEdepBD2;

    // G4int fThreshDistCluster;
    // G4double fThreshEdepForBD;
    

    // Records energy deposited in bp in one of the strands of the DNA double helix in a
    // double-nested map structure
    // map1 (key, map2) --> map2 (key, double)
    // First index specifies someting to do with variance reduction / track splitting
    // Second index specifies the bp index
    std::map<G4int, std::map<G4int, G4double> > fVEdepStrand1;
    std::map<G4int, std::map<G4int, G4double> > fVEdepStrand2;
    
    // Records energy deposited in bp in one of the strands of the DNA double helix in a
    // triple-nested map structure
    // map1 (key, map2) --> map2 (key, map3) --> map3 (key, double)
    // First index specifies the DNA fibre
    // Second index specifies someting to do with variance reduction / track splitting
    // Third index specifies the bp index
    std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand1;
    std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand2;
    
    G4int fNbOfAlgo;
    G4int fEventID;
    G4int fDNAParent; // if any
    G4int fSSB;
    G4int fDSB;
    G4int fBD;

    G4int fBasePairDepth;
    
    G4String fStrand1VolumeName;
    G4String fStrand2VolumeName;
    
    
};
#endif
