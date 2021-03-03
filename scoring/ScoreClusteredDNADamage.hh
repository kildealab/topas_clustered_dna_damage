//**************************************************************************************************
// This class is essentially a wrapper for the GeoCalculation and GeoVolume classes. Those classes
// do the actual work for creating and placing DNA volumes. This class allows public methods from
// those classes to be accessed under a shared interface.
//**************************************************************************************************

#ifndef ScoreClusteredDNADamage_hh
#define ScoreClusteredDNADamage_hh

#include "TsVNtupleScorer.hh"

#include <map>

struct DamageCluster;

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

    std::vector<G4int> RecordSimpleDamage(G4double,std::map<G4int,std::map<G4int, G4double>>);
    // std::vector<std::array<G4int,2>> RecordDSB(G4int,std::vector<G4int>,std::vector<G4int>);
    // std::vector<std::array<G4int,2>> RecordDSB();
    std::vector<G4int> RecordDSB1D();
    // std::vector<std::array<G4int,2>> RecordComplexDSB();
    std::vector<std::vector<G4int>> RecordComplexDSB();
    void RecordClusteredDamage();
    void AddDamageToCluster(DamageCluster&, G4int, G4int, G4bool);
    void RecordCluster(DamageCluster&);

    std::vector<std::array<G4int,2>> CombineSimpleDamage();

    // Methods used in validating scoring functionality
    void CreateFakeEnergyMap();
    void PrintDNADamageToConsole();
    
private:
    // G4Material* fStrand1Material;
    // G4Material* fStrand2Material;
    G4Material* fDNAMaterial;

    G4int fNucleoNum;
    G4int fBpNum;
    G4int fNumBpMagnitude;
    G4int fParserStrand;
    G4int fParserResidue;

    G4double fThresEdepForSSB;
    G4double fThresEdepForBD;
    G4int fThresDistForDSB;
    G4int fThresDistForCluster;

    // G4int fNumEdeps1;
    // G4double fTotalEdep1;
    // G4int fNumEdeps2;
    // G4double fTotalEdep2;
    // G4int fNumEdepsBD1;
    // G4double fTotalEdepBD1;
    // G4int fNumEdepsBD2;
    // G4double fTotalEdepBD2;
    
    // Records energy deposited in bp in one of the strands of the DNA double helix in a
    // double-nested map structure
    // map1 (key, map2) --> map2 (key, double)
    // First index specifies someting to do with variance reduction / track splitting
    // Second index specifies the bp index
    std::map<G4int, std::map<G4int, G4double> > fVEdepStrand1Backbone;
    std::map<G4int, std::map<G4int, G4double> > fVEdepStrand2Backbone;
    std::map<G4int, std::map<G4int, G4double> > fVEdepStrand1Base;
    std::map<G4int, std::map<G4int, G4double> > fVEdepStrand2Base;
    
    // Records energy deposited in bp in one of the strands of the DNA double helix in a
    // triple-nested map structure
    // map1 (key, map2) --> map2 (key, map3) --> map3 (key, double)
    // First index specifies the DNA fibre
    // Second index specifies someting to do with variance reduction / track splitting
    // Third index specifies the bp index
    std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand1Backbone;
    std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand2Backbone;
    std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand1Base;
    std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand2Base;
    
    G4int fNbOfAlgo;
    G4int fEventID;
    G4int fDNAParent; // if any

    G4int fTotalSSB;
    G4int fTotalBD;
    G4int fTotalDSB;
    G4int fTotalComplexDSB;
    G4int fTotalNonDSBCluster;

    G4int fBasePairDepth;
    
    // Constant variables to identify damage types
    static const G4int fIdSSB = 0;
    static const G4int fIdBD = 1;
    static const G4int fIdDSB = 2;

    // Vectors to hold indices of simple damages
    std::vector<G4int> fIndicesSSB1;
    std::vector<G4int> fIndicesSSB2;
    std::vector<G4int> fIndicesBD1;
    std::vector<G4int> fIndicesBD2;

    // Vector to hold indices of DSBs
    std::vector<std::array<G4int,2>> fIndicesDSB;
    std::vector<G4int>  fIndicesDSB1D;

    // Clustered damage handling
    std::vector<std::array<G4int,2>> fIndicesSimple; // first # is bp index, second # is 0 or 1 to represent SSB or BD
    std::vector<G4int> fComplexDSBSizes; // Vector of lengths of complex DSB (in # of bp)
    std::vector<G4int> fComplexDSBNumSSB;
    std::vector<G4int> fComplexDSBNumBD;
    std::vector<G4int> fComplexDSBNumDSB;
    std::vector<G4int> fComplexDSBNumDamage;

    std::vector<G4int> fNonDSBClusterSizes; // Vector of lengths of complex DSB (in # of bp)
    std::vector<G4int> fNonDSBClusterNumSSB;
    std::vector<G4int> fNonDSBClusterNumBD;
    std::vector<G4int> fNonDSBClusterNumDamage;
};
#endif
