//**************************************************************************************************
// Author: Logan Montgomery
// Based on code written by:
//      - J Schuemann et al. (2019). DOI:10.1667/RR15226.1
//      - A McNamara et al. (2018). DOI:10.1088/1361-6560/aad8eb
//      - E Delage et al. (2015). DOI:10.1016/j.cpc.2015.02.026
//
// This class is a custom Topas ntuple scorer that records clustered DNA damage induced in a DNA
// model (VoxelizedNuclearDNA). 
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
    //----------------------------------------------------------------------------------------------
    // Constructor. Initialize member variables using a variety of methods. Specify which data is
    // output to the main data output file.
    //----------------------------------------------------------------------------------------------
    ScoreClusteredDNADamage(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
    //--------------------------------------------------------------------------------------------------
    // Destructor
    //--------------------------------------------------------------------------------------------------
    virtual ~ScoreClusteredDNADamage();

    //----------------------------------------------------------------------------------------------
    // Read in parameters from the Topas parameter file & save values in some member variables.
    //----------------------------------------------------------------------------------------------
    void ResolveParams();
    
    //--------------------------------------------------------------------------------------------------
    // Record energy depositions in the sensitive DNA volumes.
    //--------------------------------------------------------------------------------------------------
    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    
    //----------------------------------------------------------------------------------------------
    // Optionally process energy depositions to determine DNA damage yields (event-by-event)
    //----------------------------------------------------------------------------------------------
    void UserHookForEndOfEvent();

    //----------------------------------------------------------------------------------------------
    // Optionally process energy depositions to determine DNA damage yields (aggregated over run)
    //----------------------------------------------------------------------------------------------
    void UserHookForEndOfRun();

    
private:
    //----------------------------------------------------------------------------------------------
    // Erase contents of various output files (and their corresponding header files)
    //----------------------------------------------------------------------------------------------
    void ClearOutputFiles();

    //----------------------------------------------------------------------------------------------
    // Calculate cubic volume of component attached to the scorer
    //----------------------------------------------------------------------------------------------
    G4double CalculateComponentVolume();

    //----------------------------------------------------------------------------------------------
    // Convert dose threshold to energy threshold using the cubic volume and density of the
    // geometry component. 
    //----------------------------------------------------------------------------------------------
    G4double ConvertDoseThresholdToEnergy();

    //----------------------------------------------------------------------------------------------
    // This method outputs the details of the run to a header file and a data file.
    //----------------------------------------------------------------------------------------------
    void OutputRunSummaryToFile();

    //----------------------------------------------------------------------------------------------
    // This method outputs the details of scored Complex DSBs to a header file and data file.
    //----------------------------------------------------------------------------------------------
    void OutputComplexDSBToFile();

    //----------------------------------------------------------------------------------------------
    // This method outputs the details of scored Non-DSB clusters to a header file and data file.
    //----------------------------------------------------------------------------------------------
    void OutputNonDSBClusterToFile();

    //----------------------------------------------------------------------------------------------
    // This method transfers information from worker threads to the master thread, which allows
    // results be processed on a per-run basis.
    //----------------------------------------------------------------------------------------------
    void AbsorbResultsFromWorkerScorer(TsVScorer*);

    //----------------------------------------------------------------------------------------------
    // This method transfers the contents of a map of energy depositions from the worker thread to
    // the master thread
    //----------------------------------------------------------------------------------------------
    void AbsorbMapFromWorkerScorer(std::map<G4int,std::map<G4int,std::map<G4int, G4double>>>&,
        std::map<G4int,std::map<G4int,std::map<G4int, G4double>>>&);

    //----------------------------------------------------------------------------------------------
    // Process maps of energy depositions and record DNA damage yields to member variables.
    //----------------------------------------------------------------------------------------------
    void RecordDamage(); 

    //----------------------------------------------------------------------------------------------
    // This method resets member variable values
    //----------------------------------------------------------------------------------------------
    void ResetMemberVariables();

    //----------------------------------------------------------------------------------------------
    // This method resets variables that count the yields for various types of DNA damage.
    //----------------------------------------------------------------------------------------------
    void ResetDamageCounterVariables();
    
    //----------------------------------------------------------------------------------------------
    // Record bp indices of one type of simple DNA damage (SSB or BD) in a single strand to a vector
    //----------------------------------------------------------------------------------------------
    std::vector<G4int> RecordSimpleDamage(G4double,std::map<G4int, G4double>);

    //----------------------------------------------------------------------------------------------
    // Record indices of DSBs in a 1D vector
    //----------------------------------------------------------------------------------------------
    std::vector<G4int> RecordDSB();

    //----------------------------------------------------------------------------------------------
    // Process a single sequential vector of damage indices (labelled according to damage types) to
    // enable recording of two types of clustered DNA damage: Complex DSB and Non-DSB Clusters.
    //----------------------------------------------------------------------------------------------
    void RecordClusteredDamage();

    //----------------------------------------------------------------------------------------------
    // Add a new DNA damage site to a cluster
    //----------------------------------------------------------------------------------------------
    void AddDamageToCluster(DamageCluster&, G4int, G4int, G4int, G4bool);

    //----------------------------------------------------------------------------------------------
    // Add the details of a finalized DNA damage cluster to the appropriate member variables.
    //----------------------------------------------------------------------------------------------
    void RecordCluster(DamageCluster&);

    //----------------------------------------------------------------------------------------------
    // Combine class member vectors containing various types of damages into a single, ordered,
    // vector of all damges in a DNA fibre (both strands).
    //----------------------------------------------------------------------------------------------
    std::vector<std::array<G4int,3>> CombineSimpleDamage();

    //----------------------------------------------------------------------------------------------
    // Calculate the order of magnitude (base 10) of a positive integer value.
    //----------------------------------------------------------------------------------------------
    G4int CalculateIntegerMagnitude(G4int);


    // Methods used in validating scoring functionality
    void CreateFakeEnergyMap();
    void PrintDNADamageToConsole();
    

    //----------------------------------------------------------------------------------------------
    // Member variables
    //----------------------------------------------------------------------------------------------
    // Geometry stuff
    G4Material* fDNAMaterial;
    G4int fNumNucleosomePerFiber;
    G4int fNumBpPerNucleosome;
    G4int fNumBpMagnitude;
    G4int fParserStrand;
    G4int fParserResidue;
    G4bool fBuildNucleus;
    G4int fNumVoxelsPerSide;
    G4double fVoxelSideLength;

    // Thresholds for defining DNA damage
    G4double fThresEdepForSSB;
    G4double fThresEdepForBD;
    G4int fThresDistForDSB;
    G4int fThresDistForCluster;

    // Booleans
    G4bool fScoreClusters;
    G4bool fRecordDamagePerEvent;
    G4bool fRecordDamagePerFiber;
    G4bool fOutputHeaders;

    // Running counters
    G4int fNumEvents;
    G4double fTotalEdep;

    G4double fComponentVolume;

    // Dose threshold
    G4bool fUseDoseThreshold;
    G4double fDoseThreshold;
    G4double fEnergyThreshold;

    // # of threads used
    G4int fNumberOfThreads;

    // Output file parameters
    G4String fDelimiter;
    G4String fOutHeaderExtension;
    G4String fOutFileExtension;
    G4String fFileRunSummary;
    
    // These maps record energy deposited in bp in one of the strands of the DNA double helix
    std::map<G4int, G4double> fFiberMapEdepStrand1Backbone;
    std::map<G4int, G4double> fFiberMapEdepStrand2Backbone;
    std::map<G4int, G4double> fFiberMapEdepStrand1Base;
    std::map<G4int, G4double> fFiberMapEdepStrand2Base;
    
    // These maps energy deposited in bp in one of the strands of the DNA double helix in a
    // triple-nested map structure
    // map1 (key, map2) --> map2 (key, map3) --> map3 (key, double)
    // First index specifies the voxel
    // Second index specifies the fiber
    // Third index specifies the bp index
    std::map<G4int, std::map<G4int, std::map<G4int, G4double>>> fMapEdepStrand1Backbone;
    std::map<G4int, std::map<G4int, std::map<G4int, G4double>>> fMapEdepStrand2Backbone;
    std::map<G4int, std::map<G4int, std::map<G4int, G4double>>> fMapEdepStrand1Base;
    std::map<G4int, std::map<G4int, std::map<G4int, G4double>>> fMapEdepStrand2Base;
    
    // Various IDs
    G4int fThreadID;
    G4int fEventID;
    G4int fFiberID;
    G4int fVoxelID;

    // Damage yields
    G4int fTotalSSB;
    G4int fTotalBD;
    G4int fTotalDSB;
    G4int fTotalComplexDSB;
    G4int fTotalNonDSBCluster;
    G4int fTotalNonDSBClusterTandem;
    G4int fTotalNonDSBClusterBistranded;
    
    // Constant variables to identify damage types
    static const G4int fIdSSB = 0;
    static const G4int fIdBD = 1;
    static const G4int fIdDSB = 2;

    static const G4int fParentIndexFiber = 1; // Touchable history index for accessing DNA fiber
    static const G4int fParentIndexVoxelZ = 2; // Touchable history index for accessing Z voxels
    static const G4int fParentIndexVoxelX = 3; // Touchable history index for accessing X voxels
    static const G4int fParentIndexVoxelY = 4; // Touchable history index for accessing Y voxels

    // When reporting yields if, for example, aggregating over all events, set event number = -1
    static const G4int fAggregateValueIndicator = -1; 

    G4int fNumFibers;

    // Vectors to hold indices of simple damages
    std::vector<G4int> fIndicesSSB1;
    std::vector<G4int> fIndicesSSB2;
    std::vector<G4int> fIndicesBD1;
    std::vector<G4int> fIndicesBD2;

    // Vector to hold indices of DSBs
    std::vector<G4int>  fIndicesDSB;

    // Clustered damage handling
    std::vector<std::array<G4int,3>> fIndicesSimple; // first # is bp index, second # is 0 or 1 to represent SSB or BD
    
    G4String fFileComplexDSB;
    std::vector<G4int> fComplexDSBSizes; // Vector of lengths of complex DSB (in # of bp)
    std::vector<G4int> fComplexDSBNumSSB;
    std::vector<G4int> fComplexDSBNumBD;
    std::vector<G4int> fComplexDSBNumDSB;
    std::vector<G4int> fComplexDSBNumDamage;

    G4String fFileNonDSBCluster;
    std::vector<G4int> fNonDSBClusterSizes; // Vector of lengths of complex DSB (in # of bp)
    std::vector<G4int> fNonDSBClusterNumSSB;
    std::vector<G4int> fNonDSBClusterNumBD;
    std::vector<G4int> fNonDSBClusterNumDamage;
    std::vector<G4int> fNonDSBClusterBistranded;
};
#endif
