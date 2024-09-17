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

#ifndef ScoreSDD_hh
#define ScoreSDD_hh

#include "TsVNtupleScorer.hh"

#include <map>

class G4Material;

class ScoreSDD : public TsVNtupleScorer
{
public:
    //----------------------------------------------------------------------------------------------
    // Constructor. Initialize member variables using a variety of methods. Specify which data is
    // output to the main data output file.
    //----------------------------------------------------------------------------------------------
    ScoreSDD(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

    //--------------------------------------------------------------------------------------------------
    // Destructor
    //--------------------------------------------------------------------------------------------------
    virtual ~ScoreSDD();

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


     //----------------------------------------------------------------------------------------------
    // This method outputs all events to file. (NICO)
    //----------------------------------------------------------------------------------------------

    void OutputEvent(G4int ModeOfAction, G4int PositionOnBp, G4ThreeVector position, G4int voxelID, G4int index);

    void OutputSDDHeader();


private:
    //--------------------------------------------------------------------------------------------------
    // This helper method checks whether an element is in a vector.
    //--------------------------------------------------------------------------------------------------
    G4bool IsElementInVector(G4int, std::vector<G4int>&);

    void RemoveElementFromVector(G4int, std::vector<G4int>&);

    G4bool IsDamageInflicted(G4int, G4int);


    void PrintStepInfo(G4Step*);

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

    void AbsorbDamageBlock(
    std::vector<std::pair<int, std::string>> &masterMap,
    std::vector<std::pair<int, std::string>> &workerMap);

    //--------------------------------------------------------------------------------------------------
    // This method merges and resolves duplicates of the damage yields from direct and indirect damage.
    //--------------------------------------------------------------------------------------------------
    std::vector<G4int> MergeDamageIndices(std::vector<G4int>&,std::vector<G4int>&);

    //----------------------------------------------------------------------------------------------
    // This method resets member variable values
    //----------------------------------------------------------------------------------------------
    void ResetMemberVariables();

    //----------------------------------------------------------------------------------------------
    // Record bp indices of one type of simple DNA damage (SSB or BD) in a single strand to a vector
    //----------------------------------------------------------------------------------------------
    std::vector<G4int> RecordSimpleDamage(G4double ThreshEDep, std::map<G4int, std::pair<G4double, G4ThreeVector>> mapEDep , G4int VoxelID, G4int PositionOnBp);
    std::vector<G4int> RecordBD(std::map<G4int, G4ThreeVector> IndMap , G4int VoxelID, G4int PositionOnBp);

    //----------------------------------------------------------------------------------------------
    // Record indices of DSBs in a 1D vector
    //----------------------------------------------------------------------------------------------
    std::vector<G4int> RecordDSB(G4int);
    std::vector<G4int> RecordHybridDSB();

    //----------------------------------------------------------------------------------------------
    // Process a single sequential vector of damage indices (labelled according to damage types) to
    // enable recording of two types of clustered DNA damage: Complex DSB and Non-DSB Clusters.
    //----------------------------------------------------------------------------------------------
    void RecordClusteredDamage();


    //----------------------------------------------------------------------------------------------
    // Combine class member vectors containing various types of damages into a single, ordered,
    // vector of all damges in a DNA fibre (both strands).
    //----------------------------------------------------------------------------------------------
    std::vector<std::array<G4int,3>> CombineSimpleDamage();

    //----------------------------------------------------------------------------------------------
    // Calculate the order of magnitude (base 10) of a positive integer value.
    //----------------------------------------------------------------------------------------------
    G4int CalculateIntegerMagnitude(G4int);

    void Print1DVectorContents(std::vector<G4int>);
    

    //----------------------------------------------------------------------------------------------
    // Member variables
    //----------------------------------------------------------------------------------------------

    G4int loc_in_bp;

    // Geometry stuff
    G4Material* fDNAMaterial;
    G4Material* fHistoneMaterial;
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
    G4bool fRecordDamagePerFiber;
    G4bool fOutputHeaders;
    G4bool fIncludeDirectDamage;
    G4bool fIncludeIndirectDamage;
    G4bool fHasChemistryModule;
    G4bool DsbOnly;

    // Running counters
    G4int fNumEvents;
    G4double fTotalEdep;
    G4int fNumProcessHitsCalls;

    G4double fComponentVolume;

    // Dose threshold
    G4bool fUseDoseThreshold;
    G4double fDoseThreshold;
    G4double fEnergyThreshold;

    //---Added by Nicolas---
    //G4double doseDep;
    //---------------

    // # of threads used
    G4int fNumberOfThreads;


    // Output file parameters
    G4String fDelimiter;
    G4String fOutHeaderExtension;
    G4String fOutFileExtension;
    G4String fFileRunSummary;
    //--Added by Nicolas---
    G4String fFileAllEvents;
    G4String fSDDOutputFile;
    G4String fRealDoseFile;
    G4String fMapFile;
    G4String MeanEnergy;
    G4String IncidentParticles;
    G4String PrimaryParticle;
    //------------------------

    // These maps record energy deposited in bp in one of the strands of the DNA double helix
    std::map<G4int, std::pair<G4double, G4ThreeVector>> fFiberMapEdepStrand1Backbone;
    std::map<G4int, std::pair<G4double, G4ThreeVector>> fFiberMapEdepStrand2Backbone;
    std::map<G4int, std::pair<G4double, G4ThreeVector>> fFiberMapEdepStrand1Base;
    std::map<G4int, std::pair<G4double, G4ThreeVector>> fFiberMapEdepStrand2Base;

    std::map<G4int, G4ThreeVector> fFiberMapPosStrand1Backbone;
    std::map<G4int, G4ThreeVector> fFiberMapPosStrand2Backbone;
    std::map<G4int, G4ThreeVector> fFiberMapPosStrand1Base;
    std::map<G4int, G4ThreeVector> fFiberMapPosStrand2Base;

    // These maps energy deposited in bp in one of the strands of the DNA double helix in a
    // triple-nested map structure
    // map1 (key, map2) --> map2 (key, map3) --> map3 (key, double)
    // First index specifies the voxel
    // Second index specifies the fiber
    // Third index specifies the bp index
    std::map<G4int, std::map<G4int, std::map<G4int, std::pair<G4double, G4ThreeVector>>>> fMapEdepStrand1Backbone;
    std::map<G4int, std::map<G4int, std::map<G4int, std::pair<G4double, G4ThreeVector>>>> fMapEdepStrand2Backbone;
    std::map<G4int, std::map<G4int, std::map<G4int, std::pair<G4double, G4ThreeVector>>>> fMapEdepStrand1Base;
    std::map<G4int, std::map<G4int, std::map<G4int, std::pair<G4double, G4ThreeVector>>>> fMapEdepStrand2Base;

    std::map<G4int, std::map<G4int, std::map<G4int, G4ThreeVector>>> fMapIndDamageStrand1Backbone;
    std::map<G4int, std::map<G4int, std::map<G4int, G4ThreeVector>>> fMapIndDamageStrand2Backbone;
    std::map<G4int, std::map<G4int, std::map<G4int, G4ThreeVector>>> fMapIndDamageStrand1Base;
    std::map<G4int, std::map<G4int, std::map<G4int, G4ThreeVector>>> fMapIndDamageStrand2Base;

    
    // Various IDs
    G4int fThreadID;
    G4int fEventID;
    G4int fFiberID;
    G4int fVoxelID;

    // Molecule IDs
    G4int fMoleculeID_OH;
    G4int fMoleculeID_OHm;
    G4int fMoleculeID_e_aq;
    G4int fMoleculeID_H2;
    G4int fMoleculeID_H3Op;
    G4int fMoleculeID_H;
    G4int fMoleculeID_H2O2;

    G4int fMoleculeID_HO2;
    G4int fMoleculeID_HO2m;
    G4int fMoleculeID_O2;
    G4int fMoleculeID_O2m;

    std::vector<G4int> fSpeciesToKillByDNAVolumes;
    std::vector<G4int> fspeciesToKillByHistones;
    G4bool fHistonesAsScavenger;


    // Constant variables to identify damage types
    static const G4int fIdSSB = 0;
    static const G4int fIdBD = 1;
    static const G4int fIdDSB = 2;

    static const G4int fIdDirect = 0;
    static const G4int fIdIndirect = 1;
    static const G4int fIdHybrid = 2;

    // Constant variables to identify residual DNA volume types after parsing volID in ProcessHits
    static const G4int fVolIdPhosphate = 0;
    static const G4int fVolIdDeoxyribose = 1;
    static const G4int fVolIdBase = 2;

    static const G4int fParentIndexFiber = 1; // Touchable history index for accessing DNA fiber
    static const G4int fParentIndexVoxelZ = 2; // Touchable history index for accessing Z voxels
    static const G4int fParentIndexVoxelX = 3; // Touchable history index for accessing X voxels
    static const G4int fParentIndexVoxelY = 4; // Touchable history index for accessing Y voxels

    // When reporting yields if, for example, aggregating over all events, set event number = -1
    static const G4int fAggregateValueIndicator = -1;

    G4int fNumFibers;

    // Map of moleculeID to corresponding probability of inflicting damage when molecule reacts with DNA volume
    std::map<G4int, G4float> fMoleculeDamageProb_SSB; // backbone damage
    std::map<G4int, G4float> fMoleculeDamageProb_BD; // base damage

    // Vectors to hold indices of simple damages
    std::vector<G4int>* fIndices;
    std::vector<G4int>* fIndicesSSB1;
    std::vector<G4int>* fIndicesSSB2;
    // std::vector<G4int>* fIndicesBD1;
    // std::vector<G4int>* fIndicesBD2;

    std::vector<G4int> fIndicesSSB1_merged;
    std::vector<G4int> fIndicesSSB2_merged;
    // std::vector<G4int> fIndicesBD1_merged;
    // std::vector<G4int> fIndicesBD2_merged;

    std::vector<G4int> fIndicesSSB1_direct;
    std::vector<G4int> fIndicesSSB2_direct;
    std::vector<G4int> fIndicesBD1_direct;
    std::vector<G4int> fIndicesBD2_direct;

    std::vector<G4int> fIndicesSSB1_indirect;
    std::vector<G4int> fIndicesSSB2_indirect;
    std::vector<G4int> fIndicesBD1_indirect;
    std::vector<G4int> fIndicesBD2_indirect;

    // Vector to hold indices of DSBs
    std::vector<G4int>  fIndicesDSB;
    std::vector<G4int>  fIndicesDSB_hybrid;
    std::vector<G4int>  fIndicesDSB_direct;
    std::vector<G4int>  fIndicesDSB_indirect;

    // Clustered damage handling
    std::vector<std::array<G4int,3>> fIndicesSimple; // first # is bp index, second # is 0 or 1 to represent SSB or BD


    void RecordDamageSDD();
    void ComputeDSB();
    std::vector<G4int> FindBaseDamage(std::vector<G4int> indices, G4int MinBp, G4int MaxBp);
    std::vector<std::array<G4int,3>> AllSimpleDamageSDD();
    void  SimpleClusterSDD();

    std::vector<std::array<G4int,3>> dsbBlock;

    static const G4int fIdBckBone1 = 1;
    static const G4int fIdBckBone2 = 4;


    static const G4int fIdBase1 = 2;
    static const G4int fIdBase2 = 3;

    G4int min_site;
    G4int max_site;

     std::vector<std::array<G4int,3>> fAllSimpleDamage;

    std::vector<G4int> InBlock_direct_bd1;
    std::vector<G4int> InBlock_indirect_bd1;
    std::vector<G4int> InBlock_direct_bd2;
    std::vector<G4int> InBlock_indirect_bd2;

    G4bool NewEvent;
    int field1;

    void RecordBlock(G4int EventID, std::vector<std::array<G4int,3>> block, G4int DSBcount);
    G4ThreeVector GetXYZ(G4int BpPosition, G4int PositionOnBp, G4int Mode);


    std::vector<std::pair<int, std::string>> block_pairs;

    void PrintSDDdata();

    std::map<G4int, std::vector<G4int> > ChromoMapping;
    int GetChromoMap();


};
#endif

