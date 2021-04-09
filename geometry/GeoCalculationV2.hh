//**************************************************************************************************
// Author: Logan Montgomery
// Based on code written by:
//      - J Schuemann et al. (2019). DOI:10.1667/RR15226.1
//      - A McNamara et al. (2018). DOI:10.1088/1361-6560/aad8eb
//      - S Meylan et al. (2017). DOI:10.1038/s41598-017-11851-4
//
// This class performs the calculations necessary for placing nucleotide base pairs and nucloesomes.
// The size of the individual volumes comprising nucleotides and histones are defined using the
// Initialize() function. CalculateNucleosomePosition() and CalculateDNAPosition() are used to
// place three "basis" nucleosomes and their accompanying 200 bp each (incl linker DNA).
//**************************************************************************************************

#ifndef GEOCALCULATIONV2_HH
#define GEOCALCULATIONV2_HH

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <map>

//--------------------------------------------------------------------------------------------------
// Structure containting coordinates of volumes in a DNA base pair.
//--------------------------------------------------------------------------------------------------
struct DNAPlacementData;

//--------------------------------------------------------------------------------------------------
// vector contains vector which contains DNA volume positions for one base pair
// Schema:
//--------------------------------------------------------------------------------------------------
//     Vector                SubVector containing DNAPlacementData
//     |
//     |nucleosome 1 --------| Position informations about base pair 1 (4 volumes)
//     |                     | Position informations about base pair 2 (4 volumes)
//     |                     | Position informations about base pair 3 (4 volumes)
//     |                     | etc.
//     |
//     |nucleosome 2 --------| Position informations about base pair 1 (4 volumes)
//     |                     | Position informations about base pair 2 (4 volumes)
//     |                     | Position informations about base pair 3 (4 volumes)
//     |                     | etc.
//
typedef std::vector<std::vector<DNAPlacementData> > DNAPosData;

class GeoCalculationV2
{
public:
    GeoCalculationV2(G4int verbose=0, G4double factor=1.);

    ~GeoCalculationV2();

    //----------------------------------------------------------------------------------------------
    // Initialize the GeoCalculation object by setting up all the geometrical parameters.
    //----------------------------------------------------------------------------------------------
    void Initialize();

    //----------------------------------------------------------------------------------------------
    // Getters
    //----------------------------------------------------------------------------------------------
    std::vector<G4ThreeVector>* GetNucleosomePosition(){return fPosNucleo;}
    DNAPosData* GetAllDNAVolumePositions(){return fPosDNA;}
    std::vector<DNAPlacementData>* GetDNAVolumePositionsForNucleosome(G4int nucl){return &fPosDNA->at(nucl);}
    std::map<G4ThreeVector, G4double>* GetPosAndRadiusMap(){return fPosAndRadiusMap;}
    G4double GetSugarTHFRadiusWater(){return fSugarTHFRadiusWater;}
    G4double GetSugarTMPRadiusWater(){return fSugarTMPRadiusWater;}
    G4double GetSugarTHFRadius(){return fSugarTHFRadius;}
    G4double GetSugarTMPRadius(){return fSugarTMPRadius;}
    G4double GetBaseRadiusWater(){return fBaseRadiusWater;}
    G4double GetBaseRadius(){return fBaseRadius;}
    G4double GetFiberPitch(){return fFiberPitch;}
    G4double GetFiberDeltaAngle(){return fFiberDeltaAngle;}
    G4double GetFiberNbNuclPerTurn(){return fFiberNbNuclPerTurn;}
    G4double GetBpNum(){return fBpNum;}
    G4double GetHistoneHeight(){return fHistoneHeight;}
    G4double GetHistoneRadius(){return fHistoneRadius;}

    //----------------------------------------------------------------------------------------------
    // Setters
    //----------------------------------------------------------------------------------------------
    void SetSugarTHFRadiusWater(G4double radius){fSugarTHFRadiusWater=radius;}
    void SetSugarTMPRadiusWater(G4double radius){fSugarTMPRadiusWater=radius;}
    void SetSugarTHFRadius(G4double radius){fSugarTHFRadius=radius;}
    void SetSugarTMPRadius(G4double radius){fSugarTMPRadius=radius;}
    void SetBaseRadiusWater(G4double radius){fBaseRadiusWater=radius;}
    void SetBaseRadius(G4double radius){fBaseRadius=radius;}

private:
    G4int fVerbose;
    G4double fFactor;

    // Important numbers
    G4int fNucleoNum;
    G4int fBpNum;

    // Elementary parameters
    G4double fSugarTHFRadius;
    G4double fSugarTMPRadius;
    G4double fBaseRadius;
    G4double fSugarTHFRadiusWater;
    G4double fSugarTMPRadiusWater;
    G4double fBaseRadiusWater;

    G4ThreeVector fPosSugarTMP1;
    G4ThreeVector fPosSugarTHF1;
    G4ThreeVector fPosBase1;
    G4ThreeVector fPosBase2;
    G4ThreeVector fPosSugarTHF2;
    G4ThreeVector fPosSugarTMP2;

    // Calculation parameters
    // Fiber parameters
    G4int fHistoneNum;
    G4double fHistoneRadius;
    G4double fHistoneHeight;
    G4double fFiberPitch;
    G4double fFiberCentralRadius;
    G4double fFiberNbNuclPerTurn;
    G4double fFiberDeltaAngle;
    std::vector<std::vector<G4double> > fFiberHelixMatVect;
    std::vector<std::vector<std::vector<G4double> > > fRotFiberMatVect;

    // DNA around histone parameters
    // First helix parameters
    G4int bpNumAroundHistone;//const G4int bpNumAroundHistone = 154;
    G4double angleBpAroundHistone;

    // Second helix parameters (big simple helix)
    G4double secondHelixPitch;
    G4double centralRadius;
    G4double nbBasePairPerTurn;
    G4double deltaAngle ;

    // DNA linker parameters
    G4int bpNumForLinker;//const G4int bpNumForLinker = 46;
    G4double linkerCentralRadius;
    G4double linkerHeightPerBp;
    G4double linkerArcCircleMat[3];
    G4double nbBasePairPerTurnForLinker;
    G4double deltaLinkerAngle;

    // Pos containers
    std::vector<G4ThreeVector>* fPosNucleo;
    DNAPosData* fPosDNA;
    std::map<G4ThreeVector, G4double>* fPosAndRadiusMap;

    //**********************************************************************************************
    // Methods
    //**********************************************************************************************

    //----------------------------------------------------------------------------------------------
    // Calculate positions of nucleotide volumes around the basis hisotone volumes (3). Fills and 
    // returns pointer to DNAPosData object, which is a 2D vector of DNAPlacementData objects.
    //      outer index i : spans nucleosomes (1 to 3)
    //      inner index j : spans all bp in that nucleosome
    // DNAPlacementData is a structure containing 7 G4ThreeVectors containing the coordinates of 6
    // residues & central location of a given nucleotide base pair.
    //----------------------------------------------------------------------------------------------
    DNAPosData* CalculateDNAPosition(G4int histoneNum, G4ThreeVector& posSugarTMP1, 
                                     G4ThreeVector &posSugarTHF1, G4ThreeVector &posBase1,
                                     G4ThreeVector &posBase2, G4ThreeVector &posSugarTHF2, 
                                     G4ThreeVector &posSugarTMP2);

    //----------------------------------------------------------------------------------------------
    // Calculate nucleosome positions for basis nucleosomes. Remaining nucleosome positions can be 
    // established by translating & rotating these three basis positions. Note the number of basis
    // nucleosomes (i.e. nucleoNum) is 3 for typical use. Returns a pointer to a vector of size
    // nucleoNum (3), with each element containing a G4ThreeVector of the spatial coordinates for 1
    // of the basis nucleosomes.
    //----------------------------------------------------------------------------------------------
    std::vector<G4ThreeVector>* CalculateNucleosomePosition(G4int nucleoNum);

    //----------------------------------------------------------------------------------------------
    // Create a map of the 6 volumes comprising nucleotide base pairs, across all of the basis 
    // nucleosomes. Key = G4ThreeVector of coordinates, Value = radius of the hydration shell.
    // Input parameter: pointer to a DNAPosData object created by CalculateDNAPosition().
    // Output: a map of coordinates:radius pairs that is used by 
    //     GeoVolume::CreateNucleosomeCuttedSolidsAndLogicals()
    // NOTE: it is the hydration shell radius that is used.
    // Map size = 3600 (3 nucleosomes x 200 bp/nucl x 6 volumes/bp)
    //----------------------------------------------------------------------------------------------
    std::map<G4ThreeVector, G4double>* GenerateCoordAndRadiusMap(DNAPosData *dnaPosData);

    //----------------------------------------------------------------------------------------------
    // Helper function used by CalculateDNAPosition().
    //----------------------------------------------------------------------------------------------
    G4double GetAngleToXAxis(G4ThreeVector t);
};

#endif // GEOCALCULATIONV2_HH
