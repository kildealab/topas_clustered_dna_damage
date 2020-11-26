//**************************************************************************************************
// This class handles the construction of a logical volume for a chromatin fiber using the method
// BuildLogicFiber(). Additional methods provide supporting functionality for this method (e.g. the
// generation of physical volumes and logical volumes for residues and performing the necesary
// cutting procedures to prevent geometrical overlaps).
//**************************************************************************************************

#ifndef GEOVOLUMEV2_HH
#define GEOVOLUMEV2_HH

#include <map>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

struct DNAPlacementData;

class GeoVolumeV2
{
public:
    GeoVolumeV2(G4int verbose=0, G4double factor=1.);

    ~GeoVolumeV2();

    //----------------------------------------------------------------------------------------------
    // Getters
    //----------------------------------------------------------------------------------------------
    G4double GetSugarTHFRadiusWater(){return fSugarTHFRadiusWater;}
    G4double GetSugarTMPRadiusWater(){return fSugarTMPRadiusWater;}
    G4double GetSugarTHFRadius(){return fSugarTHFRadius;}
    G4double GetSugarTMPRadius(){return fSugarTMPRadius;}
    G4double GetBaseRadiusWater(){return fBaseRadiusWater;}
    G4double GetBaseRadius(){return fBaseRadius;}
    std::map<G4String, std::vector<std::vector<double> > >* GetDNAMoleculesPositions(){return fpDnaMoleculePositions;}

    //----------------------------------------------------------------------------------------------
    // Setters
    //----------------------------------------------------------------------------------------------
    void SetSugarTHFRadiusWater(G4double radius){fSugarTHFRadiusWater=radius;}
    void SetSugarTMPRadiusWater(G4double radius){fSugarTMPRadiusWater=radius;}
    void SetSugarTHFRadius(G4double radius){fSugarTHFRadius=radius;}
    void SetSugarTMPRadius(G4double radius){fSugarTMPRadius=radius;}
    void SetBaseRadiusWater(G4double radius){fBaseRadiusWater=radius;}
    void SetBaseRadius(G4double radius){fBaseRadius=radius;}
    void SetFiberPitch(G4double fiberPitch){fFiberPitch=fiberPitch;}
    void SetFiberNbNuclPerTurn(G4double nbNuclPerTurn){fFiberNbNuclPerTurn=nbNuclPerTurn;}
    void SetFiberDeltaAngle(G4double fiberDeltaAngle){fFiberDeltaAngle=fiberDeltaAngle;}
    void SetNucleoNum(G4double nucleoNum){fNucleoNum=nucleoNum;}
    void SetBpNum(G4double bpNum){fBpNum=bpNum;}
    void SetHistoneHeight(G4double histoneHeight){fHistoneHeight=histoneHeight;}
    void SetHistoneRadius(G4double histoneRadius){fHistoneRadius=histoneRadius;}

    //----------------------------------------------------------------------------------------------
    // Create and return a logical volume for a chromatin fiber.  Within this logical volume are the
    // physical volumes for the histones, the resiudes and their hydration shells. Solids and
    // logicals are generated for the histones within this metohd directly, whereas those for the
    // residues are generated using CreateNucleosomeCuttedSolidsAndLogicals(). A map,
    // fpDnaMoleculePositions, containing the coordinates for all residues and histones is filled
    // and can be accessed using GetDNAMoleculesPositions().
    //----------------------------------------------------------------------------------------------
    G4LogicalVolume *BuildLogicFiber(std::vector<std::vector<DNAPlacementData> > *dnaVolPos,
                                     std::vector<G4ThreeVector> *posNucleo,
                                     std::map<G4ThreeVector, G4double> *posAndRadiusMap,
                                     G4bool isVisu=false);
    
private:
    G4int fVerbose;
    G4double fFactor;
    G4Material* fWater;

    // Elementary parameters
    G4double fSugarTHFRadiusWater;
    G4double fSugarTMPRadiusWater;
    G4double fSugarTHFRadius;
    G4double fSugarTMPRadius;
    G4double fBaseRadiusWater;
    G4double fBaseRadius;

    G4double fFiberPitch;
    G4double fFiberNbNuclPerTurn;
    G4double fFiberDeltaAngle;

    G4double fHistoneHeight;
    G4double fHistoneRadius;

    G4int fNucleoNum;
    G4int fBpNum;

    // moleculeName: x, y, z, copyNumber, strand
    std::map<G4String, std::vector<std::vector<double> > >* fpDnaMoleculePositions;

    //----------------------------------------------------------------------------------------------
    // Create the solid and logical volumes required to build DNA around one histone.
    // Return a map as:
    // Key: name of the volume (base1, base2, base1Water, ...). Size = 12.
    // Content: vector of corresponding logical volumes (each vector size = 200)
    //----------------------------------------------------------------------------------------------
    std::map<G4String, std::vector<G4LogicalVolume *> >* CreateNucleosomeCuttedSolidsAndLogicals(
        std::vector<DNAPlacementData> *nucleosomeVolumePositions, 
        std::map<G4ThreeVector, G4double> *posAndRadiusMap, G4bool isVisu=false);

    //----------------------------------------------------------------------------------------------
    // Cut algorithm to avoid overlaps. Idea: we must have a reference and a target. The reference
    // is the solid we are considering ( described by parameters solidOrbRef & posRef) and which
    // could be cut if an overlap is detected with the target solid. In a geometry, it implies we
    // have to go through all the target  solids for a given reference solid. Target solid info
    // (position and radius) is included in tarMap This method will return the cut spherical
    // reference solid. Note: fifth parameter, i.e. G4bool in doesn't seem to be used.
    //----------------------------------------------------------------------------------------------
    G4VSolid *CreateCutSolid(G4Orb *solidOrbRef,
                             G4ThreeVector& posRef,
                             std::map<G4ThreeVector, G4double> *tarMap,
                             G4String volName = "",
                             G4bool in = false);

    //----------------------------------------------------------------------------------------------
    // Helper function used by buildLogicFiber() in verbose mode (i.e. if want to output the mean
    // volume of cut residue and water volumes).
    //----------------------------------------------------------------------------------------------
    void CalculateMeanVol(std::map<G4String, std::vector<G4LogicalVolume *> > *logicSolidsMap);

};

#endif // GEOVOLUMEV2_HH
