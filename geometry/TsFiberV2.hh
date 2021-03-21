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

#ifndef TsFiberV2_hh
#define TsFiberV2_hh

#include "TsVGeometryComponent.hh"
#include "GeoManagerV2.hh"

#include <map>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"


struct DNAPlacementData;

class TsFiberV2 : public TsVGeometryComponent
{    
public:
	TsFiberV2(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsFiberV2();
    
    GeoManagerV2* fGeoManager;
    GeoCalculationV2* fGeoCalculation;

	G4VPhysicalVolume* Construct();    
    //void SetMaterial (G4String materialChoice);
    
    void UpdateGeometry();
    
private:
    
    G4double fWrapperRadius;
    G4double fWrapperHeight;
    G4bool fCutVolumes;
    G4bool fCheckForOverlaps;
    G4int fOverlapsResolution;
    G4bool fQuitIfOverlap;
    G4String fDNAMaterialName;
    G4Material* fDNAMaterial;

    G4bool fBuildBases;
    G4bool fBuildNucleus;
    G4int fNumVoxelsPerSide;
    G4double fVoxelSideLength;

    // GeoVolume parameters
    G4VPhysicalVolume* pFiber;

    G4bool fUseG4Volumes;

    G4int fVerbose;
    G4double fFactor;
    G4double fFiberRadius;
    G4double fFiberHalfLength;

    G4double fSugarTHFRadiusWater;
    G4double fSugarTMPRadiusWater;
    G4double fSugarTHFRadius;
    G4double fSugarTMPRadius;
    G4double fBaseRadiusWater;
    G4double fBaseRadius;

    G4double fHistoneHeight;
    G4double fHistoneRadius;

    G4int fNucleoNum;
    G4int fBpNum;

    G4double fFiberPitch;
    G4double fFiberNbNuclPerTurn;
    G4double fFiberDeltaAngle;

    G4String fWaterName;
    G4Material* fWater;

    // moleculeName: x, y, z, copyNumber, strand
    std::map<G4String, std::vector<std::vector<double> > >* fpDnaMoleculePositions;

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
                                     G4bool cutVolumes=true,
                                     G4bool checkForOverlaps=true,
                                     G4int overlapsResolution=100,
                                     G4bool quitIfOverlap=true);

    //----------------------------------------------------------------------------------------------
    // Create the solid and logical volumes required to build DNA around one histone.
    // Return a map as:
    // Key: name of the volume (base1, base2, base1Water, ...). Size = 12.
    // Content: vector of corresponding logical volumes (each vector size = 200)
    //----------------------------------------------------------------------------------------------
    std::map<G4String, std::vector<G4LogicalVolume *> >* CreateNucleosomeCuttedSolidsAndLogicals(
        std::vector<DNAPlacementData> *nucleosomeVolumePositions, 
        std::map<G4ThreeVector, G4double> *posAndRadiusMap, G4bool cutVolumes=true);

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
    //----------------------------------------------------------------------------------------------
    G4LogicalVolume *ConstructLogicalVoxel(G4LogicalVolume* logicalFiber);

    void ThrowOverlapError();
};

#endif
