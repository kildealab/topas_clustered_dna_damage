//**************************************************************************************************
// Author: Logan Montgomery
// Based on code written by:
//      - J Schuemann et al. (2019). DOI:10.1667/RR15226.1
//      - A McNamara et al. (2018). DOI:10.1088/1361-6560/aad8eb
//      - S Meylan et al. (2017). DOI:10.1038/s41598-017-11851-4
//
// This class is a custom Topas geometry component that creates a nucleus of DNA that is comprised
// of one or more voxels, each of which contains at least one chromatin fiber.
//**************************************************************************************************

#ifndef VoxelizedNuclearDNA_hh
#define VoxelizedNuclearDNA_hh

#include "TsVGeometryComponent.hh"

#include <map>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Orb.hh"
#include "GeoCalculationV2.hh"



struct DNAPlacementData;

class VoxelizedNuclearDNA : public TsVGeometryComponent
{
public:
    //----------------------------------------------------------------------------------------------
    // Constructor. Initialize member variables using a variety of methods.
    //----------------------------------------------------------------------------------------------
	VoxelizedNuclearDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);

    //----------------------------------------------------------------------------------------------
    // Destructor
    //----------------------------------------------------------------------------------------------
	~VoxelizedNuclearDNA();

    //----------------------------------------------------------------------------------------------
    // Construct the geometry.
    //----------------------------------------------------------------------------------------------
	G4VPhysicalVolume* Construct();

    //----------------------------------------------------------------------------------------------
    // Read in parameters from the Topas parameter file & save values in some member variables.
    //----------------------------------------------------------------------------------------------
    void ResolveParameters();

private:
    //----------------------------------------------------------------------------------------------
    // Create and return a logical volume for a chromatin fiber.
    //----------------------------------------------------------------------------------------------
    G4LogicalVolume *BuildLogicFiber(std::vector<std::vector<DNAPlacementData> > *dnaVolPos,
                                     std::vector<G4ThreeVector> *posNucleo,
                                     std::map<G4ThreeVector, G4double> *posAndRadiusMap);

    //----------------------------------------------------------------------------------------------
    // Create the solid and logical volumes required to build DNA around one histone.
    // Return a map as:
    // Key: name of the volume (base1, base2, base1Water, ...). Size = 12.
    // Content: vector of corresponding logical volumes (each vector size = 200)
    //----------------------------------------------------------------------------------------------
    std::map<G4String, std::vector<G4LogicalVolume *> >* CreateNucleosomeCuttedSolidsAndLogicals(
        std::vector<DNAPlacementData> *nucleosomeVolumePositions,
        std::map<G4ThreeVector, G4double> *posAndRadiusMap);

    //----------------------------------------------------------------------------------------------
    // Algorithm for cutting DNA residue solids to avoid overlaps. Return the cut spherical solid.
    //----------------------------------------------------------------------------------------------
    G4VSolid *CreateCutSolid(G4Orb *solidOrbRef,
                             G4ThreeVector& posRef,
                             std::map<G4ThreeVector, G4double> *tarMap,
                             G4String volName = "");

    //----------------------------------------------------------------------------------------------
    // Arrange identical DNA fibers in a cubic voxel. Return the logical volume of that voxel.
    //----------------------------------------------------------------------------------------------
    G4LogicalVolume *ConstructLogicalVoxel(G4LogicalVolume* logicalFiber);

    //----------------------------------------------------------------------------------------------
    // Helper function to detect if geometrical overlap & throw error if so.
    //----------------------------------------------------------------------------------------------
    void ThrowOverlapError();

    //----------------------------------------------------------------------------------------------
    // Member variables
    //----------------------------------------------------------------------------------------------
    G4double fWrapperRadius;
    G4double fWrapperHeight;

    G4bool fCutVolumes;
    G4bool fCheckForOverlaps;
    G4int fOverlapsResolution;
    G4bool fQuitIfOverlap;

    G4bool fFillFibersWithDNA;

    G4bool fBuildNucleus;
    G4int fNumVoxelsPerSide;
    G4double fVoxelSideLength;

    G4VPhysicalVolume* pFiber;

    G4bool fUseG4Volumes;

    GeoCalculationV2* fGeoCalculation;

    G4int fNumNucleosomePerFiber;
    G4int fNumBpPerNucleosome;

    G4double fFiberRadius;
    G4double fFiberHalfLength;

    G4double fFiberPitch;
    G4double fFiberNbNuclPerTurn;
    G4double fFiberDeltaAngle;

    G4double fHistoneHeight;
    G4double fHistoneRadius;

    G4double fSugarTHFRadiusWater;
    G4double fSugarTMPRadiusWater;
    G4double fSugarTHFRadius;
    G4double fSugarTMPRadius;
    G4double fBaseRadiusWater;
    G4double fBaseRadius;

    G4String fWaterName;
    G4Material* fWater;
    G4String fDNAMaterialName;
    G4Material* fDNAMaterial;
		G4String fHistoneMaterialName;
    G4Material* fHistoneMaterial;

    // This map is indexed as moleculeName: <x, y, z, copyNumber, strand>
    std::map<G4String, std::vector<std::vector<double> > >* fpDnaMoleculePositions;
};

#endif
