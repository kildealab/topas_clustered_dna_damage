// Component for VoxelizedNuclearDNA
//
//**************************************************************************************************
// Author: Logan Montgomery, James Manalad
// Based on code written by:
//      - J Schuemann et al. (2019). DOI:10.1667/RR15226.1
//      - A McNamara et al. (2018). DOI:10.1088/1361-6560/aad8eb
//      - S Meylan et al. (2017). DOI:10.1038/s41598-017-11851-4
//
// This class is a custom Topas geometry component that creates a nucleus of DNA that is comprised
// of one or more voxels, each of which contains at least one chromatin fiber.
//**************************************************************************************************

#include "VoxelizedNuclearDNA.hh"
#include "GeoCalculationV2.hh"

#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4UserLimits.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4ThreeVector.hh"

#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Exception.hh"
#include <chrono>


//--------------------------------------------------------------------------------------------------
// Structure containting coordinates of volumes in a DNA base pair.
// Identical to the definition in GeoCalculation.
//--------------------------------------------------------------------------------------------------
struct DNAPlacementData
{
    G4ThreeVector posCenterDNA;
    G4ThreeVector posSugarTMP1;
    G4ThreeVector posSugarTHF1;
    G4ThreeVector posBase1;
    G4ThreeVector posBase2;
    G4ThreeVector posSugarTHF2;
    G4ThreeVector posSugarTMP2;
};


//--------------------------------------------------------------------------------------------------
// Constructor. Initialize member variables using a variety of methods.
//--------------------------------------------------------------------------------------------------
VoxelizedNuclearDNA::VoxelizedNuclearDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters(); // initialize some member variables using Topas parameter file

    fGeoCalculation = new GeoCalculationV2(0, 1.);

    fpDnaMoleculePositions = new std::map<G4String, std::vector<std::vector<G4double> > >();

    //----------------------------------------------------------------------------------------------
    // A GeoCalculation object is used set various parameters for the configuration of DNA content
    // in single chromatin fiber.
    //----------------------------------------------------------------------------------------------
    fFiberRadius = 17.*nm;
    fFiberHalfLength = 68.*nm;

    fGeoCalculation->Initialize();

    fBaseRadius = fGeoCalculation->GetBaseRadius();
    fBaseRadiusWater = fGeoCalculation->GetBaseRadiusWater();
    fSugarTMPRadius = fGeoCalculation->GetSugarTMPRadius();
    fSugarTHFRadius = fGeoCalculation->GetSugarTHFRadius();
    fSugarTMPRadiusWater = fGeoCalculation->GetSugarTMPRadiusWater();
    fSugarTHFRadiusWater = fGeoCalculation->GetSugarTHFRadiusWater();

    fHistoneHeight = fGeoCalculation->GetHistoneHeight() ;
    fHistoneRadius = fGeoCalculation->GetHistoneRadius();
    fFiberPitch = fGeoCalculation->GetFiberPitch();
    fFiberDeltaAngle = fGeoCalculation->GetFiberDeltaAngle() ;
    fFiberNbNuclPerTurn = fGeoCalculation->GetFiberNbNuclPerTurn();

    //----------------------------------------------------------------------------------------------
    // Create modified water material to be used in DNA volumes (used to identify volumes in which
    // to score damage)
    //----------------------------------------------------------------------------------------------
    fWaterName = "G4_WATER";
    fWater = GetMaterial(fWaterName);

    fDNAMaterial = GetMaterial(fDNAMaterialName);
    fHistoneMaterial = GetMaterial(fHistoneMaterialName);
}

//--------------------------------------------------------------------------------------------------
// Destructor
//--------------------------------------------------------------------------------------------------
VoxelizedNuclearDNA::~VoxelizedNuclearDNA()
{
     delete fGeoCalculation;

     delete fpDnaMoleculePositions;
}


//--------------------------------------------------------------------------------------------------
// Read in parameters from the Topas parameter file & save values in some member variables. Default
// values are provided in case the user does not specify a value.
//--------------------------------------------------------------------------------------------------
void VoxelizedNuclearDNA::ResolveParameters() {
    if (fPm->ParameterExists(GetFullParmName("DNANumBpPerNucleosome")))
        fNumBpPerNucleosome = fPm->GetIntegerParameter(GetFullParmName("DNANumBpPerNucleosome"));
    else
        fNumBpPerNucleosome = 200;

    if (fPm->ParameterExists(GetFullParmName("DnaNumNucleosomePerFiber")))
        fNumNucleosomePerFiber = fPm->GetIntegerParameter(GetFullParmName("DnaNumNucleosomePerFiber"));
    else
        fNumNucleosomePerFiber = 90;

    if (fPm->ParameterExists(GetFullParmName("CutVolumes")))
        fCutVolumes = fPm->GetBooleanParameter(GetFullParmName("CutVolumes"));
    else
        fCutVolumes = true;

    if (fPm->ParameterExists(GetFullParmName("UseG4Volumes")))
        fUseG4Volumes = fPm->GetBooleanParameter(GetFullParmName("UseG4Volumes"));
    else
        fUseG4Volumes = false;

    if (fPm->ParameterExists(GetFullParmName("FillFibersWithDNA")))
        fFillFibersWithDNA = fPm->GetBooleanParameter(GetFullParmName("FillFibersWithDNA"));
    else
        fFillFibersWithDNA = true;

    if (fPm->ParameterExists(GetFullParmName("BuildNucleus")))
        fBuildNucleus = fPm->GetBooleanParameter(GetFullParmName("BuildNucleus"));
    else
        fBuildNucleus = true;

    if (fPm->ParameterExists(GetFullParmName("NumVoxelsPerSide")))
        fNumVoxelsPerSide = fPm->GetIntegerParameter(GetFullParmName("NumVoxelsPerSide"));
    else
        fNumVoxelsPerSide = 1;

    if (fPm->ParameterExists(GetFullParmName("VoxelSideLength")))
        fVoxelSideLength = fPm->GetDoubleParameter(GetFullParmName("VoxelSideLength"),"Length");
    else
        fVoxelSideLength = 250*nm;

    if (fPm->ParameterExists(GetFullParmName("DNAMaterialName")))
        fDNAMaterialName = fPm->GetStringParameter(GetFullParmName("DNAMaterialName"));
    else
        fDNAMaterialName = "G4_WATER_DNA";
    if (fPm->ParameterExists(GetFullParmName("HistoneMaterialName")))
        fHistoneMaterialName = fPm->GetStringParameter(GetFullParmName("HistoneMaterialName"));
    else
        fHistoneMaterialName = "G4_WATER_HISTONE";

    // Parameters not specific to this extension. Can't/don't need to use GetFullParmName()
    fCheckForOverlaps = fPm->GetBooleanParameter("Ge/CheckForOverlaps");
    fOverlapsResolution = fPm->GetIntegerParameter("Ge/CheckForOverlapsResolution");
    fQuitIfOverlap = fPm->GetBooleanParameter("Ge/QuitIfOverlapDetected");
}


//--------------------------------------------------------------------------------------------------
// Construct the geometry. "Envelope" volumes are the outer shell of this custom geometry component.
// Fill the envelope with some amount/configuration of DNA.
//--------------------------------------------------------------------------------------------------
G4VPhysicalVolume* VoxelizedNuclearDNA::Construct()
{
	BeginConstruction();

    //----------------------------------------------------------------------------------------------
    // Construct the envelope (wrapper) volume.
    //----------------------------------------------------------------------------------------------
    G4double envelopeSideLength = fNumVoxelsPerSide * fVoxelSideLength;
    G4Box* sWrapper = new G4Box("solid_wrapper", envelopeSideLength, envelopeSideLength, envelopeSideLength);
    fEnvelopeLog = CreateLogicalVolume(sWrapper);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

    //----------------------------------------------------------------------------------------------
    // Construct the logical volume for a single chromatin fiber.
    //----------------------------------------------------------------------------------------------
    G4LogicalVolume* lFiber = BuildLogicFiber(fGeoCalculation->GetAllDNAVolumePositions(),
        fGeoCalculation->GetNucleosomePosition(),
        fGeoCalculation->GetPosAndRadiusMap());

    //----------------------------------------------------------------------------------------------
    // Construct physical volume for the DNA. Either a voxelized nucleus containing many fibers or a
    // single fiber. Place in the outermost physical volume for this custom component (i.e. in
    // fEnvelopePhys).
    //----------------------------------------------------------------------------------------------
    if (fBuildNucleus) {
        G4LogicalVolume* voxelLogical = ConstructLogicalVoxel(lFiber);

        //------------------------------------------------------------------------------------------
        // Generate cubic arrangement of voxels using replica volumes. When doing nested replica
        // volumes, start construction with the outermost replicated direction and "reserve" the
        // space. Then fill that space with nested replicas as shown below.
        //------------------------------------------------------------------------------------------
        G4double nucleusSideLength = fVoxelSideLength * fNumVoxelsPerSide;

        // Outermost dimension is Y. Reserve 3D space for replicated 2D arrays of voxels.
        G4Box* solidEmptyVoxelArea = new G4Box("solid_voxel_area", nucleusSideLength, fVoxelSideLength, nucleusSideLength);
        G4LogicalVolume* logEmptyVoxelArea = CreateLogicalVolume("Voxel",fWaterName,solidEmptyVoxelArea);
        G4VPhysicalVolume* voxels_3d = CreatePhysicalVolume("ReplicaVoxels3D",logEmptyVoxelArea,fEnvelopePhys,kYAxis,fNumVoxelsPerSide,2*fVoxelSideLength);

        // Middle dimension is X. Reserve 2D space for replicated 1D arrays of voxels in xz-plane.
        G4Box* solidEmptyVoxelLength = new G4Box("solid_voxel_length", fVoxelSideLength, fVoxelSideLength, nucleusSideLength);
        G4LogicalVolume* logEmptyVoxelLength = CreateLogicalVolume("Voxel",fWaterName,solidEmptyVoxelLength);
        G4VPhysicalVolume* voxels_2d = CreatePhysicalVolume("ReplicaVoxels2D",logEmptyVoxelLength,voxels_3d,kXAxis,fNumVoxelsPerSide,2*fVoxelSideLength);

        // Innermost dimension is Z. Fill 1D array of voxels in z-dimension.
        G4VPhysicalVolume* voxels_1d = CreatePhysicalVolume("ReplicaVoxels1D",voxelLogical,voxels_2d,kZAxis,fNumVoxelsPerSide,2*fVoxelSideLength);
    }
    // Or place single fiber physical volume
    else {
        if (fUseG4Volumes){
            G4ThreeVector emptyPlacement = G4ThreeVector(0.,0.,0.);
            pFiber = new G4PVPlacement(0,emptyPlacement,lFiber,"Fiber",fEnvelopeLog,false,0);
        }
        else{
            pFiber = CreatePhysicalVolume("Fiber", lFiber, fEnvelopePhys);
        }
    }

	return fEnvelopePhys;
}


//--------------------------------------------------------------------------------------------------
// Create and return a logical volume for a chromatin fiber.
// Within this logical volume are the physical volumes for the histones, the resiudes and their
// hydration shells. Solids and logicals are generated for the histones within this metohd directly,
// whereas those for the residues are generated using CreateNucleosomeCuttedSolidsAndLogicals().
// A map, fpDnaMoleculePositions, containing the coordinates for all residues and histones is filled
// and can be accessed using GetDNAMoleculesPositions().
//--------------------------------------------------------------------------------------------------
G4LogicalVolume* VoxelizedNuclearDNA::BuildLogicFiber(std::vector<std::vector<DNAPlacementData> >* dnaVolPos,
                                            std::vector<G4ThreeVector>* posNucleo,
                                            std::map<G4ThreeVector, G4double>* posAndRadiusMap)
{
    //----------------------------------------------------------------------------------------------
    // Throw error if any of these member variables haven't been initialized correctly.
    //----------------------------------------------------------------------------------------------
    if(fSugarTHFRadius==-1||fSugarTMPRadius==-1||fBaseRadius==-1||fFiberPitch==-1||fFiberNbNuclPerTurn==-1
            ||fNumNucleosomePerFiber==-1||fNumBpPerNucleosome==-1||fHistoneHeight==-1||fHistoneRadius==-1)
    {
        G4cerr<<"FatalError: VoxelizedNuclearDNA::BuildLogicFiber. A class parameter has not been "
            << "initialized and its value is still negative"<<G4endl;
        std::exit(EXIT_FAILURE);
    }
    //----------------------------------------------------------------------------------------------
    // Create cylindrical fiber volume
    //----------------------------------------------------------------------------------------------
    G4Tubs* solidFiber = new G4Tubs("solid_fiber", 0., fFiberRadius, fFiberHalfLength, 0, 360);

    G4LogicalVolume* logicFiber;
    if (fUseG4Volumes) {
        logicFiber = new G4LogicalVolume(solidFiber,fWater,"Fiber");
    }
    else {
        logicFiber = CreateLogicalVolume("Fiber",fWaterName,solidFiber);
    }

    // If not building DNA in the fiber, return now with empty fiber (for visualization of large geometries)
    if (!fFillFibersWithDNA) {
        return logicFiber;
    }

    //----------------------------------------------------------------------------------------------
    // Create the histone volume
    //----------------------------------------------------------------------------------------------
    G4Tubs* solidHistone = new G4Tubs("solid histone", 0., fHistoneRadius, fHistoneHeight, 0, 360);

    G4LogicalVolume* logicHistone;
    if (fUseG4Volumes) {
        logicHistone = new G4LogicalVolume(solidHistone,fHistoneMaterial,"Histone");
    }
    else {
        logicHistone = CreateLogicalVolume("Histone",fHistoneMaterialName,solidHistone);
    }

    //----------------------------------------------------------------------------------------------
    // Generate logical volumes for the nucleotide base pairs.
    //----------------------------------------------------------------------------------------------
    // For the positions, only use the nucleosome #2 (index=1). It is a "middle" nucleosome and
    // thus, the two extremities will be cut to allow for  proper linking of one nucleosome to the
    // next. Note dnaVolPos is essentially the output of GeoCalculation's CalculateDNAPosition()
    // method.  I.e. all nucleotide positions around 3 basis histone complexes.
    std::vector<DNAPlacementData>* nuclVolPos = &dnaVolPos->at(1);

    // Create all the DNA volumes (solid & logical) around the histone based on nucleosome #2
    // (index=1) positions. Place this nucleosome several times to build the fiber. This is done to
    // save memory and improve speed. Logical volumes are saved a map (key = name of the volume
    // [e.g. sugar1], value = vector of corresponding logical volumes). Note posAndRadiusMap is the
    // output of GeoCalculation's GenerateCoordAndRadiusMap() method. I.e. a map of radii for 6
    // residue volumes in each of 200 bp in each of 3 basis nucleosomes (3600 volumes)
    std::map<G4String, std::vector<G4LogicalVolume*> >* volMap
            = CreateNucleosomeCuttedSolidsAndLogicals(nuclVolPos, posAndRadiusMap);
    // The resulting volMap is indexed by one of 12 entries (6 residues & 6 hydration shells). Each
    // entry has 200 elements, each corresponding to a distinct logical volume

    //----------------------------------------------------------------------------------------------
    // Save the positions of all residue volumes in nucleosome #2 in vectors
    //----------------------------------------------------------------------------------------------
    G4ThreeVector posSugarTMP1;
    G4ThreeVector posSugarTHF1;
    G4ThreeVector posBase1;
    G4ThreeVector posBase2;
    G4ThreeVector posSugarTHF2;
    G4ThreeVector posSugarTMP2;
    std::vector<G4ThreeVector> posSugarTMP1Vect;
    std::vector<G4ThreeVector> posSugarTHF1Vect;
    std::vector<G4ThreeVector> posBase1Vect;
    std::vector<G4ThreeVector> posBase2Vect;
    std::vector<G4ThreeVector> posSugarTHF2Vect;
    std::vector<G4ThreeVector> posSugarTMP2Vect;

    // Get the bp volume positions of the middle basis nucleosome (which is used to generate cut
    // solids) and save in a vector to be accessed & rotated in the following forloop
    for(int j=0;j<fNumBpPerNucleosome;++j)
    {
        posSugarTMP1 = dnaVolPos->at(1)[j].posSugarTMP1;
        posSugarTHF1 = dnaVolPos->at(1)[j].posSugarTHF1;
        posBase1 = dnaVolPos->at(1)[j].posBase1;
        posBase2 = dnaVolPos->at(1)[j].posBase2;
        posSugarTHF2 = dnaVolPos->at(1)[j].posSugarTHF2;
        posSugarTMP2 = dnaVolPos->at(1)[j].posSugarTMP2;

        posSugarTMP1Vect.push_back(posSugarTMP1);
        posSugarTHF1Vect.push_back(posSugarTHF1);
        posBase1Vect.push_back(posBase1);
        posBase2Vect.push_back(posBase2);
        posSugarTHF2Vect.push_back(posSugarTHF2);
        posSugarTMP2Vect.push_back(posSugarTMP2);
    }

    //----------------------------------------------------------------------------------------------
    // Save the first histone position
    //----------------------------------------------------------------------------------------------
    G4ThreeVector posHistone = posNucleo->at(0);

    //----------------------------------------------------------------------------------------------
    // Distance offset used to start placing volumes at one end of the fiber
    //----------------------------------------------------------------------------------------------
    G4ThreeVector minusForFiber = G4ThreeVector(0.,0.,-solidFiber->GetDz() + fHistoneHeight);

    //----------------------------------------------------------------------------------------------
    // z shift to be applied when positioning subsequent nucleosomes
    //----------------------------------------------------------------------------------------------
    G4double zShift = fFiberPitch/fFiberNbNuclPerTurn;
    G4int count = 0;

    //----------------------------------------------------------------------------------------------
    // Fill the chromatin fiber with DNA by iterating over each nucleosome & each nucleotide base
    // pair within each nucleosome. Create physical volumes using the already-created logical
    // volumes.
    //----------------------------------------------------------------------------------------------
    for(int i=0;i<fNumNucleosomePerFiber;++i)
    {
        // Rotate the nucleosome itself with respect to the z-axis, in order to align the cut
        // volumes appropriately to prevent overlaps. This is not the same as placing a nucleosome
        // at the next position around the fiber. That is done below. Our basis nucleosome is
        // nucleosome #2 (index=1). The following rotation logic accounts for this.
        // This rotation object will be applied to every physical volume placement below, in order
        // to align the cut bp volumes and prevent overlaps. This is a rotation about the
        // volume's own z-axis.
        G4RotationMatrix* rotCuts = new G4RotationMatrix();
        rotCuts->rotateZ((i-1)*-fFiberDeltaAngle);

        //------------------------------------------------------------------------------------------
        // Iterate over all bp in a nucleosome. At each iteration, generate physical volumes for all
        // residues in the current bp
        //------------------------------------------------------------------------------------------
        for(int j=0;j<fNumBpPerNucleosome;++j)
        {
            //--------------------------------------------------------------------------------------
            // The following positional rotations account for the fact that nucleosomes are
            // assembled in a helical fashion around the fibre axis (z). Note this has nothing to
            // do with rotations of the bp volumes themselves (which are handled by rotCuts obj).
            //--------------------------------------------------------------------------------------
            if(i==0)
            {
                posSugarTMP1 = posSugarTMP1Vect[j].rotateZ(-fFiberDeltaAngle);
                posSugarTHF1 = posSugarTHF1Vect[j].rotateZ(-fFiberDeltaAngle);
                posBase1 = posBase1Vect[j].rotateZ(-fFiberDeltaAngle);
                posBase2 = posBase2Vect[j].rotateZ(-fFiberDeltaAngle);
                posSugarTHF2 = posSugarTHF2Vect[j].rotateZ(-fFiberDeltaAngle);
                posSugarTMP2 = posSugarTMP2Vect[j].rotateZ(-fFiberDeltaAngle);
            }
            else
            {
                posSugarTMP1 = posSugarTMP1Vect[j].rotateZ(fFiberDeltaAngle);
                posSugarTHF1 = posSugarTHF1Vect[j].rotateZ(fFiberDeltaAngle);
                posBase1 = posBase1Vect[j].rotateZ(fFiberDeltaAngle);
                posBase2 = posBase2Vect[j].rotateZ(fFiberDeltaAngle);
                posSugarTHF2 = posSugarTHF2Vect[j].rotateZ(fFiberDeltaAngle);
                posSugarTMP2 = posSugarTMP2Vect[j].rotateZ(fFiberDeltaAngle);
            }

            //--------------------------------------------------------------------------------------
            // Apply a z shift, specific to each nucleosome, to span length of fiber
            //--------------------------------------------------------------------------------------
            posSugarTMP1 += G4ThreeVector(0.,0.,i*zShift-zShift);
            posSugarTHF1 += G4ThreeVector(0.,0.,i*zShift-zShift);
            posBase1 += G4ThreeVector(0.,0.,i*zShift-zShift);
            posBase2 += G4ThreeVector(0.,0.,i*zShift-zShift);
            posSugarTHF2 += G4ThreeVector(0.,0.,i*zShift-zShift);
            posSugarTMP2 += G4ThreeVector(0.,0.,i*zShift-zShift);

            //--------------------------------------------------------------------------------------
            // Apply shift such that fiber helix construction begins at at one end.
            //--------------------------------------------------------------------------------------
            posSugarTMP1 += minusForFiber;
            posSugarTHF1 += minusForFiber;
            posBase1 += minusForFiber;
            posBase2 += minusForFiber;
            posSugarTHF2 += minusForFiber;
            posSugarTMP2 += minusForFiber;

            //--------------------------------------------------------------------------------------
            // Place physical volumes for residues. 5 values for each physical volume
            // are recorded as a vector, specific to the current bp index, within the correct
            // map key:value pair of fpDnaMoleculePositions. Note the physical volumes themselves
            // not recorded in the map.
            // These values are: x, y, z, bp index (i.e. 1 - 200), nucleotide index (i.e. 1 or 2).
            // e.g. To get the x coordinate of the 150th sugar volume in the second DNA strand of
            // the 7th nucleosome: (*fpDnaMoleculePositions)["Desoxyribose"][(6*200*2)+(149*2)+2][0]
            //--------------------------------------------------------------------------------------
            G4int bp_index = (i*fNumBpPerNucleosome)+j;
            G4String bp_index_string = std::to_string(bp_index);

            // Phosphate 1
            //--------------------------------------------------------------------------------------
            G4String phys_name = "p_1_";
            G4VPhysicalVolume* sTMP1;
            if (fUseG4Volumes){
                sTMP1 = new G4PVPlacement(rotCuts,posSugarTMP1,volMap->at("sugarTMP1")[j],
                    phys_name,logicFiber,false,count);
            }
            else{
                sTMP1 = CreatePhysicalVolume(phys_name,count,true,volMap->at("sugarTMP1")[j],
                    rotCuts,&posSugarTMP1,logicFiber);
            }
            (*fpDnaMoleculePositions)["Phosphate"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP1.getX());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP1.getY());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP1.getZ());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(count);
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(1);

            // Sugar 1
            //--------------------------------------------------------------------------------------
            phys_name = "s_1_";
            G4VPhysicalVolume* sTHF1;
            if (fUseG4Volumes){
                sTHF1 = new G4PVPlacement(rotCuts,posSugarTHF1,volMap->at("sugarTHF1")[j],
                    phys_name,logicFiber,false,count+100000);
            }
            else{
                sTHF1 = CreatePhysicalVolume(phys_name,count+100000,true,volMap->at("sugarTHF1")[j],
                    rotCuts,&posSugarTHF1,logicFiber);
            }
            (*fpDnaMoleculePositions)["Desoxyribose"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF1.getX());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF1.getY());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF1.getZ());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(count);
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(1);

            // Base 1
            //--------------------------------------------------------------------------------------
            phys_name = "b_1_";
            G4VPhysicalVolume* base1;
            if (fUseG4Volumes){
                base1 = new G4PVPlacement(rotCuts,posBase1,volMap->at("base1")[j],phys_name,
                        logicFiber,false,count+200000);
            }
            else{
                base1 = CreatePhysicalVolume(phys_name,count+200000,true,volMap->at("base1")[j],
                    rotCuts,&posBase1,logicFiber);
            }
            (*fpDnaMoleculePositions)["Base1"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Base1"].back().push_back(posBase1.getX());
            (*fpDnaMoleculePositions)["Base1"].back().push_back(posBase1.getY());
            (*fpDnaMoleculePositions)["Base1"].back().push_back(posBase1.getZ());
            (*fpDnaMoleculePositions)["Base1"].back().push_back(count);
            (*fpDnaMoleculePositions)["Base1"].back().push_back(1);

            // Base 2
            //--------------------------------------------------------------------------------------
            phys_name = "b_2_";
            G4VPhysicalVolume* base2;
            if (fUseG4Volumes){
                base2 = new G4PVPlacement(rotCuts,posBase2,volMap->at("base2")[j],phys_name,
                                     logicFiber,false,count+1200000);
            }
            else{
                base2 = CreatePhysicalVolume(phys_name,count+1200000,true,volMap->at("base2")[j],
                    rotCuts,&posBase2,logicFiber);
            }
            (*fpDnaMoleculePositions)["Base2"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Base2"].back().push_back(posBase2.getX());
            (*fpDnaMoleculePositions)["Base2"].back().push_back(posBase2.getY());
            (*fpDnaMoleculePositions)["Base2"].back().push_back(posBase2.getZ());
            (*fpDnaMoleculePositions)["Base2"].back().push_back(count);
            (*fpDnaMoleculePositions)["Base2"].back().push_back(2);

            // Sugar 2
            //--------------------------------------------------------------------------------------
            phys_name = "s_2_";
            G4VPhysicalVolume* sTHF2;
            if (fUseG4Volumes){
                sTHF2 = new G4PVPlacement(rotCuts,posSugarTHF2,volMap->at("sugarTHF2")[j],
                    phys_name,logicFiber,false,count+1100000);
            }
            else{
                sTHF2 = CreatePhysicalVolume(phys_name,count+1100000,true,
                    volMap->at("sugarTHF2")[j],rotCuts,&posSugarTHF2,logicFiber);
            }
            (*fpDnaMoleculePositions)["Desoxyribose"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF2.getX());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF2.getY());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF2.getZ());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(count);
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(2);

            // Phosphate 2
            //--------------------------------------------------------------------------------------
            phys_name = "p_2_";
            G4VPhysicalVolume* sTMP2;
            if (fUseG4Volumes){
                sTMP2 = new G4PVPlacement(rotCuts,posSugarTMP2,volMap->at("sugarTMP2")[j],
                    phys_name,logicFiber,false,count+1000000);
            }
            else{
                sTMP2 = CreatePhysicalVolume(phys_name,count+1000000,true,
                    volMap->at("sugarTMP2")[j],rotCuts,&posSugarTMP2,logicFiber);
            }
            (*fpDnaMoleculePositions)["Phosphate"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP2.getX());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP2.getY());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP2.getZ());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(count);
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(2);

            //--------------------------------------------------------------------------------------
            // Place water volumes (containing residue placements) inside fiber volume
            //--------------------------------------------------------------------------------------
            // G4PVPlacement* sTMP1W = new G4PVPlacement(0,posSugarTMP1,
            //     volMap->at("sugarTMP1Water")[j],"sugarTMP1Hydra",logicFiber,false,count);
            // G4PVPlacement* sTHF1W = new G4PVPlacement(0,posSugarTHF1,
            //     volMap->at("sugarTHF1Water")[j],"sugarTHF1Hydra",logicFiber,false,count);
            // G4PVPlacement* b1W = new G4PVPlacement(0,posBase1,
            //     volMap->at("base1Water")[j],"base1Hydra",logicFiber,false,count);
            // G4PVPlacement* b2W = new G4PVPlacement(0,posBase2,
            //     volMap->at("base2Water")[j],"base2Hydra",logicFiber,false,count);
            // G4PVPlacement* sTHF2W = new G4PVPlacement(0,posSugarTHF2,
            //     volMap->at("sugarTHF2Water")[j],"sugarTHF2Hydra",logicFiber,false,count);
            // G4PVPlacement* sTMP2W = new G4PVPlacement(0,posSugarTMP2,
            //     volMap->at("sugarTMP2Water")[j],"sugarTMP2Hydra",logicFiber,false,count);

            //--------------------------------------------------------------------------------------
            // Check overlaps if solids have been cut.
            //--------------------------------------------------------------------------------------
            if (fCheckForOverlaps) {
                G4bool overlapDetected = false;
                if(sTMP1->CheckOverlaps(fOverlapsResolution) && fQuitIfOverlap)
                    ThrowOverlapError();
                if(sTHF1->CheckOverlaps(fOverlapsResolution) && fQuitIfOverlap)
                    ThrowOverlapError();
                if(base1->CheckOverlaps(fOverlapsResolution) && fQuitIfOverlap)
                    ThrowOverlapError();
                if(base2->CheckOverlaps(fOverlapsResolution) && fQuitIfOverlap)
                    ThrowOverlapError();
                if(sTHF2->CheckOverlaps(fOverlapsResolution) && fQuitIfOverlap)
                    ThrowOverlapError();
                if(sTMP2->CheckOverlaps(fOverlapsResolution) && fQuitIfOverlap)
                    ThrowOverlapError();
                // sTMP1W->CheckOverlaps();
                // sTHF1W->CheckOverlaps();
                // b1W->CheckOverlaps();
                // b2W->CheckOverlaps();
                // sTHF2W->CheckOverlaps();
                // sTMP2W->CheckOverlaps();
            }
            ++count;
        }

        //------------------------------------------------------------------------------------------
        // Place the histone volume
        //------------------------------------------------------------------------------------------
        G4ThreeVector posHistoneForNucleo = posHistone;
        // Rotate. Not really necessary
        posHistoneForNucleo.rotateZ(i*fFiberDeltaAngle);
        // Apply a z shift, specific to each nucleosome, to span length of fiber
        posHistoneForNucleo += G4ThreeVector(0.,0.,i*zShift);
        // Apply shift such that fiber helix construction begins at at one end.
        posHistoneForNucleo += minusForFiber;

        // Create volume
        G4String histName = "histone_" + std::to_string(i);
        G4VPhysicalVolume* pHistone;
        if (fUseG4Volumes){
            pHistone = new G4PVPlacement(0,posHistoneForNucleo,logicHistone,histName,logicFiber,true,i+2000000);
        }
        else{
            pHistone = CreatePhysicalVolume(histName,i+2000000,true,logicHistone,0,&posHistoneForNucleo,logicFiber);
        }
        (*fpDnaMoleculePositions)["Histone"].push_back(std::vector<double>());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(posHistoneForNucleo.getX());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(posHistoneForNucleo.getY());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(posHistoneForNucleo.getZ());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(i);
        (*fpDnaMoleculePositions)["Histone"].back().push_back(0);

        // Check for overlaps
        if (fCheckForOverlaps) {
            if(pHistone->CheckOverlaps(fOverlapsResolution) && fQuitIfOverlap)
                ThrowOverlapError();
        }
    }

    return logicFiber;
}


//--------------------------------------------------------------------------------------------------
// Create the solid and logical volumes required to build DNA around one histone.
// Return a map as:
// Key: name of the volume (base1, base2, base1Water, ...). Size = 12.
// Content: vector of corresponding logical volumes (each vector size = 200)
// Note: there is code to implement hydration shells around the residue volumes, but this is not
// currently implemented/tested.
//--------------------------------------------------------------------------------------------------
std::map<G4String, std::vector<G4LogicalVolume*> >* VoxelizedNuclearDNA::CreateNucleosomeCuttedSolidsAndLogicals(
    std::vector<DNAPlacementData>* nucleosomeVolumePositions, std::map<G4ThreeVector,
    G4double>* posAndRadiusMap)
{
    // This is the map to be returned
    std::map<G4String, std::vector<G4LogicalVolume*> >* logicSolidsMap = new std::map<G4String, std::vector<G4LogicalVolume*> >;

    G4int basePairNum = nucleosomeVolumePositions->size(); // 200

    //----------------------------------------------------------------------------------------------
    // Create elementary solids
    //----------------------------------------------------------------------------------------------
    // Throw error if a member variables hasn't been initialized correctly.
    if(fSugarTHFRadius==-1 || fSugarTMPRadius==-1 || fBaseRadius==-1)
    {
        G4cerr<<"************************************************************"<<G4endl;
        G4cerr<<"fSugarTHFRadius, fSugarTMPRadius or fBaseRadius were not set. Fatal error."<<G4endl;
        G4cerr<<"************************************************************"<<G4endl;
        std::exit(EXIT_FAILURE);
    }

    //----------------------------------------------------------------------------------------------
    // Create solid volumes
    //----------------------------------------------------------------------------------------------
    // residues
    G4Orb* solidSugarTHF = new G4Orb("solid_sugar_THF", fSugarTHFRadius);
    G4Orb* solidSugarTMP = new G4Orb("solid_sugar_TMP", fSugarTMPRadius);
    G4Orb* solidBase = new G4Orb("solid_base", fBaseRadius);

    // hydration shells
    // G4Orb* solidSugarTHFWater = new G4Orb("solid_sugar_THF_Water", fSugarTHFRadiusWater);
    // G4Orb* solidSugarTMPWater = new G4Orb("solid_sugar_TMP_Water", fSugarTMPRadiusWater);
    // G4Orb* solidBaseWater = new G4Orb("solid_base_Water", fBaseRadiusWater);

    //----------------------------------------------------------------------------------------------
    // Position variables for residues
    //----------------------------------------------------------------------------------------------
    G4ThreeVector posSugarTMP1;
    G4ThreeVector posSugarTHF1;
    G4ThreeVector posBase1;
    G4ThreeVector posBase2;
    G4ThreeVector posSugarTHF2;
    G4ThreeVector posSugarTMP2;

    //----------------------------------------------------------------------------------------------
    // Iterate over each base pair to generate cut solids and logical volumes.
    //----------------------------------------------------------------------------------------------
    for(int j=0;j<fNumBpPerNucleosome;++j)
    {
        //------------------------------------------------------------------------------------------
        // First: cut the solids (if requested)
        //------------------------------------------------------------------------------------------
        // Get the position
        posSugarTMP1 = nucleosomeVolumePositions->at(j).posSugarTMP1;
        posSugarTHF1 = nucleosomeVolumePositions->at(j).posSugarTHF1;
        posBase1 = nucleosomeVolumePositions->at(j).posBase1;
        posBase2 = nucleosomeVolumePositions->at(j).posBase2;
        posSugarTHF2 = nucleosomeVolumePositions->at(j).posSugarTHF2;
        posSugarTMP2 = nucleosomeVolumePositions->at(j).posSugarTMP2;

        // Variables for the cut solid volumes
        // residues
        G4VSolid* sugarTMP1;
        G4VSolid* sugarTHF1;
        G4VSolid* base1;
        G4VSolid* base2;
        G4VSolid* sugarTHF2;
        G4VSolid* sugarTMP2;

        // hydration shells
        // G4VSolid* sugarTMP1Water;
        // G4VSolid* sugarTHF1Water;
        // G4VSolid* base1Water;
        // G4VSolid* base2Water;
        // G4VSolid* sugarTHF2Water;
        // G4VSolid* sugarTMP2Water;


        // if fCutVolumes is true (i.e. need to run simulations), cut the volumes
        if(fCutVolumes)
        {
            // residues
            sugarTMP1 = CreateCutSolid(solidSugarTMP,posSugarTMP1,posAndRadiusMap, "sugarTMP");
            sugarTHF1 = CreateCutSolid(solidSugarTHF,posSugarTHF1,posAndRadiusMap, "sugarTHF");
            base1 = CreateCutSolid(solidBase,posBase1,posAndRadiusMap, "base");
            base2 = CreateCutSolid(solidBase,posBase2,posAndRadiusMap, "base");
            sugarTHF2 = CreateCutSolid(solidSugarTHF,posSugarTHF2,posAndRadiusMap, "sugarTHF");
            sugarTMP2 = CreateCutSolid(solidSugarTMP,posSugarTMP2,posAndRadiusMap, "sugarTMP");

            // hydration shells
            // sugarTMP1Water = CreateCutSolid(solidSugarTMPWater,posSugarTMP1,posAndRadiusMap);
            // sugarTHF1Water = CreateCutSolid(solidSugarTHFWater,posSugarTHF1,posAndRadiusMap);
            // base1Water = CreateCutSolid(solidBaseWater,posBase1,posAndRadiusMap);
            // base2Water = CreateCutSolid(solidBaseWater,posBase2,posAndRadiusMap);
            // sugarTHF2Water = CreateCutSolid(solidSugarTHFWater,posSugarTHF2,posAndRadiusMap);
            // sugarTMP2Water = CreateCutSolid(solidSugarTMPWater,posSugarTMP2,posAndRadiusMap);
        }
        // if fCutVolumes is false it means we just want to visualize the geometry so we do not need
        // the cutted volumes. Just use the uncut solids.
        else
        {
            //residues
            sugarTMP1 = solidSugarTMP;
            sugarTHF1 = solidSugarTHF;
            base1 = solidBase;
            base2 = solidBase;
            sugarTHF2 = solidSugarTHF;
            sugarTMP2 = solidSugarTMP;

            // hydration shells
            // sugarTMP1Water = solidSugarTMPWater;
            // sugarTHF1Water = solidSugarTHFWater;
            // base1Water = solidBaseWater;
            // base2Water = solidBaseWater;
            // sugarTHF2Water = solidSugarTHFWater;
            // sugarTMP2Water = solidSugarTMPWater;
        }

        //------------------------------------------------------------------------------------------
        // Then: create logical volumes using the cut solids
        //------------------------------------------------------------------------------------------
        // residues
        G4LogicalVolume* logicSugarTHF1;
        G4LogicalVolume* logicSugarTMP1;
        G4LogicalVolume* logicBase1;
        G4LogicalVolume* logicBase2;
        G4LogicalVolume* logicSugarTHF2;
        G4LogicalVolume* logicSugarTMP2;
        // hydration shells
        // G4LogicalVolume* logicSugarTMP1Water;
        // G4LogicalVolume* logicSugarTHF1Water;
        // G4LogicalVolume* logicBase1Water;
        // G4LogicalVolume* logicBase2Water;
        // G4LogicalVolume* logicSugarTHF2Water;
        // G4LogicalVolume* logicSugarTMP2Water;

        // Handle G4 vs Ts approach to generating logical volumes
        if (fUseG4Volumes) {
            logicSugarTMP1 = new G4LogicalVolume(sugarTMP1,fDNAMaterial,"Phosphate1");
            logicSugarTHF1 = new G4LogicalVolume(sugarTHF1,fDNAMaterial,"Sugar1");
            logicBase1 = new G4LogicalVolume(base1,fDNAMaterial,"Base1"); // PY
            logicBase2 = new G4LogicalVolume(base2,fDNAMaterial,"Base2"); // PU
            logicSugarTHF2 = new G4LogicalVolume(sugarTHF2,fDNAMaterial,"Sugar2");
            logicSugarTMP2 = new G4LogicalVolume(sugarTMP2,fDNAMaterial,"Phosphate2");
        }
        else {
            // Logical volumes for each nucleotide have the same name
            logicSugarTMP1 = CreateLogicalVolume("Phosphate1",fDNAMaterialName,sugarTMP1);
            logicSugarTHF1 = CreateLogicalVolume("Sugar1",fDNAMaterialName,sugarTHF1);
            logicBase1 = CreateLogicalVolume("Base1",fDNAMaterialName,base1);
            logicBase2 = CreateLogicalVolume("Base2",fDNAMaterialName,base2);
            logicSugarTHF2 = CreateLogicalVolume("Sugar2",fDNAMaterialName,sugarTHF2);
            logicSugarTMP2 = CreateLogicalVolume("Phosphate2",fDNAMaterialName,sugarTMP2);
        }

        // Creation of hydration shells
        // logicSugarTMP1Water = new G4LogicalVolume(sugarTMP1Water,fDNAMaterial,"logic_sugarTMP_1_hydra");
        // logicSugarTHF1Water = new G4LogicalVolume(sugarTHF1Water,fDNAMaterial,"logic_sugarTHF_1_hydra");
        // logicBase1Water = new G4LogicalVolume(base1Water, fDNAMaterial,"Base1_hydra");
        // logicBase2Water = new G4LogicalVolume(base2Water, fDNAMaterial,"Base2_hydra");
        // logicSugarTHF2Water = new G4LogicalVolume(sugarTHF2Water,fDNAMaterial,"logic_sugarTHF_2_hydra");
        // logicSugarTMP2Water = new G4LogicalVolume(sugarTMP2Water,fDNAMaterial,"logic_sugarTMP_2_hydra");

        //------------------------------------------------------------------------------------------
        // Save the logical volumes in the output map
        //------------------------------------------------------------------------------------------
        (*logicSolidsMap)["sugarTMP1"].push_back(logicSugarTMP1);
        (*logicSolidsMap)["sugarTHF1"].push_back(logicSugarTHF1);
        (*logicSolidsMap)["base1"].push_back(logicBase1);
        (*logicSolidsMap)["base2"].push_back(logicBase2);
        (*logicSolidsMap)["sugarTHF2"].push_back(logicSugarTHF2);
        (*logicSolidsMap)["sugarTMP2"].push_back(logicSugarTMP2);

        // (*logicSolidsMap)["sugarTMP1Water"].push_back(logicSugarTMP1Water);
        // (*logicSolidsMap)["sugarTHF1Water"].push_back(logicSugarTHF1Water);
        // (*logicSolidsMap)["base1Water"].push_back(logicBase1Water);
        // (*logicSolidsMap)["base2Water"].push_back(logicBase2Water);
        // (*logicSolidsMap)["sugarTHF2Water"].push_back(logicSugarTHF2Water);
        // (*logicSolidsMap)["sugarTMP2Water"].push_back(logicSugarTMP2Water);
    } // complete iterating over all bp in single nucleotide

    // Note: each vector of the logicSolidsMap has 200 elements
    return logicSolidsMap;
}


//--------------------------------------------------------------------------------------------------
// Algorithm for cutting DNA residue solids to avoid overlaps.
// Idea: we must have a reference and a target. The reference is the solid we are considering
// (described by parameters solidOrbRef & posRef) and that could be cut if an overlap is
// detected with the target solid. In a geometry, it implies we have to go through all the target
// solids for a given reference solid. Target solid info (position and radius) is included in tarMap
// This method will return the cut spherical reference solid.
//--------------------------------------------------------------------------------------------------
G4VSolid* VoxelizedNuclearDNA::CreateCutSolid(G4Orb *solidOrbRef,
                                       G4ThreeVector& posRef,
                                       std::map<G4ThreeVector,G4double>* tarMap,
                                       G4String volName)
{
    G4SubtractionSolid* solidCut(NULL); // container for the cut solid

    bool isCutted = false; // flag to indicate if volume has been cut yet
    bool isOurVol = false; // flag to indicate if target volume is our reference volume

    //----------------------------------------------------------------------------------------------
    // Define radius of the reference solid we are focusing on
    //----------------------------------------------------------------------------------------------
    G4double radiusRef (0);

    if(volName=="base") radiusRef = fBaseRadius;
    else if(volName=="sugarTHF") radiusRef = fSugarTHFRadius;
    else if(volName=="sugarTMP") radiusRef = fSugarTMPRadius;

    // Hydration shell is handled if no volName provided
    else radiusRef = solidOrbRef->GetRadius();

    //----------------------------------------------------------------------------------------------
    // iterate on all the residue volumes in the map (3600 elements), i.e. "targets"
    //----------------------------------------------------------------------------------------------
    std::map<G4ThreeVector,G4double>::iterator it;
    std::map<G4ThreeVector,G4double>::iterator ite;
    G4int count = 0;

    for(it=tarMap->begin(), ite=tarMap->end();it!=ite;++it)
    {
        G4ThreeVector posTar = it->first; // position of target
        G4double radiusTar = it->second; // radius of target
        G4double distance = std::abs( (posRef-posTar).getR() ); // 3D distance between ref & target

        //------------------------------------------------------------------------------------------
        // Check if target volume = reference volume (can only happen once)
        //------------------------------------------------------------------------------------------
        if(distance == 0 && !isOurVol)
        {
            isOurVol = true;
        }
        //------------------------------------------------------------------------------------------
        // Throw error if position of target and reference match more than once. Implies there are
        // two overlapping volumes at the same position.
        //------------------------------------------------------------------------------------------
        else if(distance == 0 && isOurVol)
        {
            G4cerr<<"DetectorConstruction::CreateCutSolid, Fatal Error. Two volumes are placed at the same position."<<G4endl;
            exit(EXIT_FAILURE);
        }
        //------------------------------------------------------------------------------------------
        // If the reference and target are overlapping, cut the reference. However, the target will
        // also be cut on another iteration of this loop. The cuts are performed at the "middle"
        // point of intersection between them. The mathematics proceed according to that for
        // sphere-sphere intersections.
        //------------------------------------------------------------------------------------------
        else if(distance <= radiusRef+radiusTar)
        {
            // Solid volume used to cut. Make size "equal" to that of larger between the reference
            // and current target.
            G4double sliceBoxSize = std::max(radiusRef,radiusTar);
            G4Box* sliceBox = new G4Box("solid box for cut", sliceBoxSize, sliceBoxSize, sliceBoxSize);

            //--------------------------------------------------------------------------------------
            // To calculate the position of the intersection center
            //--------------------------------------------------------------------------------------
            // Displacement vector between target and reference
            G4ThreeVector displacement_vector = posTar - posRef;
            // Find the middle overlap point between the target and reference
            G4double intersection = (pow(radiusRef,2)-pow(radiusTar,2)+pow(distance,2) ) / (2*distance) + sliceBox->GetZHalfLength();
            // Add small safety buffer
            intersection -= 0.001*nm;
            // Create a vector to the intersection position, where one edge of the slicing volume
            // will be placed
            G4ThreeVector posSlice = intersection *( displacement_vector/displacement_vector.getR() );

            //--------------------------------------------------------------------------------------
            // Calculate the necessary rotations.
            //--------------------------------------------------------------------------------------
            G4double phi = std::acos(posSlice.getZ()/posSlice.getR());
            G4double theta = std::acos( posSlice.getX() / ( posSlice.getR()*std::cos(M_PI/2.-phi) ) );

            if(posSlice.getY()<0) theta = -theta;

            G4ThreeVector rotAxisForPhi(1*nm,0.,0.);
            rotAxisForPhi.rotateZ(theta+M_PI/2);
            G4RotationMatrix *rotMat = new G4RotationMatrix;
            rotMat->rotate(-phi, rotAxisForPhi);

            G4ThreeVector rotZAxis(0.,0.,1*nm);
            rotMat->rotate(theta, rotZAxis);

            //--------------------------------------------------------------------------------------
            // Create the G4SubtractionSolid.
            //--------------------------------------------------------------------------------------
            if(!isCutted) solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, sliceBox, rotMat, posSlice);
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, sliceBox, rotMat, posSlice);

            isCutted = true;
        }
        count++;
    }

    if(isCutted) return solidCut;
    else return solidOrbRef;
}


//--------------------------------------------------------------------------------------------------
// This method arranges identical DNA fibers in a cubic voxel, and returns the logical volume of
// that voxel. Current implementation places 20 fibers in fractal pattern, similar to that
// described by Zhu et al. (2020). DOI: 10.1667/RR15531.1
//--------------------------------------------------------------------------------------------------
G4LogicalVolume* VoxelizedNuclearDNA::ConstructLogicalVoxel(G4LogicalVolume* lFiber) {
        // Make empty voxel
        G4Box* voxelSolid = new G4Box("solid_voxel",fVoxelSideLength,fVoxelSideLength,fVoxelSideLength);
        G4LogicalVolume* voxelLogical = CreateLogicalVolume("Voxel",fWaterName,voxelSolid);

        // Populate voxel with fibres. Arrange in fractal pattern
        // Fractal loop 1 (7 fibers)
        //------------------------------------------------------------------------------------------
        G4double posXf1 = 119*nm;
        G4double posYf1 = 119*nm;
        G4double posZf1 = 34*nm;
        G4ThreeVector fibrePlacement = G4ThreeVector(posXf1,posYf1,posZf1);
        G4RotationMatrix* fibreRotation = new G4RotationMatrix();
        CreatePhysicalVolume("Fiber",0,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        G4double posXf2 = 34*nm;
        G4double posYf2 = 119*nm;
        G4double posZf2 = -51*nm;
        fibrePlacement = G4ThreeVector(posXf2,posYf2,posZf2);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(3*CLHEP::pi/2);
        CreatePhysicalVolume("Fiber",1,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        G4double posXf3 = -51*nm;
        G4double posYf3 = 119*nm;
        G4double posZf3 = 34*nm;
        fibrePlacement = G4ThreeVector(posXf3,posYf3,posZf3);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(CLHEP::pi);
        CreatePhysicalVolume("Fiber",2,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        G4double posXf4 = -51*nm;
        G4double posYf4 = 34*nm;
        G4double posZf4 = 119*nm;
        fibrePlacement = G4ThreeVector(posXf4,posYf4,posZf4);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(CLHEP::pi);
        fibreRotation->rotateX(CLHEP::pi/2);
        CreatePhysicalVolume("Fiber",3,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        G4double posXf5 = -51*nm;
        G4double posYf5 = -51*nm;
        G4double posZf5 = 34*nm;
        fibrePlacement = G4ThreeVector(posXf5,posYf5,posZf5);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(CLHEP::pi);
        fibreRotation->rotateX(CLHEP::pi);
        CreatePhysicalVolume("Fiber",4,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        G4double posXf6 = 34*nm;
        G4double posYf6 = -51*nm;
        G4double posZf6 = -51*nm;
        fibrePlacement = G4ThreeVector(posXf6,posYf6,posZf6);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(3*CLHEP::pi/2);
        fibreRotation->rotateX(CLHEP::pi);
        CreatePhysicalVolume("Fiber",5,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        G4double posXf7 = 119*nm;
        G4double posYf7 = -51*nm;
        G4double posZf7 = 34*nm;
        fibrePlacement = G4ThreeVector(posXf7,posYf7,posZf7);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateX(CLHEP::pi);
        CreatePhysicalVolume("Fiber",6,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);


        // Fractal loop 2 (7 fibers)
        //------------------------------------------------------------------------------------------
        G4double loopshift = 34*nm;
        fibrePlacement = G4ThreeVector(posXf1-loopshift,posYf1-loopshift,posZf1-loopshift);
        fibreRotation = new G4RotationMatrix();
        CreatePhysicalVolume("Fiber",7,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        fibrePlacement = G4ThreeVector(posXf2-loopshift,posYf2-loopshift,posZf2-loopshift);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(3*CLHEP::pi/2);
        CreatePhysicalVolume("Fiber",8,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        fibrePlacement = G4ThreeVector(posXf3-loopshift,posYf3-loopshift,posZf3-loopshift);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(CLHEP::pi);
        CreatePhysicalVolume("Fiber",9,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        fibrePlacement = G4ThreeVector(posXf4-loopshift,posYf4-loopshift,posZf4-loopshift);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(CLHEP::pi);
        fibreRotation->rotateX(CLHEP::pi/2);
        CreatePhysicalVolume("Fiber",10,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        fibrePlacement = G4ThreeVector(posXf5-loopshift,posYf5-loopshift,posZf5-loopshift);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(CLHEP::pi);
        fibreRotation->rotateX(CLHEP::pi);
        CreatePhysicalVolume("Fiber",11,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        fibrePlacement = G4ThreeVector(posXf6-loopshift,posYf6-loopshift,posZf6-loopshift);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(3*CLHEP::pi/2);
        fibreRotation->rotateX(CLHEP::pi);
        CreatePhysicalVolume("Fiber",12,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        fibrePlacement = G4ThreeVector(posXf7-loopshift,posYf7-loopshift,posZf7-loopshift);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateX(CLHEP::pi);
        CreatePhysicalVolume("Fiber",13,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        // fractal loop 3 (6 fibers)
        //------------------------------------------------------------------------------------------
        loopshift *= 2;
        fibrePlacement = G4ThreeVector(posXf1-loopshift,posYf1-loopshift,posZf1-loopshift);
        fibreRotation = new G4RotationMatrix();
        CreatePhysicalVolume("Fiber",14,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        fibrePlacement = G4ThreeVector(posXf2-loopshift,posYf2-loopshift,posZf2-loopshift);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(3*CLHEP::pi/2);
        CreatePhysicalVolume("Fiber",15,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        fibrePlacement = G4ThreeVector(posXf3-loopshift,posYf3-loopshift,posZf3-loopshift);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(CLHEP::pi);
        CreatePhysicalVolume("Fiber",16,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        fibrePlacement = G4ThreeVector(posXf4-loopshift,posYf4-loopshift,posZf4-loopshift);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(CLHEP::pi);
        fibreRotation->rotateX(CLHEP::pi/2);
        CreatePhysicalVolume("Fiber",17,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        fibrePlacement = G4ThreeVector(posXf5-loopshift,posYf5-loopshift,posZf5-loopshift);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(CLHEP::pi);
        fibreRotation->rotateX(CLHEP::pi);
        CreatePhysicalVolume("Fiber",18,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        fibrePlacement = G4ThreeVector(posXf6-loopshift,posYf6-loopshift,posZf6-loopshift);
        fibreRotation = new G4RotationMatrix();
        fibreRotation->rotateY(3*CLHEP::pi/2);
        fibreRotation->rotateX(CLHEP::pi);
        CreatePhysicalVolume("Fiber",19,true,lFiber,fibreRotation,&fibrePlacement,voxelLogical);

        return voxelLogical;
}


//--------------------------------------------------------------------------------------------------
// Helper function used if a geometrical overlap is detected. Return the standard Topas message
// about geometry overlaps and exit gracefully.
//--------------------------------------------------------------------------------------------------
void VoxelizedNuclearDNA::ThrowOverlapError()
{
    G4cout << "Topas is quitting due to the above geometry overlap problem." << G4endl;
    G4cout << "If you still want the TOPAS session to continue" << G4endl;
    G4cout << "(such as to use visualization to study the overlap)," << G4endl;
    G4cout << "Set the parameter Ge/QuitIfOverlapDetected to False" << G4endl;
    exit(0);
}
