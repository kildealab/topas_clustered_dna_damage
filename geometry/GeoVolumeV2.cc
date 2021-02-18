// Extra Class for TsDNAFabric
//
//**************************************************************************************************
// This class handles the construction of a logical volume for a chromatin fiber using the method
// BuildLogicFiber(). Additional methods provide supporting functionality for this method (e.g. the
// generation of physical volumes and logical volumes for residues and performing the necesary
// cutting procedures to prevent geometrical overlaps).
//**************************************************************************************************

#include "GeoVolumeV2.hh"

#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Exception.hh"

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
// Constructor. Paramters used to set (i) verbosity of console output and (ii) a scaling factor on
// the size of all geometry components (should be 1 for most purposes).
// Member variables are not initialized here, rather they must be initialized explicitly via setter 
// methods. In practice, this is done in the initialize method of GeoManager, using information from
// a pre-existing GeoCalculation object. Arbitrary negative values are assigned here, and used as
// error checking in BuildLogicFiber().
//--------------------------------------------------------------------------------------------------
GeoVolumeV2::GeoVolumeV2(G4int verbose, G4double factor) :
    fVerbose(verbose), fFactor(factor)
{
    // Water is defined from NIST material database
    fWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

    // G4Material* H20_mod = G4NistManager::Instance()->BuildMaterialWithNewDensity("G4_WATER_MODIFIEDg","G4_WATER",1.100*g/cm3);

    // All the variables are set to -1 for checking in Build methods
    fSugarTHFRadiusWater = -1*fFactor*m;
    fSugarTMPRadiusWater = -1*fFactor*m;
    fSugarTHFRadius = -1*fFactor*m;
    fSugarTMPRadius = -1*fFactor*m;
    fBaseRadiusWater = -1*fFactor*m;
    fBaseRadius = -1*fFactor*m;
    fFiberPitch = -1*fFactor*m;
    fFiberNbNuclPerTurn = -1*fFactor*m;
    fNucleoNum = -1;
    fBpNum = -1;
    fHistoneHeight = -1*fFactor*m;
    fHistoneRadius = -1*fFactor*m;

    fpDnaMoleculePositions = new std::map<G4String, std::vector<std::vector<G4double> > >();

    //DefineMaterial();
}

//--------------------------------------------------------------------------------------------------
// Destructor. Delete variables allocated with new.
//--------------------------------------------------------------------------------------------------
GeoVolumeV2::~GeoVolumeV2()
{
    delete fpDnaMoleculePositions;
}

//--------------------------------------------------------------------------------------------------
// Create and return a logical volume for a chromatin fiber. 
// Within this logical volume are the physical volumes for the histones, the resiudes and their
// hydration shells. Solids and logicals are generated for the histones within this metohd directly,
// whereas those for the residues are generated using CreateNucleosomeCuttedSolidsAndLogicals().
// A map, fpDnaMoleculePositions, containing the coordinates for all residues and histones is filled
// and can be accessed using GetDNAMoleculesPositions().
//--------------------------------------------------------------------------------------------------
G4LogicalVolume* GeoVolumeV2::BuildLogicFiber(std::vector<std::vector<DNAPlacementData> >* dnaVolPos,
                                            std::vector<G4ThreeVector>* posNucleo,
                                            std::map<G4ThreeVector, G4double>* posAndRadiusMap,
                                            G4Material* DNAMaterial,
                                            G4bool cutVolumes,
                                            G4bool checkForOverlaps,
                                            G4int overlapsResolution,
                                            G4bool quitIfOverlap)
{
    //----------------------------------------------------------------------------------------------
    // Throw error if any of the member variables haven't been initialized correctly.
    //----------------------------------------------------------------------------------------------
    if(fSugarTHFRadius==-1||fSugarTMPRadius==-1||fBaseRadius==-1||fFiberPitch==-1||fFiberNbNuclPerTurn==-1
            ||fNucleoNum==-1||fBpNum==-1||fHistoneHeight==-1||fHistoneRadius==-1)
    {
        G4cerr<<"FatalError: GeoVolumeV2::BuildLogicFiber. A class parameter has not been " 
            << "initialized and its value is still negative"<<G4endl;
        std::exit(EXIT_FAILURE);
    }

    //----------------------------------------------------------------------------------------------
    // Create cylindrical fiber volume
    //----------------------------------------------------------------------------------------------
    G4Tubs* solidFiber = new G4Tubs("solid histone", 0., 17.*fFactor*nm, 68.*fFactor*nm, 0, 360); // radius & half length

    G4VisAttributes fiberVis(G4Colour(1.0, 1.0, 1.0, 0.1) );
    // fiberVis.SetVisibility(false);
    fiberVis.SetForceSolid(true);
    // fiberVis.SetForceWireframe(true);
    G4LogicalVolume* logicFiber = new G4LogicalVolume(solidFiber, fWater,"logic_fiber");
    logicFiber->SetVisAttributes(fiberVis);

    //----------------------------------------------------------------------------------------------
    // Create the histone volume
    //----------------------------------------------------------------------------------------------
    G4Tubs* solidHistone = new G4Tubs("solid histone", 0., fHistoneRadius, fHistoneHeight, 0, 360);

    G4VisAttributes histoneVis(G4Colour(1.0, 1.0, 1.0) );
    // histoneVis.SetVisibility(false);
    histoneVis.SetForceSolid(true);
    G4LogicalVolume* logicHistone = new G4LogicalVolume(solidHistone,fWater,"logic_histone");
    logicHistone->SetVisAttributes(histoneVis);

    //----------------------------------------------------------------------------------------------
    // Once the material for the fibre & nucleosomes has been defined (above), change the material
    // for the DNA elements (residues).
    //----------------------------------------------------------------------------------------------
    fWater = DNAMaterial;

    //----------------------------------------------------------------------------------------------
    // Generate the cut solids
    //----------------------------------------------------------------------------------------------
    // For the positions, we only use the second (number 1) nucleosome.
    // Indeed, it is a "middle" nucleosome and, thus, the two extemities will be cutted.
    // If we took the first, we would have a first volume non cutted and, thus, some overlaps.
    // Note dnaVolPos is paramater passed to BuildLogicFiber() that is essentially the output of 
    // GeoCalculation::CalculateDNAPosition(). I.e. all nucleotide positions around 3 basis
    // histone complexes.
    std::vector<DNAPlacementData>* nuclVolPos = &dnaVolPos->at(1);

    // We create here all the DNA volumes (solid & logical) around the histone based on the second 
    // nucleosome (number 1) positions. Only the volumes of the second nucleosome are generated. By 
    // placing them several times we will build the fiber. This is done to save memory and improve
    // speed. We saved the volumes in a map (key = name of the volume [e.g. sugar1], value = vector
    // of corresponding logical volumes).
    // Note posAndRadiusMap is parameter passed to BuildLogicFiber() that is essentially the output
    // of GeoCalculation::GenerateCoordAndRadiusMap(). I.e. a map of water radii for 6 volumes in
    // each of 200 bp in each of 3 basis nucleosomes (3600 volumes)
    std::map<G4String, std::vector<G4LogicalVolume*> >* volMap
            = CreateNucleosomeCuttedSolidsAndLogicals(nuclVolPos, posAndRadiusMap, cutVolumes);
    // Resulting volMap is indexed by one of 12 entries (6 residues & 6 hydration shells). Each
    // entry has 200 elements, each corresponding to a distinct logical volume

    // Print the mean volume number of the DNA volumes
    if(fVerbose>2)
    {
        CalculateMeanVol(volMap);
    }

    //----------------------------------------------------------------------------------------------
    // Save the positions of all residue volumes in the first nucleosome (not all 3) in vectors
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

    // Note: fBpNum is input in parameter file. Default = 200
    // Get the bp volume positions of the middle basis nucleosome (which is used to generate cut 
    // solids) and save in a vector to be accessed & rotated in the following forloop
    for(int j=0;j<fBpNum;++j)
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
    // Do the placements; i.e. build the nucleosome helix inside the fiber using G4PVPlacements of
    // the already-determined logical volumes
    //----------------------------------------------------------------------------------------------
    // Distance value used to place the volumes at the beginning of the fiber
    G4ThreeVector minusForFiber = G4ThreeVector(0.,0.,-solidFiber->GetDz() + fHistoneHeight);

    G4double zShift = fFiberPitch/fFiberNbNuclPerTurn; // z shift for each new nuclosome
    G4int count = 0;

    //----------------------------------------------------------------------------------------------
    // iterate on each nucleosome (fNucleoNum is input in parameter file. Default = 90).
    //----------------------------------------------------------------------------------------------
    for(int i=0;i<fNucleoNum;++i)
    {
        // Rotates the nucleosome itself with respect to the z-axis, in order to align the cut
        // volumes appropriately to prevent overlaps. This will be appliedThis is not the same as placing a nucleosome
        // at the next position around the fiber. That is done below. Our basis nucleosome is
        // actually the second one (middle). The following rotation logic accounts for this.

        // This rotation object will be applied to every physical volume placement below, in order
        // to align the cut bp volumes and prevent overlaps. This is a rotation about the
        // volume's own z-axis.
        G4RotationMatrix* rotCuts = new G4RotationMatrix();
        rotCuts->rotateZ((i-1)*-fFiberDeltaAngle);

        //------------------------------------------------------------------------------------------
        // iterate on each bp (fBpNum is input in parameter file. Default = 200)
        //------------------------------------------------------------------------------------------
        for(int j=0;j<fBpNum;++j)
        {
            //--------------------------------------------------------------------------------------
            // The following positional rotations account for the fact that nucleosomes are
            // assembled in a helical fashion around the fibre axis (z). Note this has nothing to
            // do with rotations of the bp volumes themselves (which are handled by rotCuts obj).
            //--------------------------------------------------------------------------------------
            // Since our basis nucleosome is the middle one, when i==1, there should be no rotation.
            // Before and after this, there are corresponding rotations about the z axis.
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
            // Place physical volumes. Residue volumes are placed inside of corresponding water
            // volumes (i.e. water volumes are mother volumes). 5 values for each physical volume
            // are then recorded as a vector, specific to the current bp index, within the correct
            // map key:value pair of fpDnaMoleculePositions. Note the physical volumes themselves
            // not recorded in the map.  
            // These values are: x, y, z, bp index (i.e. 1 - 200), nucleotide index (i.e. 1 or 2).
            // Subsequently (below), the water volumes are placed inside the cylindrical logicFiber 
            // volume. 
            // e.g. To get the x coordinate of the 150th sugar volume in the second DNA strand of 
            // the 7th nucleosome: (*fpDnaMoleculePositions)["Desoxyribose"][(6*200*2)+(149*2)+2][0]
            //--------------------------------------------------------------------------------------
            G4int bp_index = (i*fBpNum)+j;
            G4String bp_index_string = std::to_string(bp_index);
            // Phosphate 1
            G4String phys_name = "p_1_" + bp_index_string;
            G4PVPlacement* sTMP1 = new G4PVPlacement(rotCuts,posSugarTMP1,volMap->at("sugarTMP1")[j],
                phys_name,logicFiber,false,count);
            (*fpDnaMoleculePositions)["Phosphate"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP1.getX());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP1.getY());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP1.getZ());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(count);
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(1);
            // Sugar 1
            phys_name = "s_1_" + bp_index_string;
            G4PVPlacement* sTHF1 = new G4PVPlacement(rotCuts,posSugarTHF1,volMap->at("sugarTHF1")[j],
                phys_name,logicFiber,false,count+100000);
            (*fpDnaMoleculePositions)["Desoxyribose"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF1.getX());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF1.getY());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF1.getZ());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(count);
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(1);
            // Base 1
            phys_name = "b_1_" + bp_index_string;
            G4PVPlacement* base1 = new G4PVPlacement(rotCuts,posBase1,volMap->at("base1")[j],phys_name,
                    logicFiber,false,count+200000);
            (*fpDnaMoleculePositions)["Base1"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Base1"].back().push_back(posBase1.getX());
            (*fpDnaMoleculePositions)["Base1"].back().push_back(posBase1.getY());
            (*fpDnaMoleculePositions)["Base1"].back().push_back(posBase1.getZ());
            (*fpDnaMoleculePositions)["Base1"].back().push_back(count);
            (*fpDnaMoleculePositions)["Base1"].back().push_back(1);
            // Base 2
            phys_name = "b_2_" + bp_index_string;
            G4PVPlacement* base2 = new G4PVPlacement(rotCuts,posBase2,volMap->at("base2")[j],phys_name,
                                 logicFiber,false,count+1200000);
            (*fpDnaMoleculePositions)["Base2"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Base2"].back().push_back(posBase2.getX());
            (*fpDnaMoleculePositions)["Base2"].back().push_back(posBase2.getY());
            (*fpDnaMoleculePositions)["Base2"].back().push_back(posBase2.getZ());
            (*fpDnaMoleculePositions)["Base2"].back().push_back(count);
            (*fpDnaMoleculePositions)["Base2"].back().push_back(2);
            // Sugar 2
            phys_name = "s_2_" + bp_index_string;
            G4PVPlacement* sTHF2 = new G4PVPlacement(rotCuts,posSugarTHF2,volMap->at("sugarTHF2")[j],
                phys_name,logicFiber,false,count+1100000);
            (*fpDnaMoleculePositions)["Desoxyribose"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF2.getX());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF2.getY());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF2.getZ());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(count);
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(2);
            // Phosphate 2
            phys_name = "p_2_" + bp_index_string;
            G4PVPlacement* sTMP2 = new G4PVPlacement(rotCuts,posSugarTMP2,volMap->at("sugarTMP2")[j],
                phys_name,logicFiber,false,count+1000000);
            (*fpDnaMoleculePositions)["Phosphate"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP2.getX());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP2.getY());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP2.getZ());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(count);
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(2);
            // Place water volumes (containing residue placements) inside fiber volume
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

            // Check overlaps if solids have been cut. No need if uncut, b/c will overlap
            // Throw and overlap error if settings request as such
            if (checkForOverlaps) {
                G4bool overlapDetected = false;
                if(sTMP1->CheckOverlaps(overlapsResolution) && quitIfOverlap)
                    ThrowOverlapError();
                if(sTHF1->CheckOverlaps(overlapsResolution) && quitIfOverlap)
                    ThrowOverlapError();
                if(base1->CheckOverlaps(overlapsResolution) && quitIfOverlap)
                    ThrowOverlapError();
                if(base2->CheckOverlaps(overlapsResolution) && quitIfOverlap)
                    ThrowOverlapError();
                if(sTHF2->CheckOverlaps(overlapsResolution) && quitIfOverlap)
                    ThrowOverlapError();
                if(sTMP2->CheckOverlaps(overlapsResolution) && quitIfOverlap)
                    ThrowOverlapError();
                // sTMP1W->CheckOverlaps();
                // sTHF1W->CheckOverlaps();
                // b1W->CheckOverlaps();
                // b2W->CheckOverlaps();
                // sTHF2W->CheckOverlaps();
                // sTMP2W->CheckOverlaps();


                // 
                // Simulation results can not be trusted when any overlap exists.
                // If you still want the TOPAS session to continue
                // (such as to use visualization to study the overlap),
                // Set the parameter Ge/QuitIfOverlapDetected to False
            }

            ++count;
        }

        //------------------------------------------------------------------------------------------
        // Place the histone volume
        //------------------------------------------------------------------------------------------
        // Rotate. Not really necessary
        G4ThreeVector posHistoneForNucleo = posHistone;
        posHistoneForNucleo.rotateZ(i*fFiberDeltaAngle);
        // Apply a z shift, specific to each nucleosome, to span length of fiber
        posHistoneForNucleo += G4ThreeVector(0.,0.,i*zShift);
        // Apply shift such that fiber helix construction begins at at one end.
        posHistoneForNucleo += minusForFiber;
        // Place
        G4String histName = "histone_" + std::to_string(i);
        G4VPhysicalVolume* pHistone = new G4PVPlacement(0,posHistoneForNucleo,logicHistone,histName,logicFiber,true,i);
        (*fpDnaMoleculePositions)["Histone"].push_back(std::vector<double>());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(posHistoneForNucleo.getX());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(posHistoneForNucleo.getY());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(posHistoneForNucleo.getZ());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(i);
        (*fpDnaMoleculePositions)["Histone"].back().push_back(0);

        if (checkForOverlaps) {
            if(pHistone->CheckOverlaps(overlapsResolution) && quitIfOverlap)
                ThrowOverlapError();
        }
    }

    // logicFiber contains the placements of the histones and the water volumes, which in turn
    // contain the placements of the residues. 
    // Note: The map of positions (fpDnaMoleculePositions) is not returned explicitly but is
    // accessible via a getter method, GetDNAMoleculesPositions() GeoManager has a wrapper for this,
    // also called GetDNAMoleculesPositions()
    // CalculateMeanVol(volMap);
    return logicFiber;
}


//--------------------------------------------------------------------------------------------------
// Create the solid and logical volumes required to build DNA around one histone.
// Return a map as:
// Key: name of the volume (base1, base2, base1Water, ...). Size = 12.
// Content: vector of corresponding logical volumes (each vector size = 200)
//--------------------------------------------------------------------------------------------------
std::map<G4String, std::vector<G4LogicalVolume*> >* GeoVolumeV2::CreateNucleosomeCuttedSolidsAndLogicals(
    std::vector<DNAPlacementData>* nucleosomeVolumePositions, std::map<G4ThreeVector, 
    G4double>* posAndRadiusMap, G4bool cutVolumes)
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
        G4cerr<<"GeoVolumeV2::CreateNucleosomeCuttedSolidsAndLogicals: fSugarTHFRadius, "
            << "fSugarTMPRadius or fBaseRadius were not set. Fatal error."<<G4endl;
        G4cerr<<"************************************************************"<<G4endl;
        std::exit(EXIT_FAILURE);
    }

    //----------------------------------------------------------------------------------------------
    // Visibility attributes for nucleotide components (base, sugar, phosphate, shell)
    //----------------------------------------------------------------------------------------------
    // G4VisAttributes red(G4Colour(1.0, 0.0, 0.0) );
    // G4VisAttributes blue(G4Colour(0.0, 0.0, 1.0, 0.0) );
    // G4VisAttributes green(G4Colour(0.0, 1.0, 0.0) );
    // G4VisAttributes yellow(G4Colour(1.0, 1.0, 0.0) );
    // red.SetForceSolid(true);
    // // red.SetVisibility(false);
    // green.SetForceSolid(true);
    // yellow.SetForceSolid(true);

    G4VisAttributes visBase(G4Colour(0.92, 0.6, 0.6));    
    G4VisAttributes visSugar(G4Colour(0.43, 0.62, 0.92));    
    G4VisAttributes visPhosphate(G4Colour(0.71, 0.65, 0.84)); 
    G4VisAttributes visHydration(G4Colour(0.27, 0.82, 0.82)); 
    visBase.SetForceSolid(true);
    visSugar.SetForceSolid(true);
    visPhosphate.SetForceSolid(true);
    // visBase.SetForceWireframe(true);
    // visSugar.SetForceWireframe(true);
    // visPhosphate.SetForceWireframe(true);
    // visBase.SetVisibility(false);
    // visSugar.SetVisibility(false);
    // visPhosphate.SetVisibility(false);

    // visHydration.SetForceSolid(true);
    visHydration.SetForceWireframe(true);
    // visHydration.SetVisibility(true);
    visHydration.SetVisibility(false);

    //----------------------------------------------------------------------------------------------
    // Create solid volumes
    //----------------------------------------------------------------------------------------------
    // residues
    G4Orb* solidSugarTHF = new G4Orb("solid_sugar_THF", fSugarTHFRadius);
    G4Orb* solidSugarTMP = new G4Orb("solid_sugar_TMP", fSugarTMPRadius);
    G4Orb* solidBase = new G4Orb("solid_base", fBaseRadius);

    // hydration shells
    G4Orb* solidSugarTHFWater = new G4Orb("solid_sugar_THF_Water", fSugarTHFRadiusWater);
    G4Orb* solidSugarTMPWater = new G4Orb("solid_sugar_TMP_Water", fSugarTMPRadiusWater);
    G4Orb* solidBaseWater = new G4Orb("solid_base_Water", fBaseRadiusWater);

    //----------------------------------------------------------------------------------------------
    // Position variables of resiudes
    //---------------------------------------------------------------------------------------------- 
    G4ThreeVector posSugarTMP1;
    G4ThreeVector posSugarTHF1;
    G4ThreeVector posBase1;
    G4ThreeVector posBase2;
    G4ThreeVector posSugarTHF2;
    G4ThreeVector posSugarTMP2;

    //----------------------------------------------------------------------------------------------
    // Iterate on each base pair to generate cut solids and logical volumes.
    //----------------------------------------------------------------------------------------------
    for(int j=0;j<fBpNum;++j)
    {
        //------------------------------------------------------------------------------------------
        // First: cut the solids (if no visualization)
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


        // if cutVolumes is true it means we will run calculations so create the cut volumes
        // if(true)
        if(cutVolumes)
        {
            // residues
            sugarTMP1 = CreateCutSolid(solidSugarTMP,posSugarTMP1,posAndRadiusMap, "sugarTMP", true);
            sugarTHF1 = CreateCutSolid(solidSugarTHF,posSugarTHF1,posAndRadiusMap, "sugarTHF", true);
            base1 = CreateCutSolid(solidBase,posBase1,posAndRadiusMap, "base", true);
            base2 = CreateCutSolid(solidBase,posBase2,posAndRadiusMap, "base", true);
            sugarTHF2 = CreateCutSolid(solidSugarTHF,posSugarTHF2,posAndRadiusMap, "sugarTHF", true);
            sugarTMP2 = CreateCutSolid(solidSugarTMP,posSugarTMP2,posAndRadiusMap, "sugarTMP", true);

            // hydration shells
            // sugarTMP1Water = CreateCutSolid(solidSugarTMPWater,posSugarTMP1,posAndRadiusMap);
            // sugarTHF1Water = CreateCutSolid(solidSugarTHFWater,posSugarTHF1,posAndRadiusMap);
            // base1Water = CreateCutSolid(solidBaseWater,posBase1,posAndRadiusMap);
            // base2Water = CreateCutSolid(solidBaseWater,posBase2,posAndRadiusMap);
            // sugarTHF2Water = CreateCutSolid(solidSugarTHFWater,posSugarTHF2,posAndRadiusMap);
            // sugarTMP2Water = CreateCutSolid(solidSugarTMPWater,posSugarTMP2,posAndRadiusMap);
        }
        // if cutVolumes is false it means we just want to visualize the geometry so we do not need
        // the cutted volumes. Just assign the uncut solids to the empty variables.
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

        // G4double sugarVolume = sugarTMP1->GetCubicVolume();
        // G4cout << "**********************************************" << G4endl;
        // G4cout << "Sugar volume = " << sugarVolume/m3 << G4endl;
        // G4cout << "**********************************************" << G4endl;

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

        // Creation of residues
        logicSugarTMP1 = new G4LogicalVolume(sugarTMP1,fWater,"logic_sugar_TMP_1");
        logicSugarTHF1 = new G4LogicalVolume(sugarTHF1,fWater,"logic_sugar_THF_1");
        logicBase1 = new G4LogicalVolume(base1,fWater,"logic_base_1"); // PY
        logicBase2 = new G4LogicalVolume(base2,fWater,"logic_base_2"); // PU
        logicSugarTHF2 = new G4LogicalVolume(sugarTHF2,fWater,"logic_sugar_THF_2");
        logicSugarTMP2 = new G4LogicalVolume(sugarTMP2,fWater,"logic_sugar_TMP_2");

        // Creation of hydration shells
        // logicSugarTMP1Water = new G4LogicalVolume(sugarTMP1Water,fWater,"logic_sugarTMP_1_hydra");
        // logicSugarTHF1Water = new G4LogicalVolume(sugarTHF1Water,fWater,"logic_sugarTHF_1_hydra");
        // logicBase1Water = new G4LogicalVolume(base1Water, fWater,"logic_base_1_hydra");
        // logicBase2Water = new G4LogicalVolume(base2Water, fWater,"logic_base_2_hydra");
        // logicSugarTHF2Water = new G4LogicalVolume(sugarTHF2Water,fWater,"logic_sugarTHF_2_hydra");
        // logicSugarTMP2Water = new G4LogicalVolume(sugarTMP2Water,fWater,"logic_sugarTMP_2_hydra");

        //------------------------------------------------------------------------------------------
        // Set visualization attributes for residues & hydration shell
        //------------------------------------------------------------------------------------------
        logicSugarTMP1->SetVisAttributes(visPhosphate);
        logicSugarTHF1->SetVisAttributes(visSugar);
        logicBase1->SetVisAttributes(visBase);
        logicBase2->SetVisAttributes(visBase);
        logicSugarTHF2->SetVisAttributes(visSugar);
        logicSugarTMP2->SetVisAttributes(visPhosphate);

        // logicSugarTMP1Water->SetVisAttributes(visHydration);
        // logicSugarTHF1Water->SetVisAttributes(visHydration);
        // logicBase1Water->SetVisAttributes(visHydration);
        // logicBase2Water->SetVisAttributes(visHydration);
        // logicSugarTHF2Water->SetVisAttributes(visHydration);
        // logicSugarTMP2Water->SetVisAttributes(visHydration);

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
// Cut algorithm to avoid overlaps.
// Idea: we must have a reference and a target. The reference is the solid we are considering (
// described by parameters solidOrbRef & posRef) and which could be cut if an overlap is
// detected with the target solid. In a geometry, it implies we have to go through all the target 
// solids for a given reference solid. Target solid info (position and radius) is included in tarMap
// This method will return the cut spherical reference solid.
// Note: fifth parameter, i.e. G4bool in doesn't seem to be used.
//--------------------------------------------------------------------------------------------------
G4VSolid* GeoVolumeV2::CreateCutSolid(G4Orb *solidOrbRef,
                                               G4ThreeVector& posRef,
                                               std::map<G4ThreeVector,G4double>* tarMap,
                                               G4String volName,
                                               G4bool in)
{
    // return solidOrbRef;
    G4SubtractionSolid* solidCut(NULL); // container for the cut solid
    // G4UnionSolid* solidCut(NULL); // container for the cut solid

    bool isCutted = false; // flag to indicate if volume has been cut yet
    bool isOurVol = false; // flag to indicate if target volume is our reference volume

    //----------------------------------------------------------------------------------------------
    // Define radius of the reference solid we are focusing on
    //----------------------------------------------------------------------------------------------
    G4double radiusRef (0);
    // if(volName=="base") radiusRef = fBaseRadiusWater;
    // else if(volName=="sugarTHF") radiusRef = fSugarTHFRadiusWater;
    // else if(volName=="sugarTMP") radiusRef = fSugarTMPRadiusWater;
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
            // G4cout << "The matching index is: " << count << G4endl;
            // target volume is the reference volume
            // G4cout << "Distance is 0 at index " << count << G4endl;
            // G4cout << "Radius = " << radiusTar << G4endl;
            // G4cout << "Position = " << posTar << G4endl;
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
            intersection -= 0.001*fFactor*nm;
            // Create a vector to the intersection position, where one edge of the slicing volume
            // will be placed 
            G4ThreeVector posSlice = intersection *( displacement_vector/displacement_vector.getR() );

            //--------------------------------------------------------------------------------------
            // Calculate the necessary rotations.
            //--------------------------------------------------------------------------------------
            G4double phi = std::acos(posSlice.getZ()/posSlice.getR());
            G4double theta = std::acos( posSlice.getX() / ( posSlice.getR()*std::cos(M_PI/2.-phi) ) );

            if(posSlice.getY()<0) theta = -theta;

            G4ThreeVector rotAxisForPhi(1*fFactor*nm,0.,0.);
            rotAxisForPhi.rotateZ(theta+M_PI/2);
            G4RotationMatrix *rotMat = new G4RotationMatrix;
            rotMat->rotate(-phi, rotAxisForPhi);

            G4ThreeVector rotZAxis(0.,0.,1*fFactor*nm);
            rotMat->rotate(theta, rotZAxis);

            //--------------------------------------------------------------------------------------
            // Create the G4SubtractionSolid.
            //--------------------------------------------------------------------------------------
            if(!isCutted) solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, sliceBox, rotMat, posSlice);
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, sliceBox, rotMat, posSlice);
            // if(!isCutted) solidCut = new G4UnionSolid("solidCut", solidOrbRef, sliceBox, rotMat, posSlice);
            // else solidCut = new G4UnionSolid("solidCut", solidCut, sliceBox, rotMat, posSlice);

            isCutted = true;
        }
        count++;
    }

    if(isCutted) return solidCut;
    else return solidOrbRef;
}

//--------------------------------------------------------------------------------------------------
// Helper function used by buildLogicFiber() in verbose mode (i.e. if want to output the mean volume
// of cut residue and water volumes).
//--------------------------------------------------------------------------------------------------
void GeoVolumeV2::CalculateMeanVol(std::map<G4String, std::vector<G4LogicalVolume*> >* logicSolidsMap)
{
    G4double sugarTMP1Vol = 0;
    G4double sugarTHF1Vol = 0;
    G4double base1Vol = 0;
    G4double base2Vol = 0;
    G4double sugarTHF2Vol = 0;
    G4double sugarTMP2Vol = 0;

    G4double sugarTMP1WaterVol = 0;
    G4double sugarTHF1WaterVol = 0;

    G4double base1WaterVol = 0;
    G4double base2WaterVol = 0;
    G4double sugarTHF2WaterVol = 0;
    G4double sugarTMP2WaterVol = 0;

    // iterate on each bp
    for(int j=0;j<fBpNum;++j)
    {
        if(j==0) G4cout<<"Volume calculations\nThe fFactor value is taken into account already. We assume fFactor=1e+9."<<G4endl;
        // assert(fFactor==1e+9);
        sugarTMP1Vol += logicSolidsMap->at("sugarTMP1")[j]->GetSolid()->GetCubicVolume();
        sugarTHF1Vol += logicSolidsMap->at("sugarTHF1")[j]->GetSolid()->GetCubicVolume();
        base1Vol += logicSolidsMap->at("base1")[j]->GetSolid()->GetCubicVolume();
        base2Vol += logicSolidsMap->at("base2")[j]->GetSolid()->GetCubicVolume();
        sugarTHF2Vol += logicSolidsMap->at("sugarTHF2")[j]->GetSolid()->GetCubicVolume();
        sugarTMP2Vol += logicSolidsMap->at("sugarTMP2")[j]->GetSolid()->GetCubicVolume();

        // sugarTMP1WaterVol += logicSolidsMap->at("sugarTMPWater1")[j]->GetSolid()->GetCubicVolume();
        // sugarTHF1WaterVol += logicSolidsMap->at("sugarTHFWater1")[j]->GetSolid()->GetCubicVolume();
        // base1WaterVol += logicSolidsMap->at("baseWater1")[j]->GetSolid()->GetCubicVolume();
        // base2WaterVol += logicSolidsMap->at("baseWater2")[j]->GetSolid()->GetCubicVolume();
        // sugarTHF2WaterVol += logicSolidsMap->at("sugarTHFWater2")[j]->GetSolid()->GetCubicVolume();
        // sugarTMP2WaterVol += logicSolidsMap->at("sugarTMPWater2")[j]->GetSolid()->GetCubicVolume();
    }

    sugarTMP1Vol = sugarTMP1Vol/(fBpNum);
    sugarTHF1Vol = sugarTHF1Vol/(fBpNum);
    base1Vol = base1Vol/(fBpNum);
    base2Vol = base2Vol/(fBpNum);
    sugarTHF2Vol = sugarTHF2Vol/(fBpNum);
    sugarTMP2Vol = sugarTMP2Vol/(fBpNum);

    // sugarTMP1Vol = sugarTMP1Vol/(fBpNum*pow(fFactor,3) );
    // sugarTHF1Vol = sugarTHF1Vol/(fBpNum*pow(fFactor,3) );
    // base1Vol = base1Vol/(fBpNum*pow(fFactor,3));
    // base2Vol = base2Vol/(fBpNum*pow(fFactor,3));
    // sugarTHF2Vol = sugarTHF2Vol/(fBpNum*pow(fFactor,3));
    // sugarTMP2Vol = sugarTMP2Vol/(fBpNum*pow(fFactor,3));

    // sugarTMP1WaterVol = sugarTMP1WaterVol/(fBpNum*pow(fFactor,3) ) - sugarTMP1Vol;
    // sugarTHF1WaterVol = sugarTHF1WaterVol/(fBpNum*pow(fFactor,3) ) - sugarTHF1Vol;
    // base1WaterVol = base1WaterVol/(fBpNum*pow(fFactor,3)) - base1Vol;
    // base2WaterVol = base2WaterVol/(fBpNum*pow(fFactor,3)) - base2Vol;
    // sugarTHF2WaterVol = sugarTHF2WaterVol/(fBpNum*pow(fFactor,3)) - sugarTHF2Vol;
    // sugarTMP2WaterVol = sugarTMP2WaterVol/(fBpNum*pow(fFactor,3)) - sugarTMP2Vol;

    G4cout<<"\nsugarTMP1Vol="<<sugarTMP1Vol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"\nsugarTHF1Vol="<<sugarTHF1Vol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"base1Vol="<<base1Vol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"base2Vol="<<base2Vol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"sugarTHF2Vol="<<sugarTHF2Vol/m3*1e+27<<" nm3\n"<<G4endl;
    G4cout<<"sugarTMP2Vol="<<sugarTMP2Vol/m3*1e+27<<" nm3\n"<<G4endl;

    // G4cout<<"\nsugarTMP1WaterVol="<<sugarTMP1WaterVol/m3*1e+27<<" nm3"<<G4endl;
    // G4cout<<"\nsugarTHF1WaterVol="<<sugarTHF1WaterVol/m3*1e+27<<" nm3"<<G4endl;
    // G4cout<<"base1WaterVol="<<base1WaterVol/m3*1e+27<<" nm3"<<G4endl;
    // G4cout<<"base2WaterVol="<<base2WaterVol/m3*1e+27<<" nm3"<<G4endl;
    // G4cout<<"sugarTHF2WaterVol="<<sugarTHF2WaterVol/m3*1e+27<<" nm3\n"<<G4endl;
    // G4cout<<"sugarTMP2WaterVol="<<sugarTMP2WaterVol/m3*1e+27<<" nm3\n"<<G4endl;
}

//--------------------------------------------------------------------------------------------------
// Helper function used if a geometrical overlap is detected. Return the standard Topas message
// about geometry overlaps and exit gracefully.
//--------------------------------------------------------------------------------------------------
void GeoVolumeV2::ThrowOverlapError() 
{
    G4cout << "Topas is quitting due to the above geometry overlap problem." << G4endl;
    G4cout << "If you still want the TOPAS session to continue" << G4endl;
    G4cout << "(such as to use visualization to study the overlap)," << G4endl;
    G4cout << "Set the parameter Ge/QuitIfOverlapDetected to False" << G4endl;
    exit(0);
}
