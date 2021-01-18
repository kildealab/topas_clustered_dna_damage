// Component for TsFiberV2
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
// Code that reads in the DNA fiber geometry from the DNAFabric software package tool
// More details about the fiber geometry can be found in:
// Meylan et al. (2016) Commputer Physics Comunications, 204, 159
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:


#include "TsFiberV2.hh"
#include "GeoManagerV2.hh"

#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4UserLimits.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4NistManager.hh"

TsFiberV2::TsFiberV2(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    
    fGeoManager = new GeoManagerV2(0, 1.);
    
}


TsFiberV2::~TsFiberV2()
{
     delete fGeoManager;
}


G4VPhysicalVolume* TsFiberV2::Construct()
{
	BeginConstruction();

    G4double* envelopeDimensions = fPm->GetDoubleVector(GetFullParmName("Dimensions"),"Length");
    
    //Wrapping component for whole fiber    
    G4Tubs* sWrapper = new G4Tubs("solid_wrapper", 0., envelopeDimensions[0], envelopeDimensions[1], 0, 360);
    fEnvelopeLog = CreateLogicalVolume(sWrapper);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

    // //**********************************************************************************************
    // // Verification of cutting solid behaviour
    // //----------------------------------------------------------------------------------------------
    // // Sphere 1 & 2 solids & logicals 
    // //----------------------------------------------------------------------------------------------
    // // G4double radius1 = 10*nm;
    // // G4double radius2 = 5*nm;
    // G4double radius1 = 0.27*nm;
    // G4double radius2 = 0.29*nm;

    // G4Orb* sphere1Solid = new G4Orb("sphere1",radius1);
    // G4Orb* sphere2Solid = new G4Orb("sphere2",radius2);

    // G4LogicalVolume* sphere1Log = CreateLogicalVolume("log1", sphere1Solid);
    // G4LogicalVolume* sphere2Log = CreateLogicalVolume("log2", sphere2Solid);

    // G4VisAttributes sphere1Vis(G4Colour(0.0, 1.0, 0.0, 0.3) );
    // sphere1Vis.SetVisibility(true);
    // sphere1Vis.SetForceSolid(true);
    // sphere1Log->SetVisAttributes(sphere1Vis);

    // G4VisAttributes sphere2Vis(G4Colour(0.0, 0.0, 1.0, 0.3) );
    // sphere2Vis.SetVisibility(true);
    // sphere2Vis.SetForceSolid(true);
    // sphere2Log->SetVisAttributes(sphere2Vis);


    // //----------------------------------------------------------------------------------------------
    // // Various ways of doing positioning - uncomment diff pos2 settings to test 
    // //----------------------------------------------------------------------------------------------
    // // G4ThreeVector pos1 = G4ThreeVector(0.,0.,0.);
    // // G4ThreeVector pos2 = G4ThreeVector(0.,0.,radius1+(radius2/2));
    // // G4ThreeVector pos2 = G4ThreeVector(0.,(radius1+(radius2/2))/sqrt(2),(radius1+(radius2/2))/sqrt(2));
    // // G4ThreeVector pos2 = G4ThreeVector((radius1+(radius2/2))/sqrt(2),(radius1+(radius2/2))/sqrt(2),0.);
    // // G4ThreeVector pos2 = G4ThreeVector((radius1+(radius2/2))/sqrt(3),(radius1+(radius2/2))/sqrt(3),(radius1+(radius2/2))/sqrt(3));

    // // other octants:
    // // G4ThreeVector pos2 = G4ThreeVector(0.,0.,-1*(radius1+(radius2/2)));
    // // G4ThreeVector pos2 = G4ThreeVector(0.,-1*(radius1+(radius2/2))/sqrt(2),(radius1+(radius2/2))/sqrt(2));
    // // G4ThreeVector pos2 = G4ThreeVector(-(radius1+(radius2/2))/sqrt(3),-(radius1+(radius2/2))/sqrt(3),-(radius1+(radius2/2))/sqrt(3));

    // // other angles 
    // // G4ThreeVector pos2 = G4ThreeVector(0,radius2/2,radius1+(radius2/4));

    // // match actual sim
    // G4ThreeVector pos1 = G4ThreeVector(7.86727*nm,13.2035*nm,1.2072*nm);
    // G4ThreeVector pos2 = G4ThreeVector(7.8662*nm,13.1051*nm,1.5331*nm);


    // //----------------------------------------------------------------------------------------------
    // // offset positions from origin (include a sphere @ origin for scale)
    // //----------------------------------------------------------------------------------------------
    // G4double deltaX = 0*nm;
    // G4double deltaY = 0*nm;
    // G4double deltaZ = 0*nm;
    // G4Orb* originSolid = new G4Orb("sphere1",0.5*nm);
    // G4LogicalVolume* originLog = CreateLogicalVolume("originLog", originSolid);
    // G4VisAttributes originVis(G4Colour(1.0, 1.0, 1.0) );
    // sphere1Vis.SetVisibility(true);
    // sphere1Vis.SetForceSolid(true);
    // sphere1Log->SetVisAttributes(sphere1Vis);
    // G4RotationMatrix* rotOrigin = new G4RotationMatrix(0,0,0);
    // G4ThreeVector posOrigin = G4ThreeVector(0.,0.,0.);
    // // CreatePhysicalVolume("originPhys", originLog, rotOrigin, &posOrigin, fEnvelopePhys);

    // pos1.setX(pos1.getX()+deltaX);
    // pos1.setY(pos1.getY()+deltaY);
    // pos1.setZ(pos1.getZ()+deltaZ);

    // pos2.setX(pos2.getX()+deltaX);
    // pos2.setY(pos2.getY()+deltaY);
    // pos2.setZ(pos2.getZ()+deltaZ);

    // //----------------------------------------------------------------------------------------------
    // // rotation matrices for spheres 1 & 2 (should be zero)
    // //----------------------------------------------------------------------------------------------
    // G4RotationMatrix* rot1 = new G4RotationMatrix(0,0,0);
    // G4RotationMatrix* rot2 = new G4RotationMatrix(0,0,0);


    // //----------------------------------------------------------------------------------------------
    // // solid & logical volumes for the slicing box
    // //----------------------------------------------------------------------------------------------
    // G4Box* sliceBox = new G4Box("cutoff_box",2*radius2, 2*radius2, 2*radius2);

    // // Positioning information for sliceBox
    // G4ThreeVector displacement_vector = pos2 - pos1;
    // G4double displacement = displacement_vector.getR();
    // // G4double crossover = ((pow(radius1,2) + pow(displacement,2) - pow(radius2,2))/(2*displacement)) + sliceBox->GetZHalfLength();
    // G4double crossover = ((pow(radius1,2) + pow(displacement,2) - pow(radius2,2))/(2*displacement)) + sliceBox->GetZHalfLength() - 0.001*nm;;
    // G4ThreeVector posSlice = crossover * (displacement_vector/displacement);

    // G4LogicalVolume* sliceLog = CreateLogicalVolume("sliceLog", sliceBox);
    // G4VisAttributes sliceVis(G4Colour(1.0, 0.0, 0.0, 0.3) );
    // sliceVis.SetVisibility(true);
    // sliceVis.SetForceSolid(true);
    // sliceLog->SetVisAttributes(sliceVis);

    // //----------------------------------------------------------------------------------------------
    // // My first attempt at rotating the slicing box
    // //----------------------------------------------------------------------------------------------
    // // G4double theta = std::atan(displacement_vector.getY()/displacement_vector.getX());
    // // G4double phi = std::acos(displacement_vector.getZ()/displacement);
    // // G4RotationMatrix* rotSlice = new G4RotationMatrix(theta,phi,0);
    // // G4RotationMatrix* rotSlice = new G4RotationMatrix(displacement_vector,phi);

    // //----------------------------------------------------------------------------------------------
    // // rotation of the slicing box
    // //----------------------------------------------------------------------------------------------
    // G4double phi = std::acos(posSlice.getZ()/posSlice.getR());
    // G4double theta = std::acos( posSlice.getX() / ( posSlice.getR()*std::cos(M_PI/2.-phi) ) );

    // if(posSlice.getY()<0) theta = -theta;

    // G4ThreeVector rotAxisForPhi(1*nm,0.,0.);
    // rotAxisForPhi.rotateZ(theta+M_PI/2);
    // G4RotationMatrix *rotMat = new G4RotationMatrix;
    // rotMat->rotate(-phi, rotAxisForPhi);

    // G4ThreeVector rotZAxis(0.,0.,1*nm);
    // rotMat->rotate(theta, rotZAxis);

    // //----------------------------------------------------------------------------------------------
    // // Create a union solid
    // //----------------------------------------------------------------------------------------------
    // G4UnionSolid* unionSolid = new G4UnionSolid("union",sphere1Solid,sliceBox,rotMat,posSlice);
    // G4LogicalVolume* unionLog = CreateLogicalVolume("unionLog", unionSolid);
    // unionLog->SetVisAttributes(sphere1Vis);

    // // Dont apply positioning offset from origin for slicing box until after rotation has been done
    // // and a union solid has been made
    // posSlice.setX(posSlice.getX()+deltaX);
    // posSlice.setY(posSlice.getY()+deltaY);
    // posSlice.setZ(posSlice.getZ()+deltaZ);


    // //----------------------------------------------------------------------------------------------
    // // Second slicing box
    // //----------------------------------------------------------------------------------------------
    // // G4Box* sliceBox2 = new G4Box("cutoff_box",2*radius2, 2*radius2, 2*radius2);

    // // // Positioning information for sliceBox
    // // G4ThreeVector displacement_vector2 = pos1 - pos2;
    // // G4double displacement2 = displacement_vector2.getR();
    // // // G4double crossover = ((pow(radius1,2) + pow(displacement,2) - pow(radius2,2))/(2*displacement)) + sliceBox->GetZHalfLength();
    // // G4double crossover2 = ((pow(radius2,2) + pow(displacement2,2) - pow(radius1,2))/(2*displacement2)) + sliceBox2->GetZHalfLength() - 0.001*nm;;
    // // G4ThreeVector posSlice2 = crossover2 * (displacement_vector2/displacement2);

    // // G4LogicalVolume* sliceLog2 = CreateLogicalVolume("sliceLog2", sliceBox);
    // // G4VisAttributes sliceVis2(G4Colour(1.0, 0.0, 0.0, 0.3) );
    // // sliceVis2.SetVisibility(true);
    // // sliceVis2.SetForceSolid(true);
    // // sliceLog2->SetVisAttributes(sliceVis2);

    // // G4double phi2 = std::acos(posSlice2.getZ()/posSlice2.getR());
    // // G4double theta2 = std::acos( posSlice2.getX() / ( posSlice2.getR()*std::cos(M_PI/2.-phi) ) );

    // // if(posSlice2.getY()<0) theta2 = -theta2;

    // // G4ThreeVector rotAxisForPhi2(1*nm,0.,0.);
    // // rotAxisForPhi2.rotateZ(theta2+M_PI/2);
    // // G4RotationMatrix *rotMat2 = new G4RotationMatrix;
    // // rotMat2->rotate(-phi2, rotAxisForPhi2);

    // // G4ThreeVector rotZAxis2(0.,0.,1*nm);
    // // rotMat2->rotate(theta2, rotZAxis2);

    // // G4UnionSolid* unionSolid2 = new G4UnionSolid("union2",sphere2Solid,sliceBox2,rotMat2,posSlice2);
    // // G4LogicalVolume* unionLog2 = CreateLogicalVolume("unionLog2", unionSolid2);
    // // unionLog2->SetVisAttributes(sphere2Vis);

    // // posSlice2.setX(posSlice2.getX()+deltaX);
    // // posSlice2.setY(posSlice2.getY()+deltaY);
    // // posSlice2.setZ(posSlice2.getZ()+deltaZ);

    // //----------------------------------------------------------------------------------------------
    // // Create & place physical volumes
    // //----------------------------------------------------------------------------------------------
    // // CreatePhysicalVolume("phys1", sphere1Log, rot1, &pos1, fEnvelopePhys);
    // // CreatePhysicalVolume("phys2", sphere2Log, rot2, &pos2, fEnvelopePhys);

    // // CreatePhysicalVolume("slice", sliceLog, rotSlice, &posSlice, fEnvelopePhys);
    // // CreatePhysicalVolume("slice", sliceLog, rotMat, &posSlice, fEnvelopePhys);

    // // Topas way of doing physical volumes (1 union, 1 normal)
    // CreatePhysicalVolume("union", unionLog, rot1, &pos1, fEnvelopePhys);
    // CreatePhysicalVolume("phys2", sphere2Log, rot2, &pos2, fEnvelopePhys);

    // // G4 way of doing physical volumes instead of topas:
    // // G4PVPlacement* phys1 = new G4PVPlacement(rot1,pos1,unionLog,"phys1",fEnvelopeLog,false,0);
    // // G4PVPlacement* phys2 = new G4PVPlacement(rot2,pos2,sphere2Log,"phys2",fEnvelopeLog,false,0);

    // // Create physical volume for second union volume:
    // // CreatePhysicalVolume("union2", unionLog2, rot2, &pos2, fEnvelopePhys);
    // //**********************************************************************************************




    G4cout << "Building the Fiber" << G4endl;

    G4int numBpPerNucleosome = fPm->GetIntegerParameter(GetFullParmName("DNANumBpPerNucleosome"));
    G4int numNucleosomePerFiber = fPm->GetIntegerParameter(GetFullParmName("DnaNumNucleosomePerFiber"));
    fGeoManager->Initialize(numBpPerNucleosome,numNucleosomePerFiber);
    // fGeoManager->Initialize();
    
    G4LogicalVolume* lFiber = fGeoManager->BuildLogicFiber(false);
    
    CreatePhysicalVolume("Fiber", lFiber, fEnvelopePhys);
    
	return fEnvelopePhys;
}

G4Material * TsFiberV2::OtherMaterial(G4String materialName)
{
    G4Material * material(0);
    
    if(materialName == "G4_WATER"){
        // Water is defined from NIST material database
        G4NistManager * man = G4NistManager::Instance();
        material = man->FindOrBuildMaterial("G4_WATER");

    }
    
    return material;
}

// // G4Tubs* boxOuter = new G4Tubs("outer_box", 0., 11*nm, 20*nm, 0., 360);
// G4Box* boxOuter = new G4Box("outer_box",20*nm,20*nm,20*nm);
// // G4Orb* boxOuter = new G4Orb("outer_box",10*nm);

// G4Box* boxInner = new G4Box("inner_box",20*nm,20*nm,10*nm);
// // G4Orb* boxInner = new G4Orb("outer_box",5*nm);
// // G4Box* boxInner = new G4Box("inner_box",20*nm,20*nm,5*nm);

// G4RotationMatrix *rotMat = new G4RotationMatrix;
// G4SubtractionSolid* boxSubtract = new G4SubtractionSolid("subtraction_box",boxOuter,boxInner,0,G4ThreeVector(0.,0.,10*nm));
// // G4UnionSolid* boxSubtract = new G4UnionSolid("union_box",boxOuter,boxInner,0,G4ThreeVector(0.,0.,10*nm));
// // G4Material* fWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

// // G4LogicalVolume* logicalBox = new G4LogicalVolume(boxSubtract,fWater,"boxOuter");
// // G4PVPlacement* sTMP1 = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicalBox,
// //     "physical_box",fEnvelopeLog,false,1);


// G4LogicalVolume* logicalBox = CreateLogicalVolume("logical_box",boxSubtract);
// G4VisAttributes boxVis(G4Colour(0.0, 1.0, 0.0, 0.5) );
// boxVis.SetVisibility(true);
// boxVis.SetForceSolid(true);
// // boxVis.SetForceWireframe(true);
// logicalBox->SetVisAttributes(boxVis);

// G4double subtractVolume = boxSubtract->GetCubicVolume();
// G4double outerVolume = boxOuter->GetCubicVolume();
// G4double innerVolume = boxInner->GetCubicVolume();
// G4cout << "---------------------------------------------------------" << G4endl;
// G4cout << "Outer volume: " << outerVolume/m3 << G4endl;
// G4cout << "Inner volume: " << innerVolume/m3 << G4endl;
// G4cout << "Subtracted volume: " << subtractVolume/m3 << G4endl;
// G4cout << "---------------------------------------------------------" << G4endl;

// CreatePhysicalVolume("physical_box", logicalBox, fEnvelopePhys);
