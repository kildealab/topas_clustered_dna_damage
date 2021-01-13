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

    G4double radius1 = 10*nm;
    G4double radius2 = 5*nm;

    G4Orb* sphere1Solid = new G4Orb("sphere1",radius1);
    G4Orb* sphere2Solid = new G4Orb("sphere2",radius2);

    G4LogicalVolume* sphere1Log = CreateLogicalVolume("log1", sphere1Solid);
    G4LogicalVolume* sphere2Log = CreateLogicalVolume("log2", sphere2Solid);

    G4VisAttributes sphere1Vis(G4Colour(0.0, 1.0, 0.0, 0.3) );
    sphere1Vis.SetVisibility(true);
    sphere1Vis.SetForceSolid(true);
    sphere1Log->SetVisAttributes(sphere1Vis);

    G4VisAttributes sphere2Vis(G4Colour(0.0, 0.0, 1.0, 0.3) );
    sphere2Vis.SetVisibility(true);
    sphere2Vis.SetForceSolid(true);
    sphere2Log->SetVisAttributes(sphere2Vis);

    G4ThreeVector pos1 = G4ThreeVector(0.,0.,0.);
    // G4ThreeVector pos2 = G4ThreeVector(0.,0.,radius1+(radius2/2));
    // G4ThreeVector pos2 = G4ThreeVector(0.,(radius1+(radius2/2))/sqrt(2),(radius1+(radius2/2))/sqrt(2));
    G4ThreeVector pos2 = G4ThreeVector((radius1+(radius2/2))/sqrt(2),(radius1+(radius2/2))/sqrt(2),0.);
    // G4ThreeVector pos2 = G4ThreeVector((radius1+(radius2/2))/sqrt(3),(radius1+(radius2/2))/sqrt(3),(radius1+(radius2/2))/sqrt(3));

    G4RotationMatrix* rot1 = new G4RotationMatrix(0,0,0);
    G4RotationMatrix* rot2 = new G4RotationMatrix(0,0,0);

    // CreatePhysicalVolume("phys1", sphere1Log, fEnvelopePhys);

    CreatePhysicalVolume("phys1", sphere1Log, rot1, &pos1, fEnvelopePhys);
    CreatePhysicalVolume("phys2", sphere2Log, rot2, &pos2, fEnvelopePhys);



    G4Box* sliceBox = new G4Box("cutoff_box",radius2, radius2, radius2);

    G4ThreeVector displacement_vector = pos2 - pos1;
    G4double displacement = displacement_vector.getR();
    G4double crossover = ((pow(radius1,2) + pow(displacement,2) - pow(radius2,2))/(2*displacement)) + sliceBox->GetZHalfLength();
    G4ThreeVector posSlice = crossover * (displacement_vector/displacement);

    G4LogicalVolume* sliceLog = CreateLogicalVolume("sliceLog", sliceBox);
    G4VisAttributes sliceVis(G4Colour(1.0, 0.0, 0.0, 0.3) );
    sliceVis.SetVisibility(true);
    sliceVis.SetForceSolid(true);
    sliceLog->SetVisAttributes(sliceVis);

    G4double theta = std::atan(displacement_vector.getY()/displacement_vector.getX());
    G4double phi = std::acos(displacement_vector.getZ()/displacement);
    G4RotationMatrix* rotSlice = new G4RotationMatrix(theta,phi,0);

    CreatePhysicalVolume("slice", sliceLog, rotSlice, &posSlice, fEnvelopePhys);




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





 //    G4cout << "Building the Fiber" << G4endl;

 //    G4int numBpPerNucleosome = fPm->GetIntegerParameter(GetFullParmName("DNANumBpPerNucleosome"));
 //    G4int numNucleosomePerFiber = fPm->GetIntegerParameter(GetFullParmName("DnaNumNucleosomePerFiber"));
 //    fGeoManager->Initialize(numBpPerNucleosome,numNucleosomePerFiber);
 //    // fGeoManager->Initialize();
    
 //    G4LogicalVolume* lFiber = fGeoManager->BuildLogicFiber(false);
    
 //    CreatePhysicalVolume("Fiber", lFiber, fEnvelopePhys);
    
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

