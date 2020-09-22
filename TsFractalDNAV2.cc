// Component for TsFractalDNAV2
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

#include "TsFractalDNAV2.hh"

#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4TwoVector.hh"


#include "G4Box.hh"
#include "G4Ellipsoid.hh"

#include "G4LogicalVolume.hh"

#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"

#include "G4VisAttributes.hh"

#include "GeoManagerV2.hh"

#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <chrono>

using namespace std;

//--------------------------------------------------------------------------------------------------
// Constructor: 
// Calls ResolveParameters method to _____
// Creates a GeoManager object
//--------------------------------------------------------------------------------------------------
TsFractalDNAV2::TsFractalDNAV2(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
    fGeoManager = new GeoManagerV2(0, 1.); // parameters for verbosity & scaling ffactor
}

//--------------------------------------------------------------------------------------------------
// Destructor
//--------------------------------------------------------------------------------------------------
TsFractalDNAV2::~TsFractalDNAV2()
{
    delete fGeoManager;
    
}

//--------------------------------------------------------------------------------------------------
// Handle parameters defined in the Topas parameter file
//--------------------------------------------------------------------------------------------------
void TsFractalDNAV2::ResolveParameters(){
    // User defined parameters.
    
    // File containing fibre placement coordinates
    G4String name = GetFullParmName("FileName");
    if (!fPm->ParameterExists(name)) {
        G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
        G4cerr << "Parameter " << name << " has to be specified for this scorer." << G4endl;
        exit(1);
    }
    G4String FileName = fPm->GetStringParameter(name);
    
    // Boolean indicating whether individual bases are constructed
    fBuildBases = fPm->GetBooleanParameter(GetFullParmName("BuildBases"));
    if (fBuildBases){
        G4cout << "Topas is building Bases." << G4endl;
        G4cout << "This can slow down simulations with visualization" << G4endl;
    }
    
    //Hilbert space filling:
    const char* filename = FileName;
    std::string line = "";
    G4double x, y, z;
    
    ifstream f(filename, ios::in);
    
    // Read through file containing xyz coordinates, save to G4Double vectors fx, fy, and fz
    if (f.is_open())
    {
        while (!f.eof())
        {
            f >> x >> y >> z;
            fx.push_back(x);
            fy.push_back(y);
            fz.push_back(z);
        
        }
    }
    else
    {
        G4cout << "ERROR: Unable to open file " << FileName << G4endl;
        exit(1);
    }
    f.close();
}


//--------------------------------------------------------------------------------------------------
// Mandatory method of TSVGeometryComponent to create solids, logical volumes, and physical volumes.
// Primary Logical Volume, fEnvelopeLog, is defined as the nucleus, which may be accessed in the 
// parameter file as Ge/MyComponent.
//--------------------------------------------------------------------------------------------------
G4VPhysicalVolume* TsFractalDNAV2::Construct()
{
	BeginConstruction();

    /**************************************************************************/
    //                 Wrapper: Nucleus
    /**************************************************************************/
    
    G4Ellipsoid* fEnvelopeSolid = new G4Ellipsoid("nucl", 6*um, 4*um, 11*um);
    fEnvelopeLog= CreateLogicalVolume(fEnvelopeSolid);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

    /**************************************************************************/
    //                 Subcomponent 1: chromatin fibers
    /**************************************************************************/
    
    // Define fibre dimensions
    G4double fiberRadius = 12*nm;
    G4double length = std::sqrt( pow(fx[2]-fx[1],2) + pow(fy[2]-fy[1],2) + pow(fz[2]-fz[1],2));
    G4double scaleFactor = (170*nm)/length;  // scaleFactor used to get fiber lengths of 170*nm

    length=scaleFactor*length;

    // Create G4 solid: Only one solid fibre is created
    G4String subComponentName1 = "Chromoloop";
    G4Tubs* gLoop = new G4Tubs(subComponentName1, 0, fiberRadius, (length/2)-fiberRadius, 0*deg, 360*deg); // Why is fiberRadius subtracted from length?

    // 23 chromosome territories in humans (22 + X or Y)
    // LM - actually seems like 24 territories below, but 22 & Y are grouped together? Should be X & Y grouped together?
    // Each territory is described by its own logical volume contained in lLoop vector
    lLoop.resize(200000);
    sLoop.resize(200000);

    G4VisAttributes* LoopVisAtta = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
    LoopVisAtta->SetForceSolid(true);

    for (G4int j = 1; j < 200000; j++){
        // sLoop[j] = new G4Tubs(subComponentName1, 0, fiberRadius, (length/2)-fiberRadius, 0*deg, 360*deg);
        // lLoop[j] = CreateLogicalVolume(sLoop[j]);
        lLoop[j] = CreateLogicalVolume(gLoop);
        lLoop[j]->SetVisAttributes(LoopVisAtta);
    }

    G4LogicalVolume* myLoop = CreateLogicalVolume(gLoop);
    myLoop->SetVisAttributes(LoopVisAtta);

    // G4VisAttributes* LoopVisAtta = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
    // LoopVisAtta->SetForceWireframe(true);
    // LoopVisAtta->SetForceSolid(true);
    // lmLoop.resize(15);
    // for (G4int j = 0; j < 15; j++){
    //     lmLoop[j] = CreateLogicalVolume(gLoop);
    //     lmLoop[j]->SetVisAttributes(LoopVisAtta);
    // }
    // G4VisAttributes* LoopVisAtta = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
    // LoopVisAtta->SetForceSolid(true);
    // lmLoop[0]->SetVisAttributes(LoopVisAtta);

    // G4VisAttributes* LoopVisAttb = new G4VisAttributes(G4Colour(0.0,1.0,0.0)); //green
    // LoopVisAttb->SetForceSolid(true);
    // lmLoop[1]->SetVisAttributes(LoopVisAttb);

    // G4VisAttributes* LoopVisAttc = new G4VisAttributes(G4Colour(0.0,0.0,1.0)); //blue
    // LoopVisAttc->SetForceSolid(true);
    // lmLoop[2]->SetVisAttributes(LoopVisAttc);

    // G4VisAttributes* LoopVisAttd = new G4VisAttributes(G4Colour(1.0,1.0,0.0)); //yellow
    // LoopVisAttd->SetForceSolid(true);
    // lmLoop[3]->SetVisAttributes(LoopVisAttd);

    // G4VisAttributes* LoopVisAtte = new G4VisAttributes(G4Colour(1.0,0.0,1.0)); //magenta
    // LoopVisAtte->SetForceSolid(true);
    // lmLoop[4]->SetVisAttributes(LoopVisAtte);

    // G4VisAttributes* LoopVisAttf = new G4VisAttributes(G4Colour(0.0,1.0,1.0)); //teal
    // LoopVisAttf->SetForceSolid(true);
    // lmLoop[5]->SetVisAttributes(LoopVisAttf);

    // G4VisAttributes* LoopVisAttg = new G4VisAttributes(G4Colour(0.75,0.75,0.75)); //gray
    // LoopVisAttg->SetForceSolid(true);
    // lmLoop[6]->SetVisAttributes(LoopVisAttg);


    // TP = total number of fibres (number of lines in data file - 1)
    G4int TotalPoints = fx.size();
    const G4int TP = TotalPoints - 1;

    pLoop.resize(TP); // vector of physical volumes

    G4cout << "Number of Fibers " << TP << G4endl; // This is not true, TP is one more than # of fibres
    G4double xshift = 0, yshift = 0, zshift = 0;
    

    //----------------------------------------------------------------------------------------------
    // Build bases here (once)
    //----------------------------------------------------------------------------------------------
    // if (fBuildBases){

    auto start = std::chrono::high_resolution_clock::now();
    fGeoManager->Initialize();
    auto finish = std::chrono::high_resolution_clock::now();
    std:chrono::duration<double> elapsed = finish-start;
    G4cout << "Time to Initialize GeoManager: " << elapsed.count() << " s" << G4endl;  
    start = std::chrono::high_resolution_clock::now();  
    G4LogicalVolume* lFiber = fGeoManager->BuildLogicFiber(true);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish-start;
    G4cout << "Time to build logic fiber: " << elapsed.count() << " s" << G4endl;  
    // G4cout << elapsed.count() << G4endl;       

    // }

    // Full genome model of 172032 fibers with ~6 billion BPs

    // Fill the ellipsoid nucleus with the fibers, using a Hilbert space filling curve:
    // No fibres are created until at least 3 coordinates in .dat file
    // Last entry is never used, hence num_fibres = TP - 1
    for (G4int i = 1; i < TP; i++) {
        // G4cout << "--------------------------------" << G4endl;
        // G4cout << i << G4endl;
        auto loopstart = std::chrono::high_resolution_clock::now();
        // Don't act on the 6 fibres indicated on the next line... unsure why
        // There are discontinuities at the 6 fibre indices labeled below. Unsure why, perhaps the
        // fractal space filling used to generate the .dat file was unable to generate the 
        // longitudinal extent required for this work? The pattern in FullGenome.dat repeats at
        // line 32769, 65537, etc. (i.e. resets to x=0.48438, y=0.48438, z=0.48438) but new zshift
        // is applied as below
        // It seems like the fibres that would exist at the following indices are not generated at all.
        if ((i != 32768) || (i != 2*32768) || (i != 3*32768) || (i != 4*32768) || (i != 5*32768) || (i != 5*32768+4096)) {
            
            if ((i >= 32768) && (i<= 2*32768)) {
                zshift = 2.263*um;
                length = std::sqrt( pow(fx[i+1]-fx[i],2) + pow(fy[i+1]-fy[i],2) + pow(fz[i+1]-fz[i],2));
                scaleFactor = (170*nm)/length;
            }
            
            if ((i >= 2*32768+1) && (i<= 3*32768)) {
                zshift = -2.263*um;
                length = std::sqrt( pow(fx[i+1]-fx[i],2) + pow(fy[i+1]-fy[i],2) + pow(fz[i+1]-fz[i],2));
                scaleFactor = (170*nm)/length;
            }
            
            if ((i >= 3*32768+1) && (i<= 4*32768)) {
                zshift = 4.526*um;
                length = std::sqrt( pow(fx[i+1]-fx[i],2) + pow(fy[i+1]-fy[i],2) + pow(fz[i+1]-fz[i],2));
                scaleFactor = (170*nm)/length;
            }
            if ((i >= 4*32768+1) && (i <= 5*32768)) {
                zshift = -4.526*um;
                length = std::sqrt( pow(fx[i+1]-fx[i],2) + pow(fy[i+1]-fy[i],2) + pow(fz[i+1]-fz[i],2));
                scaleFactor = (170*nm)/length;
            }
            if ((i >= 5*32768+1) && (i <= 5*32768+4096)) {
                zshift = -8.5*um;
                length = std::sqrt( pow(fx[i+1]-fx[i],2) + pow(fy[i+1]-fy[i],2) + pow(fz[i+1]-fz[i],2));
                scaleFactor = (170*nm)/length;
            }
            if ((i >= 5*32768+4096+1) && (i <= 171102)) {
                zshift = 8.5*um;
                length = std::sqrt( pow(fx[i+1]-fx[i],2) + pow(fy[i+1]-fy[i],2) + pow(fz[i+1]-fz[i],2));
                scaleFactor = (170*nm)/length;
            }

            // Get coordinates of midpoint fibre using current & prev coordinates
            G4double midpoint_x = (fx[i]+fx[i-1])/2*scaleFactor;
            G4double midpoint_y = (fy[i]+fy[i-1])/2*scaleFactor;
            G4double midpoint_z = (fz[i]+fz[i-1])/2*scaleFactor;

            // G4cout << "Midpoint--------------------------- " << G4endl;
            // G4cout << "midx " << midpoint_x*2/scaleFactor << G4endl;
            // G4cout << "midx " << midpoint_y*2/scaleFactor << G4endl;
            // G4cout << "midx " << midpoint_z*2/scaleFactor << G4endl;

            // Get vector representing the direction the current & prev coordinates
            G4ThreeVector midpoint = G4ThreeVector(midpoint_x,midpoint_y,midpoint_z+zshift);
            G4ThreeVector direction = G4ThreeVector(fx[i]-fx[i-1],fy[i]-fy[i-1],fz[i]-fz[i-1]);

            // Only one dimension changes along length of each fibre (stays in plane of other 2
            // dimensions). Next step should be along a previously unchanged dimension.
            // I
            G4RotationMatrix rotLoop = G4RotationMatrix();
            if ((direction.x() == 0) && (direction.y() == 0)) rotLoop.rotateZ(0*deg); // Was along z axis, no rotation
            if ((direction.y() == 0) && (direction.z() == 0)) rotLoop.rotateY(90*deg); // Was along x axis, rotate zround y
            if ((direction.x() == 0) && (direction.z() == 0)) rotLoop.rotateX(90*deg); // Was along y axis, rotate around x

            start = std::chrono::high_resolution_clock::now();  
            G4Transform3D transformLoop = G4Transform3D(rotLoop, midpoint); // Defines placement of new physical volume
            finish = std::chrono::high_resolution_clock::now();
            elapsed = finish-start;
            // G4cout << "Time to do G4Transform3D " << elapsed.count() << " s" << G4endl;     
            // G4cout << elapsed.count() << G4endl;   

            start = std::chrono::high_resolution_clock::now();   

            // pLoop[i] = CreatePhysicalVolume("Fibre",
            //                                 i,
            //                                 true,
            //                                 lFiber,
            //                                 &rotLoop,
            //                                 &midpoint,
            //                                 fEnvelopePhys);

            // pLoop[i] = CreatePhysicalVolume("Fibre",
            //                                 i,
            //                                 true,
            //                                 myLoop,
            //                                 &rotLoop,
            //                                 &midpoint,
            //                                 fEnvelopePhys);

         
            // if (i < 32768 && i > 32758) {
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lmLoop[0],          //logical volume
            //                                  "Chromoloop"+std::to_string(i),      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
            // }
            // else if (i == 32768) {
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lmLoop[1],          //logical volume
            //                                  "Chromoloop"+std::to_string(i),      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
            // }
            // else if (i > 32768 && i < 32778) {
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lmLoop[2],          //logical volume
            //                                  "Chromoloop"+std::to_string(i),      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
            // }

            // else if (i < 65536 && i > 65526) {
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lmLoop[3],          //logical volume
            //                                  "Chromoloop"+std::to_string(i),      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
            // }
            // else if (i == 65536) {
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lmLoop[4],          //logical volume
            //                                  "Chromoloop"+std::to_string(i),      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
            // }
            // else if (i > 65536 && i < 65546) {
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lmLoop[5],          //logical volume
            //                                  "Chromoloop"+std::to_string(i),      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
            // }


            // pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                              myLoop,          //logical volume
            //                              "Chromoloop"+std::to_string(i),      //name
            //                              fEnvelopeLog,      //mother
            //                              false,             //many
            //                              i,                 //copy number
            //                              false);            //overlap check

            if (i <= 13733){
                pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
                                             lLoop[i],          //logical volume
                                             "Chromoloop"+std::to_string(i),      //name
                                             fEnvelopeLog,      //mother
                                             false,             //many
                                             i,                 //copy number
                                             false);            //overlap check

                // THIS SHOULD BE MOVED OUTSIDE THIS LOOP - UNNECESSARY FOR EACH FIBRE                
                // G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
                // LoopVisAtt->SetForceSolid(true);
                // lmLoop[1]->SetVisAttributes(LoopVisAtt);
            }
            // else if ((i > 13733) && (i <= 27219)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[i],          //logical volume
            //                                  "Chromoloop"+std::to_string(i),      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     // G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.556863,0.137255,0.137255)); //firebrick
            //     // LoopVisAtt->SetForceSolid(true);
            //     // lLoop[2]->SetVisAttributes(LoopVisAtt);
                
            // }
            // else if ((i > 27219) && (i <= 38300)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[3],          //logical volume
            //                                  "Chromoloop"+std::to_string(i),      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(1.0,0.25,0.0)); //orangered
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[3]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 38300) && (i <= 48926)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[4],          //logical volume
            //                                  "Chromoloop"+std::to_string(i),      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0)); //blue
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[4]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 48926) && (i <= 58972)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[5],          //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.55,0.09,0.09)); //scarlet
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[5]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 58972) && (i <= 68466)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[6],          //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.89,0.47,0.20)); //mandarinorange
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[6]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 68466) && (i <= 77290)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[7],          //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.62352,0.372549,0.623529)); //blueviolet
                
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[7]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 77290) && (i <= 85416)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[8],          //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.6,0.8,0.196078)); //yellowgreen
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[8]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 85416) && (i <= 93219)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[9],          //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.498039)); //springgreen
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[9]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 93219) && (i <= 100739)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[10],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.73,0.16,0.96)); //purple
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[10]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 100739) && (i <= 108209)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[11],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.13,0.37,0.31)); //huntergreen
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[11]->SetVisAttributes(LoopVisAtt);
                
            // }
            // else if ((i > 108209) && (i <= 115558)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[12],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.74902,0.847059,0.847059)); //lightblue
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[12]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 115558) && (i <= 121899)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[13],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0)); //cyan
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[13]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 121899) && (i <= 127808)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[14],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.22,0.69,0.87)); //summersky
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[14]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 127808) && (i <= 133382)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[15],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(1.0,0.498039,0.0)); //coral
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[15]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 133382) && (i <= 138317)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[16],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.35,0.35,0.67)); //richblue
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[16]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 138317) && (i <= 142686)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[17],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.196078,0.8,0.196078)); //magenta
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[17]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 142686) && (i <= 146915)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[18],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.419608,0.137255,0.556863)); //darkslateblue
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[18]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 146915) && (i <= 150460)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[19],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0)); //yellow
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[19]->SetVisAttributes(LoopVisAtt);
                
            // }
            // else if ((i > 150460) && (i <= 153929)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[20],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.309804,0.184314,0.317647)); //violet
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[20]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 153929) && (i <= 156680)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[21],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.0,0.1,0.0)); //green
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[21]->SetVisAttributes(LoopVisAtt);
            // }
            // else if ((i > 156680) && (i <= 159288)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[22],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(1.0,0.11,0.68)); //spicypink
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[22]->SetVisAttributes(LoopVisAtt);
            // }
            // // I think this is the Y chromosome based on bp%; unsure why is assigned to same logical
            // // volume as chromosome 22
            // else if ((i > 159288) && (i <= 162496)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[22],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(1.0,0.68,0.68)); //???
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[22]->SetVisAttributes(LoopVisAtt);
            // }
            // // I think this is the X chromosome based on bp%
            // else if ((i > 162496)){
            //     pLoop[i] = new G4PVPlacement(transformLoop,     //translation + rotation
            //                                  lLoop[23],         //logical volume
            //                                  "Chromoloop",      //name
            //                                  fEnvelopeLog,      //mother
            //                                  false,             //many
            //                                  i,                 //copy number
            //                                  false);            //overlap check
                
            //     G4VisAttributes* LoopVisAtt = new G4VisAttributes(G4Colour(0.419608,0.137255,0.556863)); //slateblue
            //     LoopVisAtt->SetForceSolid(true);
            //     lLoop[23]->SetVisAttributes(LoopVisAtt);
            // }
            finish = std::chrono::high_resolution_clock::now();
            elapsed = finish-start;
            // G4cout << "Time to create G4PVPlacement & assign vis attributes " << elapsed.count() << " s" << G4endl;
            // G4cout << elapsed.count() << G4endl;   

            if (i == 3 || i == 100 || i == 1000 || i == 10000 || i == 20000 || i == 40000 || i == 100000 || i == 171000){ 
                G4cout << "--------------------------------" << G4endl;
                G4cout << i << G4endl;
                G4cout << elapsed.count() << G4endl;
            }
            

            /**************************************************************************/
            //         Subcomponent 6: DNA (histones with wrapped double helix)
            /**************************************************************************/
            
            // if (fBuildBases){

            //     // fGeoManager->Initialize();
            
            //     // start = std::chrono::high_resolution_clock::now();  
            //     // G4LogicalVolume* lFiber = fGeoManager->BuildLogicFiber(true);
            //     // finish = std::chrono::high_resolution_clock::now();
            //     // elapsed = finish-start;
            //     // G4cout << "Time to build logic fiber: " << elapsed.count() << " s" << G4endl;     
            
            //     start = std::chrono::high_resolution_clock::now();  
            //     G4VPhysicalVolume* pFiber = CreatePhysicalVolume("Fiber", lFiber, pLoop[i]);
            //     finish = std::chrono::high_resolution_clock::now();
            //     elapsed = finish-start;
            //     // G4cout << "Time to create pfiber physical volume: " << elapsed.count() << " s" << G4endl;     
            //     G4cout << elapsed.count() << G4endl;     

            // }
        
        }
        auto loopfinish = std::chrono::high_resolution_clock::now();
        elapsed = loopfinish-loopstart;
        // G4cout << "Time to build logic fiber: " << elapsed.count() << " s" << G4endl;  
        // G4cout << elapsed.count() << G4endl;       

    }
    
    // SetTooComplexForOGLS();

	InstantiateChildren(fEnvelopePhys);

    G4cout << "Finished in TsFractalDNAV2 " << G4endl;
	return fEnvelopePhys;
}


