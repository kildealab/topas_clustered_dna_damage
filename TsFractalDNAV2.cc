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
// Constructor: Creates a GeoManager object
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
// Read fiber coordinates from a data file
//--------------------------------------------------------------------------------------------------
void TsFractalDNAV2::ResolveParameters(){
    // User defined parameters.
    // fPm is a TsParamater object, which is a member variable of TsVGeometry Component
    
    // Initialize various member variables
    fMaxNumFibers = fPm->GetIntegerParameter(GetFullParmName("MaxNumFibers"));

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
        G4cout << "Topas is building bases." << G4endl;
        G4cout << "This can slow down simulations." << G4endl;
    }
    
    //Hilbert space filling:
    const char* filename = FileName;
    std::string line = "";
    G4double x, y, z;
    
    ifstream f(filename, ios::in);
    
    // Read through file containing xyz coordinates, save to G4Double vectors fx, fy, and fz
    if (f.is_open())
    {
        std::string line;
        int i_line = 1;
        while ( getline (f,line) )
        {
            if (line.empty()) {
                break;
            }

            std::istringstream line_stream(line);
            std::string stoken; // store individual values between delimiters on a line

            // Process 3 coordinates per row
            int tokens_per_row = 3;
            int i_token = 0;
            while (getline(line_stream, stoken, ' ')) {
                i_token += 1;
                if (i_token == 1) 
                    fx.push_back(stod(stoken));
                else if (i_token == 2) 
                    fy.push_back(stod(stoken));
                else {
                    fz.push_back(stod(stoken));
                    i_token = 0;
                }
            }

            if (i_token != 0) {
                G4cerr << "Error: there is an incorrect number of coordinates on line " << i_line << " of " <<FileName << G4endl;
                exit(1);
            }
            i_line += 1;
        }
    }
    else
    {
        G4cout << "ERROR: Unable to open file " << FileName << G4endl;
        exit(1);
    }

    fNumFibers = fx.size()-1; // first fiber is defined by first 2 coordinates (i.e. one less fiber than # coordinates)
    physFibers.resize(fNumFibers); // allocate enough space for physical volumes of fibers

    f.close();
}


//--------------------------------------------------------------------------------------------------
// Mandatory method of TSVGeometryComponent to create solids, logical volumes, and physical volumes.
// Primary Logical Volume, fEnvelopeLog, is defined as the nucleus, which may be accessed in the 
// parameter file as Ge/MyDNA.
//--------------------------------------------------------------------------------------------------
G4VPhysicalVolume* TsFractalDNAV2::Construct()
{
	BeginConstruction();

    /**************************************************************************/
    //                 Wrapper: Nucleus
    /**************************************************************************/
    
    G4Ellipsoid* fEnvelopeSolid = new G4Ellipsoid("nucl", 6*um, 4*um, 11*um); // Parameterize this
    fEnvelopeLog= CreateLogicalVolume(fEnvelopeSolid);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

    /**************************************************************************/
    //                 Subcomponent: chromatin fibers
    /**************************************************************************/
    
    // Define fibre dimensions (scale xyz coordinates by desired fibre length)
    G4double fiberRadius = fPm->GetDoubleParameter(GetFullParmName("FiberRadius"), "Length"); // 12 nm
    G4double length = std::sqrt( pow(fx[2]-fx[1],2) + pow(fy[2]-fy[1],2) + pow(fz[2]-fz[1],2));
    G4double scaleFactor = (fPm->GetDoubleParameter(GetFullParmName("FiberLength"), "Length"))/length;  // 170 nm
    length=scaleFactor*length;

    // Create G4 solid: Only one solid fibre is created
    G4String subComponentName1 = "SolidFibre";
    G4Tubs* solidFiber = new G4Tubs(subComponentName1, 0, fiberRadius, (length/2)-fiberRadius, 0*deg, 360*deg);

    // Create single logical volume for the solid fiber
    G4LogicalVolume* logicalFiber = CreateLogicalVolume("Fibers",solidFiber);

    //----------------------------------------------------------------------------------------------
    // Build bases here (once)
    //----------------------------------------------------------------------------------------------
    // if (fBuildBases){
    auto start = std::chrono::high_resolution_clock::now();
    // fGeoManager->Initialize();
    auto finish = std::chrono::high_resolution_clock::now();
    std:chrono::duration<double> elapsed = finish-start;
    // G4cout << "Time to Initialize GeoManager: " << elapsed.count() << " s" << G4endl;  
    // start = std::chrono::high_resolution_clock::now();  
    // G4LogicalVolume* lFiber = fGeoManager->BuildLogicFiber(true);
    // finish = std::chrono::high_resolution_clock::now();
    // elapsed = finish-start;
    // G4cout << "Time to build logic fiber: " << elapsed.count() << " s" << G4endl;  
    // G4cout << elapsed.count() << G4endl;       
    // }


    //----------------------------------------------------------------------------------------------
    // Fill the ellipsoid nucleus with the fibers, according to the read-in fx, fy, fz data.
    //----------------------------------------------------------------------------------------------
    G4int max_num_iterations = std::min(fMaxNumFibers,fNumFibers);
    G4double xshift = 0, yshift = 0, zshift = 0;
    for (G4int i = 1; i <= max_num_iterations; i++) {

        auto loopstart = std::chrono::high_resolution_clock::now();

        //------------------------------------------------------------------------------------------
        // There are discontinuities at the 6 fibre indices labeled below. Unsure why, perhaps the
        // fractal space filling used to generate the .dat file was unable to generate the 
        // longitudinal extent required for this work? The pattern in FullGenome.dat repeats at
        // line 32769, 65537, etc. (i.e. resets to x=0.48438, y=0.48438, z=0.48438) but new zshift
        // is applied as below
        // It seems like the fibres that would exist at the following indices are not generated at all.
        //------------------------------------------------------------------------------------------
        if ((i != 32768) || (i != 2*32768) || (i != 3*32768) || (i != 4*32768) || (i != 5*32768) || (i != 5*32768+4096)) {
            
            //--------------------------------------------------------------------------------------
            // Apply z shifts at specific intervals along placement of the fibers.
            //--------------------------------------------------------------------------------------
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

            //--------------------------------------------------------------------------------------
            // Get coordinates of midpoint fibre using current & prev coordinates
            //--------------------------------------------------------------------------------------
            G4double midpoint_x = (fx[i]+fx[i-1])/2*scaleFactor;
            G4double midpoint_y = (fy[i]+fy[i-1])/2*scaleFactor;
            G4double midpoint_z = (fz[i]+fz[i-1])/2*scaleFactor;
            G4ThreeVector midpoint = G4ThreeVector(midpoint_x,midpoint_y,midpoint_z+zshift);

            //--------------------------------------------------------------------------------------
            // Establish the rotation for the current fiber by evaluating the direction between
            // the current data point and previous data point.
            //--------------------------------------------------------------------------------------
            G4ThreeVector direction = G4ThreeVector(fx[i]-fx[i-1],fy[i]-fy[i-1],fz[i]-fz[i-1]);

            G4RotationMatrix rotLoop = G4RotationMatrix();
            if ((direction.x() == 0) && (direction.y() == 0)) rotLoop.rotateZ(0*deg); // If should be placed along z, no rotation needed
            if ((direction.y() == 0) && (direction.z() == 0)) rotLoop.rotateY(90*deg); // If should be placed along x, rotate about y
            if ((direction.x() == 0) && (direction.z() == 0)) rotLoop.rotateX(90*deg); // If should be placed along y, rotate about z

            //--------------------------------------------------------------------------------------
            // Create a physical volume to place the new fiber, using the translation (midpoint) &
            // rotation information
            //--------------------------------------------------------------------------------------
            // start = std::chrono::high_resolution_clock::now();  
            G4Transform3D transformLoop = G4Transform3D(rotLoop, midpoint); // Defines placement of new physical volume
            // finish = std::chrono::high_resolution_clock::now();
            // elapsed = finish-start;
            // G4cout << "Time to do G4Transform3D " << elapsed.count() << " s" << G4endl;     
            // G4cout << elapsed.count() << G4endl;   

            start = std::chrono::high_resolution_clock::now();   

            G4String phys_name = "fibre_" + std::to_string(i);
            physFibers[i] = CreatePhysicalVolume(phys_name,
                                            i,
                                            true,
                                            logicalFiber,
                                            transformLoop,
                                            fEnvelopePhys);
         
            finish = std::chrono::high_resolution_clock::now();
            elapsed = finish-start;

            // G4cout << "Time to create G4PVPlacement & assign vis attributes " << elapsed.count() << " s" << G4endl;
            // G4cout << elapsed.count() << G4endl;   
            if (i == 3 || i == 100 || i == 1000 || i == 10000 || i == 20000 || i == 40000 || i == 60000 || i == 80000 || i == 100000 || i == 171000){ 
                G4cout << "--------------------------------" << G4endl;
                G4cout << i << G4endl;
                G4cout << elapsed.count() << G4endl;
            }
            

            //--------------------------------------------------------------------------------------
            // If need to build bases within each loop. Shouldn't be necessary --> do only once,
            // outside of loop
            //--------------------------------------------------------------------------------------
            
            // if (fBuildBases){

            //     // fGeoManager->Initialize();
            
            //     // start = std::chrono::high_resolution_clock::now();  
            //     // G4LogicalVolume* lFiber = fGeoManager->BuildLogicFiber(true);
            //     // finish = std::chrono::high_resolution_clock::now();
            //     // elapsed = finish-start;
            //     // G4cout << "Time to build logic fiber: " << elapsed.count() << " s" << G4endl;     
            
            //     start = std::chrono::high_resolution_clock::now();  
            //     G4VPhysicalVolume* pFiber = CreatePhysicalVolume("Fiber", lFiber, physFibers[i]);
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