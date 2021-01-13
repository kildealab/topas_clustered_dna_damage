//**************************************************************************************************
// Create a nuclear volume containing chromatin fibres arranged in a fractal pattern.
//**************************************************************************************************

#ifndef TsFractalDNAV2_hh
#define TsFractalDNAV2_hh

#include "TsVGeometryComponent.hh"
#include "GeoManagerV2.hh"


class TsFractalDNAV2 : public TsVGeometryComponent
{    
public:
	TsFractalDNAV2(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM,
                   TsGeometryManager* gM, TsVGeometryComponent* parentComponent, 
                   G4VPhysicalVolume* parentVolume, G4String& name);
	
    ~TsFractalDNAV2();
    
    GeoManagerV2* fGeoManager;
	
    //--------------------------------------------------------------------------------------------------
    // Handle parameters defined in the Topas parameter file
    // Read fiber coordinates from a data file
    //--------------------------------------------------------------------------------------------------
    void ResolveParameters();

    //----------------------------------------------------------------------------------------------
    // Mandatory method of TSVGeometryComponent to create solids, logical volumes, and physical
    // volumes. Primary Logical Volume, fEnvelopeLog, is defined as the nucleus, which may be
    // accessed in the  parameter file as Ge/MyDNA.
    //----------------------------------------------------------------------------------------------
	G4VPhysicalVolume* Construct();
    
    std::vector<G4ThreeVector> GetSugar1Info();
    
private:
    
    std::vector<G4VPhysicalVolume*> physFibers;
    G4VPhysicalVolume* physFiber;
    
    std::vector<G4double> fx;
    std::vector<G4double> fy;
    std::vector<G4double> fz;
    
    // std::vector<G4Tubs*> sLoop; // Vector of solids containing cylindrical fibers
    std::vector<G4LogicalVolume*> lLoop; // Vector of logical volumes containing chromosome territories
    std::vector<G4LogicalVolume*> lmLoop; // Vector of logical volumes containing chromosome territories
    
    G4bool fBuildBases;
    G4bool fGraphics;
    
    G4int fMaxNumFibers;
    G4int fNumFibers;  
};

#endif
