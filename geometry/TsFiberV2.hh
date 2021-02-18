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


class TsFiberV2 : public TsVGeometryComponent
{    
public:
	TsFiberV2(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsFiberV2();
    
    GeoManagerV2* fGeoManager;
	G4VPhysicalVolume* Construct();
    
    //void SetMaterial (G4String materialChoice);
    
    void UpdateGeometry();
    
private:
    
    G4Material *OtherMaterial(G4String materialName);
    G4double fWrapperRadius;
    G4double fWrapperHeight;
    G4bool fCutVolumes;
    G4bool fCheckForOverlaps;
    G4int fOverlapsResolution;
    G4bool fQuitIfOverlap;
    G4Material* fDNAMaterial;

    void ThrowOverlapError();
};

#endif
