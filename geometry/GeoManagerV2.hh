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

#ifndef GEOMANAGERV2_HH
#define GEOMANAGERV2_HH

#include "GeoCalculationV2.hh"
#include "GeoVolumeV2.hh"


class GeoManagerV2
{
public:
    GeoManagerV2(G4int verbose=0, G4double factor=1.);
    ~GeoManagerV2();

    void Initialize();

    G4LogicalVolume* BuildLogicFiber(G4bool isVisu=false);

    //void RemoveMoleculeFromChemStage(const G4String& molName, G4int copyNum, G4int strand)
    //{geoDnaMolecules->RemoveMoleculeFromChemStage(molName, copyNum, strand);}
    //void ResetRemoveMoleculeList(){geoDnaMolecules->ResetRemoveMoleculeList();}

    // Getters
    G4int GetVerbose(){return fVerbose;}
    G4double GetFactor(){return fFactor;}
    std::map<G4String, std::vector<std::vector<double> > >* GetDNAMoleculesPositions();

    // Setters
    void SetVerbose(G4int verbose){fVerbose = verbose;}
    void SetFactor(G4double factor){fFactor = factor;}

private:
    G4int fVerbose;
    G4double fFactor;

    GeoCalculationV2* geoCalculation;
    GeoVolumeV2* geoVolume;



};

#endif // GEOMANAGERV2_HH
