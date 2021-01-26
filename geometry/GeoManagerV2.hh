//**************************************************************************************************
// This class is essentially a wrapper for the GeoCalculation and GeoVolume classes. Those classes
// do the actual work for creating and placing DNA volumes. This class allows public methods from
// those classes to be accessed under a shared interface.
//**************************************************************************************************

#ifndef GEOMANAGERV2_HH
#define GEOMANAGERV2_HH

#include "GeoCalculationV2.hh"
#include "GeoVolumeV2.hh"


class GeoManagerV2
{
public:
    GeoManagerV2(G4int verbose=0, G4double factor=1.);
    ~GeoManagerV2();

    //----------------------------------------------------------------------------------------------
    // Getters
    //----------------------------------------------------------------------------------------------
    G4int GetVerbose(){return fVerbose;}
    G4double GetFactor(){return fFactor;}

    //----------------------------------------------------------------------------------------------
    // Setters
    //----------------------------------------------------------------------------------------------
    void SetVerbose(G4int verbose){fVerbose = verbose;}
    void SetFactor(G4double factor){fFactor = factor;}

    //----------------------------------------------------------------------------------------------
    // Initialize parameters in the GeoCalculation object and GeoVolume object. GeoCalculation is
    // done first, and is subsequently used in the initialization of the GeoVolume object. Default
    // values for parameters provided in header file.
    //----------------------------------------------------------------------------------------------
    void Initialize( G4int BpNum=2, G4int NucleoNum=12);

    //----------------------------------------------------------------------------------------------
    // Wrapper for the BuildLogicFiber method of GeoVolume. 
    //----------------------------------------------------------------------------------------------
    G4LogicalVolume* BuildLogicFiber(G4bool cutVolumes=true, G4bool checkForOverlaps=true, 
        G4int overlapsResolution=100, G4bool quitIfOverlap=true);

    //----------------------------------------------------------------------------------------------
    // Wrapper for the GetDNAMoleculesPositions method of GeoVolume. 
    //----------------------------------------------------------------------------------------------
    std::map<G4String, std::vector<std::vector<double> > >* GetDNAMoleculesPositions();


    //void RemoveMoleculeFromChemStage(const G4String& molName, G4int copyNum, G4int strand)
    //{geoDnaMolecules->RemoveMoleculeFromChemStage(molName, copyNum, strand);}
    //void ResetRemoveMoleculeList(){geoDnaMolecules->ResetRemoveMoleculeList();}

private:
    G4int fVerbose;
    G4double fFactor;

    GeoCalculationV2* geoCalculation;
    GeoVolumeV2* geoVolume;



};

#endif // GEOMANAGERV2_HH
