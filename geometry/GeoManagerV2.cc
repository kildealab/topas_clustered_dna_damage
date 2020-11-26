// Extra Class for TsFiber
//
//**************************************************************************************************
// This class is essentially a wrapper for the GeoCalculation and GeoVolume classes. Those classes
// do the actual work for creating and placing DNA volumes. This class allows public methods from
// those classes to be accessed under a shared interface.
//**************************************************************************************************

#include "GeoManagerV2.hh"

//--------------------------------------------------------------------------------------------------
// Constructor: Creates a GeoManager object, thereby instantiating the class members.
//--------------------------------------------------------------------------------------------------
GeoManagerV2::GeoManagerV2(G4int verbose, G4double factor)
{
    fVerbose = verbose;
    fFactor = factor;

    geoCalculation = new GeoCalculationV2(fVerbose, fFactor);
    geoVolume = new GeoVolumeV2(fVerbose, fFactor);

}

//--------------------------------------------------------------------------------------------------
// Destructor
//--------------------------------------------------------------------------------------------------
GeoManagerV2::~GeoManagerV2()
{
    // clean the memory
    delete geoCalculation;
    delete geoVolume;
}

//--------------------------------------------------------------------------------------------------
// Initialize parameters in the GeoCalculation object and GeoVolume object. GeoCalculation is done
// first, and is subsequently used in the initialization of the GeoVolume object.
// Default values for parameters provided in header file.
//--------------------------------------------------------------------------------------------------
void GeoManagerV2::Initialize(G4int BpNum, G4int NucleoNum)
{
    // do the calculations using GeoCalculation
    geoCalculation->Initialize();

    // Set the parameters of the GeoVolume class.
    // This must be done before calling any of the Build method in the class.
    geoVolume->SetBaseRadius(geoCalculation->GetBaseRadius());
    geoVolume->SetBaseRadiusWater(geoCalculation->GetBaseRadiusWater());
    geoVolume->SetSugarTMPRadius(geoCalculation->GetSugarTMPRadius());
    geoVolume->SetSugarTHFRadius(geoCalculation->GetSugarTHFRadius());
    geoVolume->SetSugarTMPRadiusWater(geoCalculation->GetSugarTMPRadiusWater());
    geoVolume->SetSugarTHFRadiusWater(geoCalculation->GetSugarTHFRadiusWater());

    geoVolume->SetHistoneHeight(geoCalculation->GetHistoneHeight());
    geoVolume->SetHistoneRadius(geoCalculation->GetHistoneRadius());

    geoVolume->SetNucleoNum(NucleoNum);
    geoVolume->SetBpNum(BpNum);

    geoVolume->SetFiberPitch(geoCalculation->GetFiberPitch());
    geoVolume->SetFiberDeltaAngle(geoCalculation->GetFiberDeltaAngle());
    geoVolume->SetFiberNbNuclPerTurn(geoCalculation->GetFiberNbNuclPerTurn());
}

//--------------------------------------------------------------------------------------------------
// Wrapper for the BuildLogicFiber method of GeoVolume. 
//--------------------------------------------------------------------------------------------------
G4LogicalVolume* GeoManagerV2::BuildLogicFiber(G4bool isVisu)
{
    G4LogicalVolume* lFiber = geoVolume->BuildLogicFiber(geoCalculation->GetAllDNAVolumePositions(),
                                                             geoCalculation->GetNucleosomePosition(),
                                                             geoCalculation->GetPosAndRadiusMap(),
                                                             isVisu);

    return lFiber;
}

//--------------------------------------------------------------------------------------------------
// Wrapper for the GetDNAMoleculesPositions method of GeoVolume. 
//--------------------------------------------------------------------------------------------------
std::map<G4String, std::vector<std::vector<double> > > *GeoManagerV2::GetDNAMoleculesPositions()
{
    return geoVolume->GetDNAMoleculesPositions();
}
