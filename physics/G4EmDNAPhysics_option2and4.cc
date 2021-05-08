// Physics Module for g4em-dna_opt2and4
//
//**************************************************************************************************
// Author: Logan Montgomery
// Based on code written by:
//      - CM Lund et al. (2020). DOI:10.1016/j.ejmp.2020.04.001
//
// This class is a custom Topas physics module that combines the physics models available in both
// the G4EmDNAPhysics_option2 and G4EmDNAPhysics_option4 modules. Specifically: 
//    - the models included in option4 are used in the range: 10 eV to 10 keV
//    - the models included in option2 are used in the range: 10 keV to 1 MeV
//**************************************************************************************************

#include "G4EmDNAPhysics_option2and4.hh"

#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"

#include "G4DNAGenericIonsManager.hh"

// *** Processes and models for Geant4-DNA

#include "G4DNAElectronSolvation.hh"
#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"

#include "G4DNAIonisation.hh"
#include "G4DNABornIonisationModel.hh"
#include "G4DNAEmfietzoglouIonisationModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNABornExcitationModel.hh"
#include "G4DNAEmfietzoglouExcitationModel.hh"

#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"

#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"

// particles

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

// Warning : the following is needed in order to use EM Physics builders
// e+
#include "G4Positron.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
// gamma
#include "G4Gamma.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"

#include "G4EmParameters.hh"
// end of warning

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"

#include "G4EmConfigurator.hh"
#include "G4VEmModel.hh"
#include "G4DummyModel.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option2and4);


G4EmDNAPhysics_option2and4::G4EmDNAPhysics_option2and4(TsParameterManager* pM)
  : G4VPhysicsConstructor("G4EmStandard"), fPm(pM)
{
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();

    if (fPm->ParameterExists("Ph/Default/DeexcitationIgnoreCut"))
        fDeexcitationIgnoreCut = fPm->GetBooleanParameter("Ph/Default/DeexcitationIgnoreCut");
    else
        fDeexcitationIgnoreCut = false;

    G4cout << "fDeexcitationIgnoreCut = " << fDeexcitationIgnoreCut << G4endl;

    SetPhysicsType(bElectromagnetic);
}


G4EmDNAPhysics_option2and4::G4EmDNAPhysics_option2and4(G4int ver)
: G4VPhysicsConstructor("G4EmStandard"), verbose(ver)
{
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();
    param->SetFluo(true);  
    param->SetAuger(true);  
    param->SetAugerCascade(true);  
    param->SetDeexcitationIgnoreCut(true);
    fDeexcitationIgnoreCut = true;
    param->SetVerbose(verbose);
    SetPhysicsType(bElectromagnetic);
}


G4EmDNAPhysics_option2and4::~G4EmDNAPhysics_option2and4()
{}


void G4EmDNAPhysics_option2and4::ConstructParticle()
{
  // bosons
    G4Gamma::Gamma();

  // leptons
    G4Electron::Electron();
    G4Positron::Positron();
    
  // baryons
    G4Proton::Proton();

    G4GenericIon::GenericIonDefinition();

    G4DNAGenericIonsManager * genericIonsManager;
    genericIonsManager=G4DNAGenericIonsManager::Instance();
    genericIonsManager->GetIon("alpha++");
    genericIonsManager->GetIon("alpha+");
    genericIonsManager->GetIon("helium");
    genericIonsManager->GetIon("hydrogen");
}


void G4EmDNAPhysics_option2and4::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() )
  {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "e-") {
      
      // *** Solvation ***
      
      G4DNAElectronSolvation* solvation =
       new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
      G4DNAOneStepThermalizationModel* therm =
       new G4DNAOneStepThermalizationModel();
      therm->SetHighEnergyLimit(10.*eV);
      solvation->SetEmModel(therm);
      ph->RegisterProcess(solvation, particle);
      
      // *** Elastic scattering (two alternative models available) ***
      
      G4DNAElastic* theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");

      // To combine multiple models, start with a dummy model
      theDNAElasticProcess->SetEmModel(new G4DummyModel(), 1);

      ph->RegisterProcess(theDNAElasticProcess, particle);

      // *** Excitation ***
      G4DNAExcitation* theDNAExcitationProcess = new G4DNAExcitation("e-_G4DNAExcitation");

      // Only valid up to 10 keV
      theDNAExcitationProcess->SetEmModel(new G4DummyModel(), 1);
      ph->RegisterProcess(theDNAExcitationProcess, particle);

      // *** Ionisation ***
      G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("e-_G4DNAIonisation");

      // As above
      theDNAIonisationProcess->SetEmModel(new G4DNAEmfietzoglouIonisationModel());

      G4DNABornIonisationModel* higherEnergyIonisation = new G4DNABornIonisationModel(); //initiate Born Ionisation model
      higherEnergyIonisation->SetLowEnergyLimit(10.*keV); // set the min limit to 10 keV ( the maximum limit of Emfietzoglou)
      theDNAIonisationProcess->SetEmModel(higherEnergyIonisation);

      ph->RegisterProcess(theDNAIonisationProcess, particle);

      // *** Vibrational excitation ***
      ph->RegisterProcess(new G4DNAVibExcitation("e-_G4DNAVibExcitation"), particle);
      
      // *** Attachment ***
      ph->RegisterProcess(new G4DNAAttachment("e-_G4DNAAttachment"), particle); 
    
    } else if ( particleName == "proton" ) {
      ph->RegisterProcess(new G4DNAElastic("proton_G4DNAElastic"), particle);
      ph->RegisterProcess(new G4DNAExcitation("proton_G4DNAExcitation"), particle);
      ph->RegisterProcess(new G4DNAIonisation("proton_G4DNAIonisation"), particle);
      ph->RegisterProcess(new G4DNAChargeDecrease("proton_G4DNAChargeDecrease"), particle);

    } else if ( particleName == "hydrogen" ) {
      ph->RegisterProcess(new G4DNAElastic("hydrogen_G4DNAElastic"), particle);
      ph->RegisterProcess(new G4DNAExcitation("hydrogen_G4DNAExcitation"), particle);
      ph->RegisterProcess(new G4DNAIonisation("hydrogen_G4DNAIonisation"), particle);
      ph->RegisterProcess(new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease"), particle);

    } else if ( particleName == "alpha" ) {
      ph->RegisterProcess(new G4DNAElastic("alpha_G4DNAElastic"), particle);
      ph->RegisterProcess(new G4DNAExcitation("alpha_G4DNAExcitation"), particle);
      ph->RegisterProcess(new G4DNAIonisation("alpha_G4DNAIonisation"), particle);
      ph->RegisterProcess(new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease"), particle);

    } else if ( particleName == "alpha+" ) {
      ph->RegisterProcess(new G4DNAElastic("alpha+_G4DNAElastic"), particle);
      ph->RegisterProcess(new G4DNAExcitation("alpha+_G4DNAExcitation"), particle);
      ph->RegisterProcess(new G4DNAIonisation("alpha+_G4DNAIonisation"), particle);
      ph->RegisterProcess(new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease"), particle);
      ph->RegisterProcess(new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease"), particle);

    } else if ( particleName == "helium" ) {
      ph->RegisterProcess(new G4DNAElastic("helium_G4DNAElastic"), particle);
      ph->RegisterProcess(new G4DNAExcitation("helium_G4DNAExcitation"), particle);
      ph->RegisterProcess(new G4DNAIonisation("helium_G4DNAIonisation"), particle);
      ph->RegisterProcess(new G4DNAChargeIncrease("helium_G4DNAChargeIncrease"), particle);
    
    // Extension to HZE proposed by Z. Francis

    } else if ( particleName == "GenericIon" ) {
      ph->RegisterProcess(new G4DNAIonisation("GenericIon_G4DNAIonisation"), particle);
    }
    
    // Warning : the following particles and processes are needed by EM Physics builders
    // They are taken from the default Livermore Physics list
    // These particles are currently not handled by Geant4-DNA
    
      // e+
      
    else if (particleName == "e+") {

      // Identical to G4EmStandardPhysics_option3
      
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);      

      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);

    } else if (particleName == "gamma") {
    
      // photoelectric effect - Livermore model only
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
      ph->RegisterProcess(thePhotoElectricEffect, particle);

      // Compton scattering - Livermore model only
      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
      ph->RegisterProcess(theComptonScattering, particle);

      // gamma conversion - Livermore model below 80 GeV
      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel());
      ph->RegisterProcess(theGammaConversion, particle);

      // default Rayleigh scattering is Livermore
      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      ph->RegisterProcess(theRayleigh, particle);
    }
    
   // Warning : end of particles and processes are needed by EM Physics builders 
    
  }

  // Deexcitation
  //
  if (fDeexcitationIgnoreCut) {
      G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
      G4LossTableManager::Instance()->SetAtomDeexcitation(de);
  }


  G4EmConfigurator* em_config = 
      G4LossTableManager::Instance()->EmConfigurator();

  G4VEmModel* mod;

  // Use Champion Elastic above 10 keV
  mod = new G4DNAChampionElasticModel();
  mod->SetActivationLowEnergyLimit(10.*keV);
  em_config->SetExtraEmModel(
          "e-",              // Particle name
          "e-_G4DNAElastic", // Process name
          mod,               // G4VEmModel pointer
          "World",          // Region name
          10.*keV,           // Lower energy limit
          1.*MeV             // Upper energy limit
          );


  // Use Uehara below 10 keV
  mod = new G4DNAUeharaScreenedRutherfordElasticModel();
  em_config->SetExtraEmModel(
          "e-",              // Particle name
          "e-_G4DNAElastic", // Process name
          mod,               // G4VEmModel pointer
          "World",          // Region name
          0.,                // Lower energy limit
          10.*keV            // Upper energy limit
          );

  // Same for Born (above 10 keV), Emfietzoglou (below) for both excitation and ionization
  mod = new G4DNABornExcitationModel();
  mod->SetActivationLowEnergyLimit(10.*keV);
  em_config->SetExtraEmModel(
          "e-",
          "e-_G4DNAExcitation",
          mod,
          "World",
          10.*keV,
          1.*MeV
          );

  mod = new G4DNAEmfietzoglouExcitationModel();
  em_config->SetExtraEmModel(
          "e-",
          "e-_G4DNAExcitation",
          mod,
          "World",
          0.,
          10.*keV
          );

  mod = new G4DNABornIonisationModel();
  mod->SetActivationLowEnergyLimit(10.*keV);
  em_config->SetExtraEmModel(
          "e-",
          "e-_G4DNAIonisation",
          mod,
          "World",
          10.*keV,
          1.*MeV
          );

  mod = new G4DNAEmfietzoglouIonisationModel();
  em_config->SetExtraEmModel(
          "e-",
          "e-_G4DNAIonisation",
          mod,
          "World",
          0.,
          10.*keV
          );
}