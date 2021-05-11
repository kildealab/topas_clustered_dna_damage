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

#ifndef G4EmDNAPhysics_option2and4_h
#define G4EmDNAPhysics_option2and4_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class TsParameterManager;

class G4EmDNAPhysics_option2and4 : public G4VPhysicsConstructor
{
public:
  G4EmDNAPhysics_option2and4(TsParameterManager* pM);
  G4EmDNAPhysics_option2and4(G4int ver = 0);

  virtual ~G4EmDNAPhysics_option2and4();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  G4int  verbose;
	
	TsParameterManager* fPm;
    G4bool fDeexcitationIgnoreCut;
};

#endif