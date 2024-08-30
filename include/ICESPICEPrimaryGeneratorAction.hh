#ifndef ICESPICEPrimaryGeneratorAction_h
#define ICESPICEPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ICESPICEPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  ICESPICEPrimaryGeneratorAction();    
  ~ICESPICEPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event*);
  
private:
  G4GeneralParticleSource* particleSource;
};

#endif


