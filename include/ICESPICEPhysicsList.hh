#ifndef ICESPICEPhysicsList_h
#define ICESPICEPhysicsList_h 1
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "G4EmConfigurator.hh"

class ICESPICEPhysicsList: public G4VModularPhysicsList {
  public:
    ICESPICEPhysicsList();
    ~ICESPICEPhysicsList();

    virtual void ConstructParticle();
    void AddPhysicsList(const G4String& name);
    virtual void ConstructProcess();

    void AddRadioactiveDecay();

  private:
    G4String                             fEmName;
    G4VPhysicsConstructor*               fEmPhysicsList;
    G4VPhysicsConstructor*               fDecPhysicsList;

};
#endif



