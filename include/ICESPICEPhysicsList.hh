#ifndef ICESPICEPhysicsList_h
#define ICESPICEPhysicsList_h 1
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "G4EmConfigurator.hh"

class ICESPICEPhysicsList: public G4VModularPhysicsList {
  public:
    ICESPICEPhysicsList();
    ~ICESPICEPhysicsList() override = default;
    void SetCuts() override;
};
#endif



