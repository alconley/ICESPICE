#include "ICESPICEActionInitializer.hh"
#include "ICESPICEDetectorConstruction.hh"
#include "ICESPICEPrimaryGeneratorAction.hh"
#include "ICESPICERunAction.hh"
#include "ICESPICEEventAction.hh"
#include "ICESPICETrackingAction.hh"
#include "ICESPICESteppingAction.hh"
#include "ICESPICESteppingVerbose.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ICESPICEActionInitializer::ICESPICEActionInitializer() : 
  G4VUserActionInitialization()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ICESPICEActionInitializer::Build() const 
{
  const ICESPICEDetectorConstruction* detector = 
        static_cast<const ICESPICEDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction()); 

  SetUserAction(new ICESPICEPrimaryGeneratorAction());

  //Optional user classes
  SetUserAction(new ICESPICERunAction());
  SetUserAction(new ICESPICEEventAction());
  SetUserAction(new ICESPICETrackingAction()); 

  auto eventAction = new ICESPICEEventAction;
  SetUserAction(eventAction);
  SetUserAction(new ICESPICESteppingAction(detector,eventAction));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ICESPICEActionInitializer::BuildForMaster() const
{
  SetUserAction(new ICESPICERunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VSteppingVerbose* ICESPICEActionInitializer::InitializeSteppingVerbose() const
{
  // Verbose output class
  return new ICESPICESteppingVerbose();
}

