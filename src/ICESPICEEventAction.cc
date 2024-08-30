#include "ICESPICEEventAction.hh"
#include "ICESPICERunAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICEEventAction::ICESPICEEventAction()
  : G4UserEventAction(),
    fEnergySilicon(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICEEventAction::~ICESPICEEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ICESPICEEventAction::BeginOfEventAction(const G4Event* evt)
{ 
  fEnergySilicon = 0.;

  G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) 
   G4cout << "\n---> Begin Of Event: " << evtNb << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ICESPICEEventAction::EndOfEventAction(const G4Event* evt)
{  
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(0, fEnergySilicon);

  // analysisManager->FillNtupleDColumn(0, fEnergySilicon);
  // analysisManager->AddNtupleRow(); 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
