#include "ICESPICERunAction.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4AnalysisManager.hh"
#include "Randomize.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICERunAction::ICESPICERunAction()
  : G4UserRunAction()
  {   
    // set printing event number per each event
    // G4RunManager::GetRunManager()->SetPrintProgress(1);  

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;

    analysisManager->SetDefaultFileType("root");

    // Default settings
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFileName("ICESPICE");

    // analysisManager->SetNtupleMerging(true);

    analysisManager->CreateH1("Esil","Edep in silicon", 2000, 0., 2000.0*keV);

    // analysisManager->CreateNtuple("ICESPICE", "Edep");
    // analysisManager->CreateNtupleDColumn("Esil");
    // analysisManager->FinishNtuple();
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICERunAction::~ICESPICERunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ICESPICERunAction::BeginOfRunAction(const G4Run* aRun)
{  
    // Get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    analysisManager->Reset();

  // Open an output file 
  // it can be overwritten in a macro
    analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ICESPICERunAction::EndOfRunAction(const G4Run* aRun)
{      
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile(true);      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
