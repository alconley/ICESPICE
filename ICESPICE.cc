#include "G4Types.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4AnalysisManager.hh"
#include "ICESPICEDetectorConstruction.hh"
#include "ICESPICEPhysicsList.hh"
#include "ICESPICEActionInitializer.hh"

int main(int argc,char** argv) {

  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

 // Construct the default run manager
 //
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  G4int nThreads = 12;
  runManager->SetNumberOfThreads(nThreads);

  // set mandatory initialization classes
  runManager->SetUserInitialization(new ICESPICEDetectorConstruction);
  runManager->SetUserInitialization(new ICESPICEPhysicsList);
  runManager->SetUserInitialization(new ICESPICEActionInitializer());

  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  //Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc==1)   // Define UI session for interactive mode.
    {
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute vis.mac");
      if (ui->IsGUI()) {
        UImanager->ApplyCommand("/control/execute gui.mac");
      }
      ui->SessionStart();
      delete ui;
    }
  else           // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  // // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();

  // job termination
  delete visManager;

  delete runManager;

  return 0;
}

