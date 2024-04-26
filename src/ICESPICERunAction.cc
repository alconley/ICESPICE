//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Code developed by:
//  S.Larsson and modified by Alex Conley
//
//    ******************************
//    *                            *
//    *    ICESPICERunAction.cc     *
//    *                            *
//    ******************************
//
//

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

    // Default settings
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFileName("ICESPICE.root");

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

  //inform the runManager to save random number seed
  // G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  
  if (IsMaster())    
    G4cout << "---> Run " << aRun->GetRunID() << " (master) start." 
	   << G4endl;
  else
    G4cout << "---> Run " << aRun->GetRunID() << " (worker) start." 
	   << G4endl;      

    // Get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file 
  // it can be overwritten in a macro
    analysisManager->OpenFile();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ICESPICERunAction::EndOfRunAction(const G4Run* aRun)
{      
  if (IsMaster())    
    G4cout << "Total number of event = " << aRun->GetNumberOfEvent() << G4endl;
  else
    G4cout << "Partial number of event in this worker = " 
	   << aRun->GetNumberOfEvent() << G4endl;
 
       
  if (IsMaster())
    {
      
    }

  // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(0) ) {      
    G4cout << " ESil : mean = " 
      << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy") 
      // get the counts greator than the first bin
      << " counts = " << analysisManager->GetH1(0)->entries()

      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;

  }

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
