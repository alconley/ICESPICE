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
//  S.Larsson
//
//    ******************************
//    *                            *
//    *    MiniOrangeRunAction.cc     *
//    *                            *
//    ******************************
//
//

#include "MiniOrangeRunAction.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4AnalysisManager.hh"
#include "Randomize.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MiniOrangeRunAction::MiniOrangeRunAction()
  : G4UserRunAction()
  {   
    // set printing event number per each event
    G4RunManager::GetRunManager()->SetPrintProgress(1);  

    auto analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;

    analysisManager->SetVerboseLevel(1);
    analysisManager->SetNtupleMerging(true);

    analysisManager->CreateH1("Esil","Edep in silicon", 1000, 0., 2.0*MeV);
    analysisManager->CreateH1("Lsil","trackL in silicon", 1000, 0., 1*mm);

    analysisManager->CreateNtuple("MiniOrange", "Edep and TrackL");
    analysisManager->CreateNtupleDColumn("Esil");
    analysisManager->CreateNtupleDColumn("Lsil");
    analysisManager->FinishNtuple();
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MiniOrangeRunAction::~MiniOrangeRunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MiniOrangeRunAction::BeginOfRunAction(const G4Run* aRun)
{  
  //Analysis must be handled by master and workers
  Book();  
  
  if (IsMaster())    
    G4cout << "---> Run " << aRun->GetRunID() << " (master) start." 
	   << G4endl;
  else
    G4cout << "---> Run " << aRun->GetRunID() << " (worker) start." 
	   << G4endl;      

    // Get analysis manager
    auto analysisManager = G4AnalysisManager::Instance();

    // Open an output file
    G4String fileName = "MiniOrange";
    analysisManager->OpenFile(fileName);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MiniOrangeRunAction::EndOfRunAction(const G4Run* aRun)
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
    auto analysisManager = G4AnalysisManager::Instance();
    if ( analysisManager->GetH1(1) ) {
      G4cout << G4endl << " ----> print histograms statistic ";
      if(isMaster) {
        G4cout << "for the entire run " << G4endl << G4endl; 
      }
      else {
        G4cout << "for the local thread " << G4endl << G4endl; 
      }
      
      G4cout << " ESil : mean = " 
        << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy") 
        << " rms = " 
        << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;
      G4cout << " LSil : mean = " 
        << G4BestUnit(analysisManager->GetH1(1)->mean(), "Length") 
        << " rms = " 
        << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Length") << G4endl;
    }

    // save histograms & ntuple
    //
    analysisManager->Write();
    analysisManager->CloseFile();            
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MiniOrangeRunAction::Book() 
{
  // Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->SetDefaultFileType("csv");
 
  // Open an output file
  man->OpenFile("MiniOrange"); 
  man->SetFirstNtupleId(1);
  
  // Create 1st ntuple (id = 1)
  //    
  man->CreateNtuple("n101", "Electron");
  man->CreateNtupleDColumn("ex");
  man->CreateNtupleDColumn("ey");
  man->CreateNtupleDColumn("ez");
  man->CreateNtupleDColumn("ee");
  man->CreateNtupleDColumn("epx");
  man->CreateNtupleDColumn("epy");
  man->CreateNtupleDColumn("epz");
  man->FinishNtuple();
  G4cout << "Ntuple-1 created" << G4endl;

  // Create 2nd ntuple (id = 2)
  // //    
  // man->CreateNtuple("n102", "Gamma");
  // man->CreateNtupleDColumn("gx");
  // man->CreateNtupleDColumn("gy");
  // man->CreateNtupleDColumn("gz");
  // man->CreateNtupleDColumn("ge");
  // man->CreateNtupleDColumn("gpx");
  // man->CreateNtupleDColumn("gpy");
  // man->CreateNtupleDColumn("gpz");
  // man->FinishNtuple();
  // G4cout << "Ntuple-2 created" << G4endl;
 
  // // Create 3rd ntuple (id = 3)
  // //
  // man->CreateNtuple("n103", "Positron");
  // man->CreateNtupleDColumn("px");
  // man->CreateNtupleDColumn("py");
  // man->CreateNtupleDColumn("pz");
  // man->CreateNtupleDColumn("pe");
  // man->CreateNtupleDColumn("ppx");
  // man->CreateNtupleDColumn("ppy");
  // man->CreateNtupleDColumn("ppz");
  // man->FinishNtuple();
  // G4cout << "Ntuple-3 created" << G4endl;

  return;
}



