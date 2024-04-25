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
//    ********************************
//    *                              *
//    *    ICESPICEEventAction.cc     *
//    *                              *
//    ********************************
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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
