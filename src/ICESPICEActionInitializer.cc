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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

