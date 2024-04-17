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

#include "MiniOrangeActionInitializer.hh"
#include "MiniOrangeDetectorConstruction.hh"
#include "MiniOrangePrimaryGeneratorAction.hh"
#include "MiniOrangeRunAction.hh"
#include "MiniOrangeEventAction.hh"
#include "MiniOrangeTrackingAction.hh"
#include "MiniOrangeSteppingAction.hh"
#include "MiniOrangeSteppingVerbose.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MiniOrangeActionInitializer::MiniOrangeActionInitializer() : 
  G4VUserActionInitialization()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MiniOrangeActionInitializer::Build() const 
{
  const MiniOrangeDetectorConstruction* detector = 
        static_cast<const MiniOrangeDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction()); 

  SetUserAction(new MiniOrangePrimaryGeneratorAction());

  //Optional user classes
  SetUserAction(new MiniOrangeRunAction());
  SetUserAction(new MiniOrangeEventAction());
  SetUserAction(new MiniOrangeTrackingAction()); 
  SetUserAction(new MiniOrangeSteppingAction(detector));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MiniOrangeActionInitializer::BuildForMaster() const
{
  SetUserAction(new MiniOrangeRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VSteppingVerbose* MiniOrangeActionInitializer::InitializeSteppingVerbose() const
{
  // Verbose output class
  return new MiniOrangeSteppingVerbose();
}

