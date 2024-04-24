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
//    ********************************************
//    *                                          *
//    *    MiniOrangePrimaryGeneratorAction.cc     *
//    *                                          *
//    ********************************************
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "MiniOrangePrimaryGeneratorAction.hh"

#include "MiniOrangeDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"

#include "Randomize.hh" // Include this header for random number generation


//Print position of primaries.
#define POSITION 0

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MiniOrangePrimaryGeneratorAction::MiniOrangePrimaryGeneratorAction()
  :rndmVertex(false)
{
  //default kinematic
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  G4ParticleDefinition* particle
    = G4Electron::Definition();

  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleEnergy(1.*MeV);

  //Momentum Direction
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MiniOrangePrimaryGeneratorAction::~MiniOrangePrimaryGeneratorAction()
{
  delete particleGun;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MiniOrangePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  //Start position of primaries
  G4double z0 = 50.*mm;
  G4double x0 = 0.*mm;
  G4double y0 = 0.*mm;
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  // Set a random direction within the lower half-sphere
  G4double phi = CLHEP::twopi * G4UniformRand();  // Random azimuthal angle between 0 and 2pi
  G4double theta = 0.5 * CLHEP::pi * G4UniformRand() + 0.5 * CLHEP::pi;  // Random polar angle between pi/2 and pi
  
  G4double ux = std::sin(theta) * std::cos(phi);
  G4double uy = std::sin(theta) * std::sin(phi);
  G4double uz = std::cos(theta);  // Positive values only, pointing in the positive z-direction
  
  particleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));

  // particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));

  
  particleGun->GeneratePrimaryVertex(anEvent);

#if POSITION  
  G4cout << "\n----Particle Gun--------------------------------------------\n";
  G4cout << "\n ---> SetParticlePosition(G4ThreeVector(x0,y0,z0)) ";
  G4cout << "\n ---> x0 = " << x0 << " "<< G4BestUnit(x0,"Length");
  G4cout << "\n ---> y0 = " << y0 << " "<< G4BestUnit(y0,"Length");
  G4cout << "\n ---> z0 = " << z0 << " "<< G4BestUnit(z0,"Length");
  G4cout << "\n-----------------------------------------------------------\n";
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


