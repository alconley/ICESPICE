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
//    ********************************************
//    *                                          *
//    *    ICESPICEPrimaryGeneratorAction.cc     *
//    *                                          *
//    ********************************************
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "ICESPICEPrimaryGeneratorAction.hh"

#include "ICESPICEDetectorConstruction.hh"

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
#define RANDOM_DIRECTION 0

#define DIFFERENT_DIRECTION_SAME_ENERGY 0
#define SAME_DIRECTION_DIFFERENT_ENERGY 1


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICEPrimaryGeneratorAction::ICESPICEPrimaryGeneratorAction()
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

  //this function is called at the begining of event
  //Start position of primaries
  G4double z0 = 70.*mm;
  G4double x0 = 0.*mm;
  G4double y0 = 0.*mm;
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICEPrimaryGeneratorAction::~ICESPICEPrimaryGeneratorAction()
{
  delete particleGun;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ICESPICEPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    // Function to convert degrees to radians
  auto DegToRad = [](G4double angleInDegrees) {
    return angleInDegrees * CLHEP::pi / 180.0;
  };

  #if RANDOM_DIRECTION
    // Set a random direction within the lower half-sphere
    G4double phi = CLHEP::twopi * G4UniformRand();  // Random azimuthal angle between 0 and 2pi
    G4double theta = 0.5 * CLHEP::pi * G4UniformRand() + 0.5 * CLHEP::pi;  // Random polar angle between pi/2 and pi
    
    G4double ux = std::sin(theta) * std::cos(phi);
    G4double uy = std::sin(theta) * std::sin(phi);
    G4double uz = std::cos(theta);  // Positive values only, pointing in the positive z-direction
    
    particleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));
  #endif

  #if DIFFERENT_DIRECTION_SAME_ENERGY
    // Set a fixed energy of 1 MeV
    particleGun->SetParticleEnergy(1.*MeV);

    // Set the initial position of the particles to (0 mm, 0 mm, +70 mm)
    G4double x0 = 0.*mm;
    G4double y0 = 0.*mm;
    G4double z0 = 70.*mm;
    particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));

    // Define the angle range for the lower y-z plane in degrees
    const G4double minAngle = 20.0;
    const G4double maxAngle = 0.0;

    // Convert angle range to radians
    G4double minTheta = DegToRad(90.0 + minAngle);
    G4double maxTheta = DegToRad(90.0 + maxAngle);

    // Generate a random direction within the specified range on the lower y-z plane
    G4double theta = minTheta + (maxTheta - minTheta) * G4UniformRand();  // Random angle between min and max
    G4double phi = 0.0;  // No azimuthal variation for the y-z plane

    G4double uy = std::cos(theta);  // Direction component in y
    G4double uz = -std::sin(theta);  // Direction component in z (negative to point towards -z)

    particleGun->SetParticleMomentumDirection(G4ThreeVector(0., uy, uz));
  #endif

  #if SAME_DIRECTION_DIFFERENT_ENERGY
    // Set a fixed direction

    const G4double angle = 17.0;
    const G4double theta = DegToRad(90.0 + angle);

    G4double uy = std::cos(theta);  // Direction component in y
    G4double uz = -std::sin(theta);  // Direction component in z (negative to point towards -z)

    particleGun->SetParticleMomentumDirection(G4ThreeVector(0., uy, uz));

    // the energy can be 100 to 1500 keV in steps of 100 keV
    G4double energy = 100.0 * (G4UniformRand() * 14 + 1);

    particleGun->SetParticleEnergy(energy * keV);
  #endif

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


