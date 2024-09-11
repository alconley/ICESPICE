#include "ICESPICETrackingAction.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"

ICESPICETrackingAction::ICESPICETrackingAction() : G4UserTrackingAction() {}
ICESPICETrackingAction::~ICESPICETrackingAction() {}

void ICESPICETrackingAction::PreUserTrackingAction(const G4Track* track)
{
    // // Check if the particle is a secondary (ParentID > 0)
    // if (track->GetParentID() > 0) {
    //     G4String particleName = track->GetDefinition()->GetParticleName();
    //     G4double kineticEnergy = track->GetKineticEnergy() / CLHEP::keV;  // Kinetic energy in keV
    //     G4int parentID = track->GetParentID();
    //     G4int trackID = track->GetTrackID();
        
    //     // Get process that created this secondary particle
    //     G4String creationProcess = track->GetCreatorProcess()->GetProcessName();
        
    //     // Get particle's current position
    //     G4ThreeVector position = track->GetPosition();
    //     G4double posX = position.x() / CLHEP::mm;  // Position in mm
    //     G4double posY = position.y() / CLHEP::mm;
    //     G4double posZ = position.z() / CLHEP::mm;
        
    //     // Get momentum direction
    //     G4ThreeVector momentumDir = track->GetMomentumDirection();
        
    //     // Get the current volume to see if the particle is in the specific region
    //     G4VPhysicalVolume* currentVolume = track->GetVolume();
    //     G4String volumeName = currentVolume->GetName();

    //     // Check if the kinetic energy is in the desired range and the particle was created by Radioactivation
    //     if (kineticEnergy > 489 * CLHEP::keV && kineticEnergy < 490 * CLHEP::keV && creationProcess == "Radioactivation") {

    //     // if the volume is in the detector print
    //         // Print the information about the secondary particle
    //         G4cout << "Secondary Particle: " << particleName
    //             << " | Kinetic Energy: " << kineticEnergy << " keV"
    //             << " | Parent ID: " << parentID
    //             << " | Track ID: " << trackID
    //             << " | Created by process: " << creationProcess
    //             << " | Position (x, y, z): (" << posX << ", " << posY << ", " << posZ << ") mm"
    //             << " | Momentum Direction: (" << momentumDir.x() << ", " << momentumDir.y() << ", " << momentumDir.z() << ")"
    //             << " | Current Volume: " << volumeName
    //             << G4endl;
    //     }
    // }
}

void ICESPICETrackingAction::PostUserTrackingAction(const G4Track* track)
{
}








