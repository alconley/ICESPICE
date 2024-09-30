#include "ICESPICESteppingAction.hh"
#include "ICESPICERunAction.hh"
#include "ICESPICEEventAction.hh"
#include "ICESPICEDetectorConstruction.hh"

#include "G4SteppingManager.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4AnalysisManager.hh"
#include "G4Track.hh" 
#include "G4SystemOfUnits.hh"

#define STOPPARTICLES 0

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICESteppingAction::ICESPICESteppingAction(
                      const ICESPICEDetectorConstruction* detectorConstruction,
                      ICESPICEEventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICESteppingAction::~ICESPICESteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ICESPICESteppingAction::UserSteppingAction(const G4Step* aStep)
  
{ 
    // get volume of the current step
	auto volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

	// // energy deposit
	auto edep = aStep->GetTotalEnergyDeposit();
			
    if (volume == fDetConstruction->GetSiliconPV()) {
        G4Track* track = aStep->GetTrack();

        // Print track information
        G4int trackID = track->GetTrackID();
        G4int parentID = track->GetParentID();
        G4String particleName = track->GetDefinition()->GetParticleName();
        G4double kineticEnergy = track->GetKineticEnergy();
        G4double globalTime = track->GetGlobalTime();
        G4double localTime = track->GetLocalTime();
        G4double trackLength = track->GetTrackLength();
        G4ThreeVector position = track->GetPosition();
        G4ThreeVector momentumDirection = track->GetMomentumDirection();
        G4double velocity = track->GetVelocity();
        G4VPhysicalVolume* currentVolume = track->GetVolume();
        G4int stepNumber = track->GetCurrentStepNumber();
        const G4ParticleDefinition* particleDefinition = track->GetDefinition();

        const G4VProcess* creatorProcess = track->GetCreatorProcess();
        
        // For radioactive decay sources
        // if (creatorProcess->GetProcessName() == "Radioactivation") {
            // make it only electrons
            // fEventAction->AddSil(edep);
            
            // if ( kineticEnergy > 489 * keV && kineticEnergy < 490 * keV 
            // || kineticEnergy > 984 * keV && kineticEnergy < 985 * keV           
            //    ) {
            //     // Print details
            //     G4cout << "Track ID: " << trackID << G4endl;
            //     G4cout << "Parent ID: " << parentID << G4endl;
            //     G4cout << "Particle: " << particleName << G4endl;
            //     G4cout << "Kinetic Energy: " << kineticEnergy / keV << " MeV" << G4endl;
            //     // G4cout << "Global Time: " << globalTime / ns << " ns" << G4endl;
            //     // G4cout << "Local Time: " << localTime / ns << " ns" << G4endl;
            //     // G4cout << "Track Length: " << trackLength / mm << " mm" << G4endl;
            //     // G4cout << "Position: " << position << G4endl;
            //     // G4cout << "Momentum Direction: " << momentumDirection << G4endl;
            //     // G4cout << "Velocity: " << velocity / (m/s) << " m/s" << G4endl;
            //     if (currentVolume) {
            //         G4String volumeName = currentVolume->GetName();
            //         G4cout << "Current Volume: " << volumeName << G4endl;
            //     }
            //     // G4cout << "Step Number: " << stepNumber << G4endl;
            //     // G4cout << "Particle Charge: " << particleDefinition->GetPDGCharge() << G4endl;
            //     // G4cout << "Particle Mass: " << particleDefinition->GetPDGMass() / MeV << " MeV" << G4endl;

            // }


        // }



        // Comment out if not using radioactive decay sources
        fEventAction->AddSil(edep);

    }

    #if STOPPARTICLES

        if (volume == fDetConstruction->GetAttenuatorPV() ) {
            G4Track* track = aStep->GetTrack();
            track->SetTrackStatus(fStopAndKill);
        }

        if ( volume == fDetConstruction->GetDetectorHousingPV() ) {
            G4Track* track = aStep->GetTrack();
            track->SetTrackStatus(fStopAndKill);
        }

        if ( volume == fDetConstruction->GetDetectorWindowPV() ) {
            G4Track* track = aStep->GetTrack();
            track->SetTrackStatus(fStopAndKill);
        }

    #endif
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....