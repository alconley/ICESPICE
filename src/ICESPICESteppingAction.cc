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
    // // get volume of the current step
	// auto volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

	// // energy deposit
	auto edep = aStep->GetTotalEnergyDeposit();

    // // Get the secondaries produced in this step
    // // const G4TrackVector* secondaries = aStep->GetSecondary();

    // // Check if there are any secondaries
    // // if (secondaries->size() > 0) {
    // //     G4cout << "Secondaries produced in this step: " << G4endl;
        
    // //     // Loop through the secondaries and print details
    // //     for (size_t i = 0; i < secondaries->size(); ++i) {
    // //         G4Track* secondary = (*secondaries)[i];
    // //         G4ParticleDefinition* particleDef = secondary->GetDefinition();
    // //         G4String particleName = particleDef->GetParticleName();
    // //         G4double energy = secondary->GetKineticEnergy();
    // //         G4cout << "  Particle: " << particleName 
    // //             << "  Energy: " << energy / keV << " keV"
    // //             << G4endl;
    // //     }
    // // }
			
    // if (volume == fDetConstruction->GetSiliconPV()) {
    //     fEventAction->AddSil(edep);
    // }

    auto prePV = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (prePV->GetLogicalVolume() == fDetConstruction->GetSiliconLV()) {
        fEventAction->AddSil(aStep->GetTotalEnergyDeposit());
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