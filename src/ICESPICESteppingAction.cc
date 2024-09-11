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
        const G4VProcess* creatorProcess = track->GetCreatorProcess();
        
        // For radioactive decay sources
        if (creatorProcess && creatorProcess->GetProcessName() == "Radioactivation") {
            fEventAction->AddSil(edep);
        }

        // Comment out if not using radioactive decay sources
        // fEventAction->AddSil(edep);

    }


    // // Check if the volume is the silicon detector
    // if (volume == fDetConstruction->GetSiliconPV()) {
    //     G4double edep = aStep->GetTotalEnergyDeposit();
        
    //     // Get the current track (i.e., the particle)
    //     G4Track* track = aStep->GetTrack();

    //     // Get the parent ID
    //     G4int parentID = track->GetParentID();

    //     // Check if the particle was created by radioactive decay
    //     const G4VProcess* creatorProcess = track->GetCreatorProcess();
    //     // G4cout << "Creator process: " << creatorProcess->GetProcessName() << G4endl;

    //     // make sure the particle is an electron
    //     // if (track->GetDefinition() == G4Electron::ElectronDefinition()) {
    //         if (edep > 0 ){
    //             if (creatorProcess->GetProcessName() == "Radioactivation") {
    //                 fEventAction->AddSil(edep);  // Fill histogram (or handle energy deposit) as per your action
    //             }
    //         }
    //     // }

    // }

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