#include "ICESPICESteppingAction.hh"
#include "ICESPICERunAction.hh"
#include "ICESPICEEventAction.hh"
#include "ICESPICEDetectorConstruction.hh"

#include "G4SteppingManager.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4AnalysisManager.hh"

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
	
	// energy deposit
	auto edep = aStep->GetTotalEnergyDeposit();
			
	if ( volume == fDetConstruction->GetSiliconPV() ) {
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
