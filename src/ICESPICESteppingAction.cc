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
//    ***********************************
//    *                                 *
//    *    ICESPICESteppingAction.cc     *
//    *                                 *
//    ***********************************
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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
