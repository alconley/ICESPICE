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
//    ********************************
//    *                              *
//    *    MiniOrangePhysicsList.cc     *
//    *                              *
//    ********************************
//

#include "MiniOrangePhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"              
#include "G4EmStandardPhysics_option4.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4DecayPhysics.hh"

MiniOrangePhysicsList::MiniOrangePhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1*micrometer;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  cutForProton    = defaultCutValue;

  fEmPhysicsList = new G4EmStandardPhysics_option4();
  fDecPhysicsList = new G4DecayPhysics();
  SetVerboseLevel(1);
}

MiniOrangePhysicsList::~MiniOrangePhysicsList()
{
 delete fDecPhysicsList;
 delete fEmPhysicsList;
}

void MiniOrangePhysicsList::ConstructParticle()
{
 fDecPhysicsList -> ConstructParticle();
} 

void MiniOrangePhysicsList::ConstructProcess()
{
  AddTransportation();
  fEmPhysicsList -> ConstructProcess();
  
// Deexcitation
// Both Fluorescence and Auger e- emission activated
//G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
//G4LossTableManager::Instance()->SetAtomDeexcitation(de);
//de -> SetFluo(true);
//de -> SetAuger(true);
}

void MiniOrangePhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "MiniOrangePhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  
  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton 
  SetCutValue(cutForProton, "proton");
  SetCutValue(cutForProton, "anti_proton");
  
  //  SetCutValueForOthers(defaultCutValue);
  
  if (verboseLevel>0) DumpCutValuesTable();
}

void MiniOrangePhysicsList::SetGammaLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "MiniOrangePhysicsList::SetCuts:";
    G4cout << "Gamma cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  // G4Gamma::SetEnergyRange(lowcut,1e5); 
  SetGELowLimit(lowcut);
}

void MiniOrangePhysicsList::SetElectronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    
    G4cout << "MiniOrangePhysicsList::SetCuts:";
    G4cout << "Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  // G4Electron::SetEnergyRange(lowcut,1e5);
  SetGELowLimit(lowcut);
}

void MiniOrangePhysicsList::SetPositronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    
    G4cout << "MiniOrangePhysicsList::SetCuts:";
    G4cout << "Positron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  G4cerr << "MiniOrangePhysicsList::SetPositronLowLimit: Not currently able to set Positron LowLimit." << G4endl;
  G4Exception("MiniOrangePhysicsList::SetPositronLowLimit()","PurMag001",
	      FatalException,"Positron Low Limit: not implemented in MiniOrangePhysicsList"); 
  //
  // G4Positron::SetEnergyRange(lowcut,1e5);
}


void MiniOrangePhysicsList::SetProtonLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    
    G4cout << "MiniOrangePhysicsList::SetCuts:";
    G4cout << "Proton cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;  
  }  

  G4cerr << "MiniOrangePhysicsList::SetProtonLowLimit: Not currently able to set Proton LowLimit." << G4endl;
  G4Exception("MiniOrangePhysicsList::SetProtonLowLimit()","PurMag002",
	      FatalException,"Proton Low Limit: not implemented in MiniOrangePhysicsList"); 
  //
  // G4Proton::SetEnergyRange(lowcut,1e5);
  // G4AntiProton::SetEnergyRange(lowcut,1e5);
}

void MiniOrangePhysicsList::SetGEPLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "MiniOrangePhysicsList::SetGEPLowLimit:";
    G4cout << "Gamma and Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  // G4Gamma::SetEnergyRange(lowcut,1e5);
  // G4Electron::SetEnergyRange(lowcut,1e5);
  // G4Positron::SetEnergyRange(lowcut,1e5);
  this->SetGELowLimit(lowcut); 

  G4cerr << " SetGEPLowLimit : Uncertain whether setting Positron low limit " << G4endl;
}

void MiniOrangePhysicsList::SetGELowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "MiniOrangePhysicsList::SetGELowLimit:";
    G4cout << "Gamma and Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  
 
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowcut,1e5);
}
void MiniOrangePhysicsList::SetGammaCut(G4double val)
{
  cutForGamma = val;
}

void MiniOrangePhysicsList::SetElectronCut(G4double val)
{
  cutForElectron = val;
}

void MiniOrangePhysicsList::SetPositronCut(G4double val)
{
  cutForPositron = val;
}

void MiniOrangePhysicsList::SetProtonCut(G4double val)
{
  cutForProton = val;
}






