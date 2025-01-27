#include "ICESPICEPhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"
#include "G4ios.hh"       

#include "G4EmStandardPhysics.hh"       
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmStandardPhysicsSS.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4DecayPhysics.hh"

#include "G4RadioactiveDecayPhysics.hh"

#include "G4GenericIon.hh"

#include "G4NuclideTable.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"

ICESPICEPhysicsList::ICESPICEPhysicsList() 
: G4VModularPhysicsList() {
  
  // protected member of G4VModularPhysicsList
  verboseLevel = 0;

  fEmPhysicsList = new G4EmStandardPhysics_option4(verboseLevel);

  // Deacy physics and all particles
  fDecPhysicsList = new G4DecayPhysics(verboseLevel);

  // mandatory for G4NuclideTable
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(1.0*picosecond);
  G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);
   
  //read new PhotonEvaporation data set 
  G4DeexPrecoParameters* deex = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  deex->SetCorrelatedGamma(false);
  deex->SetStoreAllLevels(true);
  deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()
                /std::log(2.));

  G4LossTableManager::Instance();
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(0.01*eV, 1*GeV);
  SetDefaultCutValue(1*micrometer);

  G4EmParameters* params = G4EmParameters::Instance();
  params->SetDeexcitationIgnoreCut(false);
  params->SetApplyCuts(true);
  params->Dump();
}

ICESPICEPhysicsList::~ICESPICEPhysicsList() {
  delete fEmPhysicsList;
  delete fDecPhysicsList;
}

void ICESPICEPhysicsList::ConstructParticle()
{
  fDecPhysicsList->ConstructParticle();
}

void ICESPICEPhysicsList::ConstructProcess()
{
  AddTransportation();

  fEmPhysicsList->ConstructProcess();

  fDecPhysicsList->ConstructProcess();

  AddRadioactiveDecay();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4PhysicsListHelper.hh"
#include "G4RadioactiveDecay.hh"
#include "G4Radioactivation.hh"
#include "G4GenericIon.hh"

void ICESPICEPhysicsList::AddRadioactiveDecay()
{  
  G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();
  
  // G4Radioactivation* radioactiveDecay = new G4Radioactivation();

  radioactiveDecay->SetARM(true);        //Atomic Rearangement
  radioactiveDecay->SetThresholdForVeryLongDecayTime(1e+60);

  // Initialize atomic deexcitation
  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* deex = man->AtomDeexcitation();
  if (!deex) {
     G4EmParameters::Instance()->SetFluo(true);
     G4EmParameters::Instance()->SetAugerCascade(true);
     deex = new G4UAtomicDeexcitation();
     deex->InitialiseAtomicDeexcitation();
     deex->SetAuger(false);
     deex->SetFluo(false);
     deex->SetPIXE(false);
     deex->SetVerboseLevel(1);
     man->SetAtomDeexcitation(deex);
  } else {
      G4EmParameters::Instance()->SetFluo(false);
      deex = new G4UAtomicDeexcitation();
      deex->InitialiseAtomicDeexcitation();
      deex->SetFluo(false);
      deex->SetVerboseLevel(1);
      man->SetAtomDeexcitation(deex);
  }

  // Register radioactive decay for relevant particles
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
  ph->RegisterProcess(radioactiveDecay, G4Electron::Electron());
  ph->RegisterProcess(radioactiveDecay, G4Gamma::Gamma());
}