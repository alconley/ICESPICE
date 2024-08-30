#include "ICESPICEPhysicsList.hh"
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
#include "G4RadioactiveDecayPhysics.hh"
#include "G4GenericIon.hh"

ICESPICEPhysicsList::ICESPICEPhysicsList():  G4VUserPhysicsList()
{
  RegisterPhysics(new G4DecayPhysics());
	RegisterPhysics(new G4EmStandardPhysics_option4());
  RegisterPhysics(new G4RadioactiveDecayPhysics());
}

ICESPICEPhysicsList::~ICESPICEPhysicsList()
{
}