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
//    *****************************************
//    *                                       *
//    *    ICESPICEDetectorConstruction.cc     *
//    *                                       *
//    *****************************************
//
//
#include "ICESPICEDetectorConstruction.hh"
#include "ICESPICETabulatedField3D.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EqMagElectricField.hh"

#include "G4ChordFinder.hh"
#include "G4UniformMagField.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4NistManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "G4Tubs.hh"

#include "G4GenericMessenger.hh"
#include "G4RunManager.hh"

#include "CADMesh.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Possibility to turn off (0) magnetic field and measurement volume. 
#define MAG 0          // Magnetic field grid
#define MAGNETS 0      // N42 1"X1"x1/8"
#define ATTENUATOR 0   // AC: Volume for attenuator 
#define DETECTOR 1     // AC: Volume for detector

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICEDetectorConstruction::ICESPICEDetectorConstruction()
  :physiWorld(NULL), logicWorld(NULL), solidWorld(NULL),
    physiDetector(NULL), logicDetector(NULL), // AC
    WorldMaterial(NULL), 
    AttenuatorMaterial(NULL), // AC
    MagnetMaterial(NULL), // AC
    DetectorMaterial(NULL), // AC
    DetectorHousingMaterial(NULL), // AC
    fMessenger(0)  
{
  fField.Put(0);
  WorldSizeXY=WorldSizeZ=0;
  DetectorPosition=-30.*mm; // AC
  DefineCommands();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICEDetectorConstruction::~ICESPICEDetectorConstruction()
{
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* ICESPICEDetectorConstruction::Construct()

{
  DefineMaterials();
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ICESPICEDetectorConstruction::DefineMaterials()
{ 
  //This function illustrates the possible ways to define materials.
  //Density and mass per mole taken from Physics Handbook for Science
  //and engineering, sixth edition. This is a general material list
  //with extra materials for other examples.
  
  G4String name, symbol;             
  G4double density;            
  
  G4int ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;
  
  // Define Elements  
  // Example: G4Element* Notation  = new G4Element ("Element", "Notation", z, a);
  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   N  = new G4Element ("Nitrogen", "N", 7., 14.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   Ar = new G4Element ("Argon" , "Ar", 18., 39.948*g/mole );

  // AC 
  G4Element*   Ta = new G4Element ("Tantalum" , "Ta", 73., 180.9479*g/mole ); // mini-orange attenuator
  density = 16.65*g/cm3 ;
  G4Material* Tantalum = new G4Material(name="Tantalum", density, ncomponents=1);
  Tantalum->AddElement(Ta, 1);
  
  // Define Elements
  G4Element* Nd = new G4Element("Neodymium", "Nd", 60., 144.24*g/mole);
  G4Element* Fe = new G4Element("Iron", "Fe", 26., 55.85*g/mole);
  G4Element* B = new G4Element("Boron", "B", 5., 10.81*g/mole);

  // Define N42 Magnet material (Neodymium Iron Boron)
  density = 7.5*g/cm3;
  G4Material* NdFeB = new G4Material(name="N42_Magnet", density, ncomponents=3);
  NdFeB->AddElement(Nd, natoms=2);
  NdFeB->AddElement(Fe, natoms=14);
  NdFeB->AddElement(B, natoms=1);

  G4NistManager* nist = G4NistManager::Instance();
  G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");
  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");

  // Define materials from elements.
  
  // Case 1: chemical molecule  
  // Water 
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="H2O"  , density, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  
  // Case 2: mixture by fractional mass.
  // Air
  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  // Vacuum
  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;         //from PhysicalConstants.h
  G4Material* vacuum = new G4Material(name="vacuum", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
  vacuum->AddMaterial(Air, fractionmass=1.);


  // Laboratory vacuum: Dry air (average composition)
  density = 1.7836*mg/cm3 ;       // STP
  G4Material* Argon = new G4Material(name="Argon", density, ncomponents=1);
  Argon->AddElement(Ar, 1);
  
  density = 1.25053*mg/cm3 ;       // STP
  G4Material* Nitrogen = new G4Material(name="N2", density, ncomponents=1);
  Nitrogen->AddElement(N, 2);
  
  density = 1.4289*mg/cm3 ;       // STP
  G4Material* Oxygen = new G4Material(name="O2", density, ncomponents=1);
  Oxygen->AddElement(O, 2);
  
  
  density  = 1.2928*mg/cm3 ;       // STP
  density *= 1.0e-8 ;              // pumped vacuum
  
  temperature = STP_Temperature;
  pressure = 1.0e-8*STP_Pressure;

  G4Material* LaboratoryVacuum = new G4Material(name="LaboratoryVacuum",
						density,ncomponents=3,
						kStateGas,temperature,pressure);
  LaboratoryVacuum->AddMaterial( Nitrogen, fractionmass = 0.7557 ) ;
  LaboratoryVacuum->AddMaterial( Oxygen,   fractionmass = 0.2315 ) ;
  LaboratoryVacuum->AddMaterial( Argon,    fractionmass = 0.0128 ) ;
  
  G4cout << G4endl << *(G4Material::GetMaterialTable()) << G4endl;

  // Default materials in setup.
  WorldMaterial = LaboratoryVacuum;
  AttenuatorMaterial = Tantalum; // AC
  MagnetMaterial = NdFeB; // AC
  DetectorMaterial = silicon; // AC
  DetectorHousingMaterial = aluminum; // AC

  G4cout << "end material"<< G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* ICESPICEDetectorConstruction::ConstructCalorimeter()
{  
  //The World
  WorldSizeXY  = 100.*mm;  // Cube
  WorldSizeZ   = 300.*mm;

  zOffset = 0.0*mm;  // Offset of the magnetic field grid

// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  solidWorld = new G4Box("World",				       //its name
			   WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2);  //its size
  
  logicWorld = new G4LogicalVolume(solidWorld,	        //its solid
				   WorldMaterial,	//its material
				   "World");		//its name
  
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  // Visualization attributes
  G4VisAttributes* simpleWorldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  simpleWorldVisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(simpleWorldVisAtt);
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#if DETECTOR

  PIPS1000Detector();

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#if ATTENUATOR

  auto attenuator = CADMesh::TessellatedMesh::FromPLY("./cad_files/tantalum_5_slot_attenuator.PLY");
  auto solidAttenuator = attenuator->GetSolid();
  auto logicAttenuator = new G4LogicalVolume(solidAttenuator, AttenuatorMaterial, "Attenuator");
  auto physiAttenuator = new G4PVPlacement(0,			             //no rotation
          G4ThreeVector(0.,0.,0.), //at (0,0,0)
                                "Attenuator",		             //its name
                                logicAttenuator,		             //its logical volume
                                physiWorld,			             //its mother  volume
                                false,			                     //no boolean operation
                                0);			                     //copy number

  // Visualization attributes
  G4VisAttributes* simpleAttenuatorVisAtt= new G4VisAttributes(G4Colour(0.25, 0.25, 0.25)); //grey
  simpleAttenuatorVisAtt->SetVisibility(true);
  simpleAttenuatorVisAtt->SetForceSolid(true);
  logicAttenuator->SetVisAttributes(simpleAttenuatorVisAtt);

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#if MAGNETS

  auto magnet = CADMesh::TessellatedMesh::FromPLY("./cad_files/1x1x1_8in_square_magnet.PLY");
  auto solidMagnet = magnet->GetSolid();
  auto logicMagnet = new G4LogicalVolume(solidMagnet, MagnetMaterial, "Magnet");

  // Calculate the placement and rotation for each magnet
  G4double placementRadius = 3.5*mm;  // Adjusting for the corner to be at 3.5mm
  G4int numMagnets = 5;
  G4double angleStep = 360.0*deg / numMagnets;

  for (int i = 0; i < numMagnets; i++) {
      G4double angle = i * angleStep;
      G4ThreeVector pos(placementRadius * std::sin(angle), placementRadius * std::cos(angle), 0);
      G4RotationMatrix* rot = new G4RotationMatrix();
      rot->rotateZ(angle); // Rotation to spread magnets around the origin

      new G4PVPlacement(rot,          // rotation
                        pos,          // position
                        logicMagnet,  // logical volume
                        "Magnet",     // name
                        logicWorld,   // mother volume
                        false,        // no boolean operations
                        i);           // copy number
  }

  // Set visualization attributes to make it look shiny
  G4VisAttributes* MagnetVisAtt = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));  // light grey color
  MagnetVisAtt->SetVisibility(true);
  MagnetVisAtt->SetForceSolid(true);
  logicMagnet->SetVisAttributes(MagnetVisAtt);
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  return physiWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ICESPICEDetectorConstruction::ConstructSDandField()
{
//  Magnetic Field
//
#if MAG
  
  if (fField.Get() == 0)
    {
      //Field grid in A9.TABLE. File must be in accessible from run urn directory. 
      G4MagneticField* ICESPICEField= new ICESPICETabulatedField3D("ICESPICE3D.TABLE", zOffset);
      fField.Put(ICESPICEField);
      
      //This is thread-local
      G4FieldManager* pFieldMgr = 
	G4TransportationManager::GetTransportationManager()->GetFieldManager();
           
      G4cout<< "DeltaStep "<<pFieldMgr->GetDeltaOneStep()/mm <<"mm" <<G4endl;
      //G4ChordFinder *pChordFinder = new G4ChordFinder(ICESPICEField);

      pFieldMgr->SetDetectorField(fField.Get());
      pFieldMgr->CreateChordFinder(fField.Get());
      
    }
#endif
}

void ICESPICEDetectorConstruction::DefineCommands() {

    fMessenger = new G4GenericMessenger(this, 
                                        "/ICESPICE/Detector/", 
                                        "Detector control");

    // Detector Position
    G4GenericMessenger::Command& detectorPosition
      = fMessenger->DeclareMethodWithUnit("Position","mm",
                                  &ICESPICEDetectorConstruction::SetDetectorPosition, 
                                  "Set the position of the detector");
    detectorPosition.SetParameterName("position", true);
    detectorPosition.SetRange("position>-100. && position<=0.");
    detectorPosition.SetDefaultValue("-30.");

    G4GenericMessenger::Command& detectorType
      = fMessenger->DeclareMethod("SetDetector",
                                  &ICESPICEDetectorConstruction::SetDetector, 
                                  "Set the detector type");
    detectorType.SetParameterName("type", true);
    detectorType.SetDefaultValue("PIPS1000");
}

void ICESPICEDetectorConstruction::SetDetectorPosition(G4double val) {
    DetectorPosition = val;
    physiDetector->SetTranslation(G4ThreeVector(0, 0, DetectorPosition));
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void ICESPICEDetectorConstruction::SetDetector(G4String name) {

  if (name == "PIPS100") {
    PIPS100Detector();
  } else if (name == "PIPS300") {
    PIPS300Detector();
  } else if (name == "PIPS500") {
    PIPS500Detector();
  } else if (name == "PIPS1000") {
    PIPS1000Detector();
  } else {
    G4cerr << "Unknown detector type: " << name << G4endl;
  }

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void ICESPICEDetectorConstruction::PIPS100Detector() {
    // Assuming that the detector window and housing are positioned relative to the detector's dimensions.

    auto detector = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips100/active_area.PLY");

    auto solidDetector = detector->GetSolid();
    logicDetector = new G4LogicalVolume(solidDetector,
                                          DetectorMaterial,
                                          "DetectorWindow");

    G4double DetectorActiveArea = 50.*mm2; // Active area of the detector
    G4double DetectorWindowThickness = 50.*nanometer; // Thickness of the detector window
    G4double DetectorRadius = std::sqrt(DetectorActiveArea / 3.14);

    G4Tubs* solidDetectorWindow = new G4Tubs("DetectorWindow",
                                      0,  // Inner radius
                                      std::sqrt(DetectorActiveArea / 3.14),  // Outer radius
                                      DetectorWindowThickness / 2.,  // Half-height
                                      0.*deg,  // Start angle
                                      360.*deg);  // Spanning angle

    G4LogicalVolume* logicDetectorWindow = new G4LogicalVolume(solidDetectorWindow,
                                        DetectorMaterial,
                                        "DetectorWindow");

    // Recalculate the position if it's dependent on the detector's thickness
    G4double windowZPosition = - DetectorWindowThickness / 2.;

    new G4PVPlacement(nullptr,  // No rotation
                G4ThreeVector(0, 0, windowZPosition),  // Position in the detector
                logicDetectorWindow,
                "DetectorWindow",
                logicDetector,  // Parent volume
                false,  // No boolean operation
                0);  // Copy number

    // Create the outer housing for the detector

    auto detectorHousing = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips100/detector_housing.PLY");
    auto solidDetectorHousing = detectorHousing->GetSolid();
    auto logicDetectorHousing = new G4LogicalVolume(solidDetectorHousing,
                                              DetectorHousingMaterial,
                                              "DetectorHousing");
    // Place the detector within the housing
    new G4PVPlacement(nullptr,  // No rotation
                      G4ThreeVector(0, 0, 0),  // Position relative to housing center
                      logicDetectorHousing,
                      "Detector",
                      logicDetector,  // Parent volume
                      false,  // No boolean operation
                      0);  // Copy number
                
    physiDetector = new G4PVPlacement(nullptr,  // no rotation
            G4ThreeVector(0, 0, DetectorPosition),  // position in world
            logicDetector,  // logical volume to place
            "Detector",  // name
            logicWorld,  // parent volume (world)
            false,  // no boolean operation
            0);  // copy number


    // Visualization attributes for various components
      G4VisAttributes* visAttributesDetector = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));  // Green for the detector
      visAttributesDetector->SetVisibility(true);
      visAttributesDetector->SetForceSolid(true);
      logicDetector->SetVisAttributes(visAttributesDetector);

      G4VisAttributes* visAttributesWindow = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));  // Red for the window
      visAttributesWindow->SetVisibility(true);
      visAttributesWindow->SetForceSolid(true);
      logicDetectorWindow->SetVisAttributes(visAttributesWindow);

      G4VisAttributes* visAttributesHousing = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));  // Gray for the housing
      visAttributesHousing->SetVisibility(true);
      visAttributesHousing->SetForceSolid(false);
      logicDetectorHousing->SetVisAttributes(visAttributesHousing);

    }

void ICESPICEDetectorConstruction::PIPS300Detector() {
    // Assuming that the detector window and housing are positioned relative to the detector's dimensions.

    auto detector = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips300/active_area.PLY");

    auto solidDetector = detector->GetSolid();
    logicDetector = new G4LogicalVolume(solidDetector,
                                          DetectorMaterial,
                                          "DetectorWindow");

    G4double DetectorActiveArea = 50.*mm2; // Active area of the detector
    G4double DetectorWindowThickness = 50.*nanometer; // Thickness of the detector window
    G4double DetectorRadius = std::sqrt(DetectorActiveArea / 3.14);

    G4Tubs* solidDetectorWindow = new G4Tubs("DetectorWindow",
                                      0,  // Inner radius
                                      std::sqrt(DetectorActiveArea / 3.14),  // Outer radius
                                      DetectorWindowThickness / 2.,  // Half-height
                                      0.*deg,  // Start angle
                                      360.*deg);  // Spanning angle

    G4LogicalVolume* logicDetectorWindow = new G4LogicalVolume(solidDetectorWindow,
                                        DetectorMaterial,
                                        "DetectorWindow");

    // Recalculate the position if it's dependent on the detector's thickness
    G4double windowZPosition = - DetectorWindowThickness / 2.;

    new G4PVPlacement(nullptr,  // No rotation
                G4ThreeVector(0, 0, windowZPosition),  // Position in the detector
                logicDetectorWindow,
                "DetectorWindow",
                logicDetector,  // Parent volume
                false,  // No boolean operation
                0);  // Copy number

    // Create the outer housing for the detector

    auto detectorHousing = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips300/detector_housing.PLY");
    auto solidDetectorHousing = detectorHousing->GetSolid();
    auto logicDetectorHousing = new G4LogicalVolume(solidDetectorHousing,
                                              DetectorHousingMaterial,
                                              "DetectorHousing");
    // Place the detector within the housing
    new G4PVPlacement(nullptr,  // No rotation
                      G4ThreeVector(0, 0, 0),  // Position relative to housing center
                      logicDetectorHousing,
                      "Detector",
                      logicDetector,  // Parent volume
                      false,  // No boolean operation
                      0);  // Copy number
                
      physiDetector = new G4PVPlacement(nullptr,  // no rotation
                  G4ThreeVector(0, 0, DetectorPosition),  // position in world
                  logicDetector,  // logical volume to place
                  "Detector",  // name
                  logicWorld,  // parent volume (world)
                  false,  // no boolean operation
                  0);  // copy number

    // Visualization attributes for various components
      G4VisAttributes* visAttributesDetector = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));  // Green for the detector
      visAttributesDetector->SetVisibility(true);
      visAttributesDetector->SetForceSolid(true);
      logicDetector->SetVisAttributes(visAttributesDetector);

      G4VisAttributes* visAttributesWindow = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));  // Red for the window
      visAttributesWindow->SetVisibility(true);
      visAttributesWindow->SetForceSolid(true);
      logicDetectorWindow->SetVisAttributes(visAttributesWindow);

      G4VisAttributes* visAttributesHousing = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));  // Gray for the housing
      visAttributesHousing->SetVisibility(true);
      visAttributesHousing->SetForceSolid(false);
      logicDetectorHousing->SetVisAttributes(visAttributesHousing);

    }

void ICESPICEDetectorConstruction::PIPS500Detector() {
    // Assuming that the detector window and housing are positioned relative to the detector's dimensions.

    auto detector = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips500/active_area.PLY");

    auto solidDetector = detector->GetSolid();
    logicDetector = new G4LogicalVolume(solidDetector,
                                          DetectorMaterial,
                                          "DetectorWindow");

    G4double DetectorActiveArea = 50.*mm2; // Active area of the detector
    G4double DetectorWindowThickness = 50.*nanometer; // Thickness of the detector window
    G4double DetectorRadius = std::sqrt(DetectorActiveArea / 3.14);

    G4Tubs* solidDetectorWindow = new G4Tubs("DetectorWindow",
                                      0,  // Inner radius
                                      std::sqrt(DetectorActiveArea / 3.14),  // Outer radius
                                      DetectorWindowThickness / 2.,  // Half-height
                                      0.*deg,  // Start angle
                                      360.*deg);  // Spanning angle

    G4LogicalVolume* logicDetectorWindow = new G4LogicalVolume(solidDetectorWindow,
                                        DetectorMaterial,
                                        "DetectorWindow");

    // Recalculate the position if it's dependent on the detector's thickness
    G4double windowZPosition = - DetectorWindowThickness / 2.;

    new G4PVPlacement(nullptr,  // No rotation
                G4ThreeVector(0, 0, windowZPosition),  // Position in the detector
                logicDetectorWindow,
                "DetectorWindow",
                logicDetector,  // Parent volume
                false,  // No boolean operation
                0);  // Copy number

    // Create the outer housing for the detector

    auto detectorHousing = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips500/detector_housing.PLY");
    auto solidDetectorHousing = detectorHousing->GetSolid();
    auto logicDetectorHousing = new G4LogicalVolume(solidDetectorHousing,
                                              DetectorHousingMaterial,
                                              "DetectorHousing");
    // Place the detector within the housing
    new G4PVPlacement(nullptr,  // No rotation
                      G4ThreeVector(0, 0, 0),  // Position relative to housing center
                      logicDetectorHousing,
                      "Detector",
                      logicDetector,  // Parent volume
                      false,  // No boolean operation
                      0);  // Copy number
                
      physiDetector = new G4PVPlacement(nullptr,  // no rotation
                  G4ThreeVector(0, 0, DetectorPosition),  // position in world
                  logicDetector,  // logical volume to place
                  "Detector",  // name
                  logicWorld,  // parent volume (world)
                  false,  // no boolean operation
                  0);  // copy number

    // Visualization attributes for various components
      G4VisAttributes* visAttributesDetector = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));  // Green for the detector
      visAttributesDetector->SetVisibility(true);
      visAttributesDetector->SetForceSolid(true);
      logicDetector->SetVisAttributes(visAttributesDetector);

      G4VisAttributes* visAttributesWindow = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));  // Red for the window
      visAttributesWindow->SetVisibility(true);
      visAttributesWindow->SetForceSolid(true);
      logicDetectorWindow->SetVisAttributes(visAttributesWindow);

      G4VisAttributes* visAttributesHousing = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));  // Gray for the housing
      visAttributesHousing->SetVisibility(true);
      visAttributesHousing->SetForceSolid(false);
      logicDetectorHousing->SetVisAttributes(visAttributesHousing);

    }

void ICESPICEDetectorConstruction::PIPS1000Detector() {
    // Assuming that the detector window and housing are positioned relative to the detector's dimensions.

    auto detector = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips1000/active_area.PLY");

    auto solidDetector = detector->GetSolid();
    logicDetector = new G4LogicalVolume(solidDetector,
                                          DetectorMaterial,
                                          "DetectorWindow");

    G4double DetectorActiveArea = 50.*mm2; // Active area of the detector
    G4double DetectorWindowThickness = 50.*nanometer; // Thickness of the detector window
    G4double DetectorRadius = std::sqrt(DetectorActiveArea / 3.14);

    G4Tubs* solidDetectorWindow = new G4Tubs("DetectorWindow",
                                      0,  // Inner radius
                                      std::sqrt(DetectorActiveArea / 3.14),  // Outer radius
                                      DetectorWindowThickness / 2.,  // Half-height
                                      0.*deg,  // Start angle
                                      360.*deg);  // Spanning angle

    G4LogicalVolume* logicDetectorWindow = new G4LogicalVolume(solidDetectorWindow,
                                        DetectorMaterial,
                                        "DetectorWindow");

    // Recalculate the position if it's dependent on the detector's thickness
    G4double windowZPosition = - DetectorWindowThickness / 2.;

    new G4PVPlacement(nullptr,  // No rotation
                G4ThreeVector(0, 0, windowZPosition),  // Position in the detector
                logicDetectorWindow,
                "DetectorWindow",
                logicDetector,  // Parent volume
                false,  // No boolean operation
                0);  // Copy number

    // Create the outer housing for the detector

    auto detectorHousing = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips1000/detector_housing.PLY");
    auto solidDetectorHousing = detectorHousing->GetSolid();
    auto logicDetectorHousing = new G4LogicalVolume(solidDetectorHousing,
                                              DetectorHousingMaterial,
                                              "DetectorHousing");
    // Place the detector within the housing
    new G4PVPlacement(nullptr,  // No rotation
                      G4ThreeVector(0, 0, 0),  // Position relative to housing center
                      logicDetectorHousing,
                      "Detector",
                      logicDetector,  // Parent volume
                      false,  // No boolean operation
                      0);  // Copy number
                
    physiDetector = new G4PVPlacement(nullptr,  // no rotation
                G4ThreeVector(0, 0, DetectorPosition),  // position in world
                logicDetector,  // logical volume to place
                "Detector",  // name
                logicWorld,  // parent volume (world)
                false,  // no boolean operation
                0);  // copy number

    // Visualization attributes for various components
    G4VisAttributes* visAttributesDetector = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));  // Green for the detector
    visAttributesDetector->SetVisibility(true);
    visAttributesDetector->SetForceSolid(true);
    logicDetector->SetVisAttributes(visAttributesDetector);

    G4VisAttributes* visAttributesWindow = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));  // Red for the window
    visAttributesWindow->SetVisibility(true);
    visAttributesWindow->SetForceSolid(true);
    logicDetectorWindow->SetVisAttributes(visAttributesWindow);

    G4VisAttributes* visAttributesHousing = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));  // Gray for the housing
    visAttributesHousing->SetVisibility(true);
    visAttributesHousing->SetForceSolid(false);
    logicDetectorHousing->SetVisAttributes(visAttributesHousing);

    }
