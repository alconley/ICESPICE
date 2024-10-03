#include "ICESPICEDetectorConstruction.hh"
#include "ICESPICETabulatedField3D.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EqMagElectricField.hh"
#include "G4ChordFinder.hh"
#include "G4UniformMagField.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4GenericMessenger.hh"
#include "G4RunManager.hh"
#include "CADMesh.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Possibility to turn off (0) magnetic field and measurement volume. 
#define MAG 0         // Magnetic field grid
#define MAGNETS 0      // N42 1"X1"x1/8"
#define ATTENUATOR 0   // AC: Volume for attenuator 
#define DETECTOR 1     // AC: Volume for detector
#define DETECTORHOLDER 1 // AC: Volume for detector holder
#define MAGNETHOLDER 0 // AC: Volume for magnet holder/mounting rings
#define SOURCEBACKING 1 // AC: Volume for source backing

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICEDetectorConstruction::ICESPICEDetectorConstruction()
  :physiWorld(NULL), logicWorld(NULL), solidWorld(NULL),
    physiDetector(NULL), logicDetector(NULL), solidDetector(NULL),
    physiDetectorWindow(NULL), logicDetectorWindow(NULL), solidDetectorWindow(NULL),
    physiDetectorHousing(NULL), logicDetectorHousing(NULL), solidDetectorHousing(NULL),
    physiDetectorHolder(NULL), logicDetectorHolder(NULL), solidDetectorHolder(NULL),
    physiAttenuator(NULL), logicAttenuator(NULL), solidAttenuator(NULL),
    physiSourceBacking(NULL), logicSourceBacking(NULL), solidSourceBacking(NULL),
    WorldMaterial(NULL), 
    AttenuatorMaterial(NULL), 
    MagnetMaterial(NULL), 
    DetectorMaterial(NULL), 
    DetectorHousingMaterial(NULL), 
    DetectorWindowMaterial(NULL),
    SourceBackingMaterial(NULL),
    fMessenger(0)  
{
  fField.Put(0);
  WorldSizeXY=WorldSizeZ=0;
  DetectorPosition=-30.*mm; // AC
  SourcePosition=70.*mm; // AC
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

  G4NistManager* nist = G4NistManager::Instance();

  // Define Elements  
  G4Element* H = nist->FindOrBuildElement("H");
  G4Element* N = nist->FindOrBuildElement("N");
  G4Element* O = nist->FindOrBuildElement("O");
  G4Element* Ar = nist->FindOrBuildElement("Ar");
  G4Element* Ta = nist->FindOrBuildElement("Ta");
  G4Element* Fe = nist->FindOrBuildElement("Fe");
  G4Element* Cr = nist->FindOrBuildElement("Cr");
  G4Element* Ni = nist->FindOrBuildElement("Ni");
  G4Element* B = nist->FindOrBuildElement("B");
  G4Element* Nd = nist->FindOrBuildElement("Nd");
  G4Element* Si = nist->FindOrBuildElement("Si");
  G4Element* C = nist->FindOrBuildElement("C");
  
  // Define Tantalum for the attenuator
  G4Material* Tantalum = new G4Material(name="Tantalum", 16.65*g/cm3, ncomponents=1);
  Tantalum->AddElement(Ta, 1);
  
  // Define N42 Magnet material (Neodymium Iron Boron) for magnets
  G4Material* NdFeB = new G4Material(name="N42_Magnet", 7.5*g/cm3, ncomponents=3);
  NdFeB->AddElement(Nd, natoms=2);
  NdFeB->AddElement(Fe, natoms=14);
  NdFeB->AddElement(B, natoms=1);

  G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");
  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");

  // create SiO2 material for detector window (maybe) otherwise its silicon
  G4Material* SiO2 = new G4Material("SiO2", 2.65*g/cm3, 2);
  SiO2->AddElement(Si, 1);
  SiO2->AddElement(O, 2);

  // Define stainless steel for detector housing
  G4Material* StainlessSteel = new G4Material("StainlessSteel", 8.00*g/cm3, 3);
  StainlessSteel->AddElement(Fe, 70 * perCent);  // 70% Iron
  StainlessSteel->AddElement(Cr, 18 * perCent);  // 18% Chromium
  StainlessSteel->AddElement(Ni, 12 * perCent);  // 12% Nickel

  // Define Acetal for the detector holder
  G4Material* Acetal = new G4Material("Acetal", 1.41*g/cm3, ncomponents = 3);
  Acetal->AddElement(C, 1); // 1 carbon atom
  Acetal->AddElement(H, 2); // 2 hydrogen atoms
  Acetal->AddElement(O, 1); // 1 oxygen atom

  // Define Galactic Vacuum for the world volume
  G4Material* gal_vac = nist->FindOrBuildMaterial("G4_Galactic");

  // Define nickel for the 207Bi source backing
  G4Material* Nickel = new G4Material("Nickel", 8.908*g/cm3, 1);
  Nickel->AddElement(Ni, 1);

  // G4cout << G4endl << *(G4Material::GetMaterialTable()) << G4endl;

  // Default materials in setup.
  WorldMaterial = gal_vac;
  AttenuatorMaterial = Tantalum; // AC
  MagnetMaterial = NdFeB; // AC
  DetectorMaterial = silicon; // AC
  DetectorHousingMaterial = StainlessSteel; // AC
  DetectorWindowMaterial = silicon; // AC
  DetectorHolderMaterial = Acetal; // AC
  SourceBackingMaterial = StainlessSteel; // AC

  // G4cout << "end material"<< G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* ICESPICEDetectorConstruction::ConstructCalorimeter()
{  
  //The World
  WorldSizeXY  = 300.*mm;  // Cube
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

#if SOURCEBACKING
  Bi207SourceBacking();
#endif

ICESPICE();


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  return physiWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ICESPICEDetectorConstruction::ICESPICE()
{
  // Add the attenuator at the origin
  #if ATTENUATOR
    auto attenuator = CADMesh::TessellatedMesh::FromPLY("./cad_files/tantalum_5_slot_attenuator.PLY");
    solidAttenuator = attenuator->GetSolid();
    logicAttenuator = new G4LogicalVolume(solidAttenuator, AttenuatorMaterial, "Attenuator");
    // rotate 180 degrees to match the CAD file
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateY(180*deg);
    physiAttenuator = new G4PVPlacement(rot,			             //no rotation
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
  // Add the magnets around the origin with the same offset as in solidworks
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

  #if MAGNETHOLDER
      auto magnetHolder = CADMesh::TessellatedMesh::FromPLY("./cad_files/5N42_1x1x1_8in_magnets_mount.stl");
      auto solidMagnetHolder = magnetHolder->GetSolid();
      auto logicMagnetHolder = new G4LogicalVolume(solidMagnetHolder, DetectorHousingMaterial, "MagnetHolder");

      // Place the magnet holder at the origin
      auto physiMagnetHolder = new G4PVPlacement(0,			             //no rotation
              G4ThreeVector(0.,0.,0.), //at (0,0,0)
                                    "MagnetHolder",		             //its name
                                    logicMagnetHolder,		             //its logical volume
                                    physiWorld,			             //its mother  volume
                                    false,			                     //no boolean operation
                                    0);			                     //copy number

      // Visualization attributes
      G4VisAttributes* simpleMagnetHolderVisAtt= new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); //grey 
      simpleMagnetHolderVisAtt->SetVisibility(true);
      simpleMagnetHolderVisAtt->SetForceSolid(true);
      logicMagnetHolder->SetVisAttributes(simpleMagnetHolderVisAtt);
  #endif
}

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
      G4FieldManager* pFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
           
      G4cout<< "DeltaStep "<<pFieldMgr->GetDeltaOneStep()/mm <<"mm" <<G4endl;
      // G4ChordFinder *pChordFinder = new G4ChordFinder(ICESPICEField);

      pFieldMgr->SetDetectorField(fField.Get());
      pFieldMgr->CreateChordFinder(fField.Get());
      
    }
#endif
}

void ICESPICEDetectorConstruction::DefineCommands() {

    fMessenger = new G4GenericMessenger(this, 
                                        "/ICESPICE/", 
                                        "Detector control");

    // Detector Position
    G4GenericMessenger::Command& detectorPosition
      = fMessenger->DeclareMethodWithUnit("DetectorPosition","mm",
                                  &ICESPICEDetectorConstruction::SetDetectorPosition, 
                                  "Set the position of the detector face");
    detectorPosition.SetParameterName("position", true);
    detectorPosition.SetRange("position>-100. && position<=0.");
    detectorPosition.SetDefaultValue("-30.");


    // Source Position
    G4GenericMessenger::Command& sourcePosition
      = fMessenger->DeclareMethodWithUnit("SourcePosition","mm",
                                  &ICESPICEDetectorConstruction::SetSourcePosition, 
                                  "Set the position of the source");
    sourcePosition.SetParameterName("position", true);
    sourcePosition.SetDefaultValue("70.");

}

void ICESPICEDetectorConstruction::SetDetectorPosition(G4double val) {
    DetectorPosition = val;
    physiDetector->SetTranslation(G4ThreeVector(0, 0, DetectorPosition-DetectorThickness/2.));
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void ICESPICEDetectorConstruction::PIPS1000Detector() {
  // Assuming that the detector window and housing are positioned relative to the detector's dimensions.
  // auto detector = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips1000/active_area.PLY");

  // solidDetector = detector->GetSolid();
  // logicDetector = new G4LogicalVolume(solidDetector,
  //                                       DetectorMaterial,
  //                                       "Detector");


  DetectorActiveArea = 49.9*mm2; // Active area of the detector
  DetectorThickness = 1000.*micrometer; // Thickness of the detector
  DetectorWindowThickness = 50.*nanometer; // Thickness of the detector window
  G4double DetectorRadius = std::sqrt(DetectorActiveArea / 3.14);

    // Create the cylindrical detector (G4Tubs)
  solidDetector = new G4Tubs("Detector",
                             0.,                   // Inner radius
                             DetectorRadius,        // Outer radius
                             DetectorThickness/2.,  // Half thickness
                             0.*deg,               // Starting angle
                             360.*deg);            // Spanning angle
  
  // Create the logical volume for the detector
  logicDetector = new G4LogicalVolume(solidDetector,
                                      DetectorMaterial,
                                      "Detector");

  solidDetectorWindow = new G4Tubs("DetectorWindow",
                                    0,  // Inner radius
                                    std::sqrt(DetectorActiveArea / 3.14),  // Outer radius
                                    DetectorWindowThickness / 2.,  // Half-height
                                    0.*deg,  // Start angle
                                    360.*deg);  // Spanning angle

  logicDetectorWindow = new G4LogicalVolume(solidDetectorWindow,
                                      DetectorWindowMaterial,
                                      "DetectorWindow");


  // Recalculate the position if it's dependent on the detector's thickness
  G4double windowZPosition = - DetectorWindowThickness / 2. + DetectorThickness / 2.;

  physiDetectorWindow = new G4PVPlacement(nullptr,  // No rotation
              G4ThreeVector(0, 0, windowZPosition),  // Position in the detector
              logicDetectorWindow,
              "DetectorWindow",
              logicDetector,  // Parent volume
              false,  // No boolean operation
              0);  // Copy number

  // Create the outer housing for the detector
  auto detectorHousing = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips1000/detector_housing.PLY");
  solidDetectorHousing = detectorHousing->GetSolid();
  logicDetectorHousing = new G4LogicalVolume(solidDetectorHousing,
                                            AttenuatorMaterial,
                                            "DetectorHousing");

                        // Place the detector within the housing
  physiDetectorHousing = new G4PVPlacement(nullptr,  // No rotation
                    G4ThreeVector(0, 0, DetectorThickness/2.),  // Position relative to housing center
                    logicDetectorHousing,
                    "DetectorHousing",
                    logicDetector,  // Parent volume
                    false,  // No boolean operation
                    0);  // Copy number

  // create the detector holder 
  #if DETECTORHOLDER
    auto detectorHolder = CADMesh::TessellatedMesh::FromPLY("./cad_files/PIPS_holder.PLY");
    solidDetectorHolder = detectorHolder->GetSolid();
    logicDetectorHolder = new G4LogicalVolume(solidDetectorHolder,
                                              DetectorHolderMaterial,
                                              "DetectorHolder");

    // Place the holder at the origin of the detector volume
    physiDetectorHolder = new G4PVPlacement(nullptr,  // No rotation
                    G4ThreeVector(0, 0, 1.0*mm),  // Position relative to housing center
                    logicDetectorHolder,
                    "DetectorHolder",
                    logicDetector,  // Parent volume
                    false,  // No boolean operation
                    0);  // Copy number
  #endif

  physiDetector = new G4PVPlacement(nullptr,  // no rotation
              G4ThreeVector(0, 0, DetectorPosition-DetectorThickness/2.),  // position in world
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

  #if DETECTORHOLDER
    G4VisAttributes* visAttributesHolder = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));  // blue for the holder
    visAttributesHolder->SetVisibility(false);
    visAttributesHolder->SetForceSolid(false);
    logicDetectorHolder->SetVisAttributes(visAttributesHolder);
  #endif

  }

void ICESPICEDetectorConstruction::Bi207SourceBacking() {
  // Create the source backing
  SourceBackingThickness = 0.25*25.4*mm; // 0.5 inches;
  G4double SourceBackingDiameter = 1.0*25.4*mm; // 1 inch;

  // create a cylinder that is 1 inch in diameter and 0.5 in thick

  solidSourceBacking=  new G4Tubs("SourceBacking",
                                  0,  // Inner radius
                                  SourceBackingDiameter/2.0,  // Outer radius
                                  SourceBackingThickness/2.0,  // Half-height
                                  0.*deg,  // Start angle
                                  360.*deg);  // Spanning angle

  logicSourceBacking = new G4LogicalVolume(solidSourceBacking, SourceBackingMaterial, "SourceBacking");

  // // add an ring with a hole in the center with a diameter of 0.75 in and and outer diameter of 1 in
  G4double ringThickness = 0.5*mm; // 0.125 inches
  G4double ringInnerDiameter = 0.75*25.4*mm; // 0.75 inches
  G4double ringOuterDiameter = 1.0*25.4*mm; // 1 inch

  G4Tubs* solidRing = new G4Tubs("SourceRing",
                                  ringInnerDiameter/2.0,  // Inner radius
                                  ringOuterDiameter/2.0,  // Outer radius
                                  ringThickness/2.0,  // Half-height
                                  0.*deg,  // Start angle
                                  360.*deg);  // Spanning angle

  G4LogicalVolume* logicRing = new G4LogicalVolume(solidRing, SourceBackingMaterial, "SourceRing");

  // place the ring at the center of the source backing
  new G4PVPlacement(0,  // No rotation
              G4ThreeVector(0, 0, - SourceBackingThickness/2.0 - ringThickness/2.0 + 1.0*micrometer ),  // Position in the world
              logicRing,
              "SourceRing",
              logicSourceBacking,  // Parent volume
              false,  // No boolean operation
              0,  // Copy number
              false);

  physiSourceBacking = new G4PVPlacement(0,  // No rotation
              G4ThreeVector(0, 0, SourcePosition + SourceBackingThickness/2.0),  // Position in the world
              logicSourceBacking,
              "SourceBacking",
              logicWorld,  // Parent volume
              false,  // No boolean operation
              0);  // Copy number

  // visualization attributes
  G4VisAttributes* simpleSourceBackingVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));  // grey
  simpleSourceBackingVisAtt->SetVisibility(true);
  simpleSourceBackingVisAtt->SetForceSolid(true);
  logicSourceBacking->SetVisAttributes(simpleSourceBackingVisAtt);

  G4VisAttributes* simpleRingVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));  // grey
  simpleRingVisAtt->SetVisibility(true);
  simpleRingVisAtt->SetForceSolid(true);
  logicRing->SetVisAttributes(simpleRingVisAtt);

}

void ICESPICEDetectorConstruction::SetSourcePosition(G4double val) {

    G4double halfTragetThickness = 500.0*nanometer; 
    SourcePosition = val - halfTragetThickness;

    // Now, use the same value to update the GPS position
    G4UImanager* UI = G4UImanager::GetUIpointer();
    std::ostringstream command;
    command << "/gps/pos/centre 0 0 " << SourcePosition << " mm";  // Assuming the GPS is along the z-axis
    UI->ApplyCommand(command.str());

    G4cout << "Source position set to: " << SourcePosition << G4endl;
    G4cout << "GPS center set to: (0, 0, " << SourcePosition << ") mm" << G4endl;

    G4double BackingPosition = val + SourceBackingThickness/2.0 + 0.01*micrometer;

    physiSourceBacking->SetTranslation(G4ThreeVector(0, 0, BackingPosition));
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}