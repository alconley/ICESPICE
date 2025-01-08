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
#define MAG 1

#define ICESPICE_5N42_1x1x1_8in_FLAG 1 
#define ICESPICE_3N42_1x1x1_16in_FLAG 0 
#define ICESPICE_5N42_1x1x1_16in_FLAG 0 
#define ICESPICE_6N42_1x1x1_16in_FLAG 0 

#define PIPS1000 1
#define PIPS500 0
#define PIPS300 0
#define PIPS100 0

#define DETECTORHOLDER 1 // AC: Volume for detector holder

#define BI207SOURCEBACKING 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICEDetectorConstruction::ICESPICEDetectorConstruction()
  :physiWorld(NULL), logicWorld(NULL), solidWorld(NULL),

    physiDetector(NULL), logicDetector(NULL), solidDetector(NULL),
    physiDetectorWindow(NULL), logicDetectorWindow(NULL), solidDetectorWindow(NULL),
    physiDetectorHousing(NULL), logicDetectorHousing(NULL), solidDetectorHousing(NULL),
    physiDetectorHolder(NULL), logicDetectorHolder(NULL), solidDetectorHolder(NULL),
    physiSourceBacking(NULL), logicSourceBacking(NULL), solidSourceBacking(NULL),

    Aluminum(NULL), Silicon(NULL), StainlessSteel(NULL), Tantalum(NULL), Nylon(NULL),
    NdFeB(NULL), Acetal(NULL), Vacuum(NULL), SiO2(NULL),

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
  
  // Define N42 Magnet material (Neodymium Iron Boron) for magnets
  G4Material* Nd_magnet = new G4Material(name="N42_Magnet", 7.5*g/cm3, ncomponents=3);
    Nd_magnet->AddElement(Nd, natoms=2);
    Nd_magnet->AddElement(Fe, natoms=14);
    Nd_magnet->AddElement(B, natoms=1);
  NdFeB = Nd_magnet;

  Aluminum = nist->FindOrBuildMaterial("G4_Al");
  Silicon = nist->FindOrBuildMaterial("G4_Si");

  G4Material* stainlesssteel = new G4Material("StainlessSteel", 8.00*g/cm3, 3);
    stainlesssteel->AddElement(Fe, 70 * perCent);  // 70% Iron
    stainlesssteel->AddElement(Cr, 18 * perCent);  // 18% Chromium
    stainlesssteel->AddElement(Ni, 12 * perCent);  // 12% Nickel
    StainlessSteel = stainlesssteel;

  Tantalum = nist->FindOrBuildMaterial("G4_Ta");


  G4Material* nylon = new G4Material("Nylon", 1.14 * g / cm3, ncomponents = 4);
    // Add elements to Nylon
    nylon->AddElement(C, 6);  // 6 Carbon atoms
    nylon->AddElement(H, 11); // 11 Hydrogen atoms
    nylon->AddElement(N, 1);  // 1 Nitrogen atom
    nylon->AddElement(O, 1);  // 1 Oxygen atom
  Nylon = nylon;

  // Define Acetal for the detector holder
  G4Material* acetal = new G4Material("Acetal", 1.41*g/cm3, ncomponents = 3);
  acetal->AddElement(C, 1); // 1 carbon atom
  acetal->AddElement(H, 2); // 2 hydrogen atoms
  acetal->AddElement(O, 1); // 1 oxygen atom
  Acetal = acetal;

  // Define Galactic Vacuum for the world volume
  G4Material* gal_vac = nist->FindOrBuildMaterial("G4_Galactic");
  Vacuum = gal_vac;

  // SiO2 for the detector window
  G4Material* sio2 = new G4Material("SiO2", 2.2*g/cm3, ncomponents = 2);
  sio2->AddElement(Si, 1);
  sio2->AddElement(O, 2);
  SiO2 = sio2;

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
				   Vacuum,	//its material
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
#if PIPS1000
  PIPS1000Detector();
#endif

#if PIPS500
  PIPS500Detector();
#endif

#if PIPS300
  PIPS300Detector();
#endif

#if PIPS100
  PIPS100Detector();
#endif

#if SOURCEBACKING
  Bi207SourceBacking();
#endif

#if ICESPICE_5N42_1x1x1_8in_FLAG
  ICESPICE_5N42_1x1x1_8in();
#endif

#if ICESPICE_3N42_1x1x1_16in_FLAG
  ICESPICE_3N42_1x1x1_16in();
#endif

#if ICESPICE_5N42_1x1x1_16in_FLAG
  ICESPICE_5N42_1x1x1_16in();
#endif

#if ICESPICE_6N42_1x1x1_16in_FLAG
  ICESPICE_6N42_1x1x1_16in();
#endif

#if BI207SOURCEBACKING
  Bi207SourceBacking();
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

      #if ICESPICE_5N42_1x1x1_8in_FLAG
        // G4MagneticField* ICESPICEField= new ICESPICETabulatedField3D("comsol_output_5N42_1x1x8in_x50_y50_z70_res1_2mm.mag", zOffset);
        G4MagneticField* ICESPICEField= new ICESPICETabulatedField3D("comsol_output_5N42_1x1x8in_with_mounts_x50_y50_z70_res1_2mm.mag", zOffset);

      #endif

      #if ICESPICE_3N42_1x1x1_16in_FLAG
        G4MagneticField* ICESPICEField= new ICESPICETabulatedField3D("comsol_output_3N42_1x1x16in_x50_y50_z70_res1_2mm.mag", zOffset);
      #endif

      #if ICESPICE_5N42_1x1x1_16in_FLAG
        G4MagneticField* ICESPICEField= new ICESPICETabulatedField3D("comsol_output_5N42_1x1x16in_x50_y50_z70_res1_2mm.mag", zOffset);
      #endif

      #if ICESPICE_6N42_1x1x1_16in_FLAG
        G4MagneticField* ICESPICEField= new ICESPICETabulatedField3D("comsol_output_6N42_1x1x16in_x50_y50_z70_res1_2mm.mag", zOffset);
      #endif
      
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

    if (physiSourceBacking) {
        G4double BackingPosition = val + SourceBackingThickness/2.0 + 0.01*micrometer;
        physiSourceBacking->SetTranslation(G4ThreeVector(0, 0, BackingPosition));
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
        G4RunManager::GetRunManager()->ReinitializeGeometry();
    }
}

void ICESPICEDetectorConstruction::PIPS1000Detector() {
  DetectorActiveArea = 50.0*mm2; // Active area of the detector
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
                                      Silicon,
                                      "Detector");

  solidDetectorWindow = new G4Tubs("DetectorWindow",
                                    0,  // Inner radius
                                    std::sqrt(DetectorActiveArea / 3.14),  // Outer radius
                                    DetectorWindowThickness / 2.,  // Half-height
                                    0.*deg,  // Start angle
                                    360.*deg);  // Spanning angle

  logicDetectorWindow = new G4LogicalVolume(solidDetectorWindow,
                                      SiO2,
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
                                            StainlessSteel,
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
                                              Acetal,
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
  visAttributesHousing->SetForceSolid(true);
  logicDetectorHousing->SetVisAttributes(visAttributesHousing);

  #if DETECTORHOLDER
    G4VisAttributes* visAttributesHolder = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));  // blue for the holder
    visAttributesHolder->SetVisibility(true);
    visAttributesHolder->SetForceSolid(true);
    logicDetectorHolder->SetVisAttributes(visAttributesHolder);
  #endif

  }

void ICESPICEDetectorConstruction::PIPS500Detector() {
  // Assuming that the detector window and housing are positioned relative to the detector's dimensions.
  // auto detector = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips1000/active_area.PLY");

  // solidDetector = detector->GetSolid();
  // logicDetector = new G4LogicalVolume(solidDetector,
  //                                       DetectorMaterial,
  //                                       "Detector");


  DetectorActiveArea = 49.9*mm2; // Active area of the detector
  DetectorThickness = 500.*micrometer; // Thickness of the detector
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
                                      Silicon,
                                      "Detector");

  solidDetectorWindow = new G4Tubs("DetectorWindow",
                                    0,  // Inner radius
                                    std::sqrt(DetectorActiveArea / 3.14),  // Outer radius
                                    DetectorWindowThickness / 2.,  // Half-height
                                    0.*deg,  // Start angle
                                    360.*deg);  // Spanning angle

  logicDetectorWindow = new G4LogicalVolume(solidDetectorWindow,
                                      Silicon,
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
  auto detectorHousing = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips500/detector_housing.PLY");
  solidDetectorHousing = detectorHousing->GetSolid();
  logicDetectorHousing = new G4LogicalVolume(solidDetectorHousing,
                                            StainlessSteel,
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
                                              Acetal,
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
  visAttributesHousing->SetForceSolid(true);
  logicDetectorHousing->SetVisAttributes(visAttributesHousing);

  #if DETECTORHOLDER
    G4VisAttributes* visAttributesHolder = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));  // blue for the holder
    visAttributesHolder->SetVisibility(true);
    visAttributesHolder->SetForceSolid(true);
    logicDetectorHolder->SetVisAttributes(visAttributesHolder);
  #endif

  }

void ICESPICEDetectorConstruction::PIPS300Detector() {
  // Assuming that the detector window and housing are positioned relative to the detector's dimensions.
  // auto detector = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips1000/active_area.PLY");

  // solidDetector = detector->GetSolid();
  // logicDetector = new G4LogicalVolume(solidDetector,
  //                                       DetectorMaterial,
  //                                       "Detector");


  DetectorActiveArea = 49.9*mm2; // Active area of the detector
  DetectorThickness = 300.*micrometer; // Thickness of the detector
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
                                      Silicon,
                                      "Detector");

  solidDetectorWindow = new G4Tubs("DetectorWindow",
                                    0,  // Inner radius
                                    std::sqrt(DetectorActiveArea / 3.14),  // Outer radius
                                    DetectorWindowThickness / 2.,  // Half-height
                                    0.*deg,  // Start angle
                                    360.*deg);  // Spanning angle

  logicDetectorWindow = new G4LogicalVolume(solidDetectorWindow,
                                      Silicon,
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
  auto detectorHousing = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips300/detector_housing.PLY");
  solidDetectorHousing = detectorHousing->GetSolid();
  logicDetectorHousing = new G4LogicalVolume(solidDetectorHousing,
                                            StainlessSteel,
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
                                              Acetal,
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
  visAttributesHousing->SetForceSolid(true);
  logicDetectorHousing->SetVisAttributes(visAttributesHousing);

  #if DETECTORHOLDER
    G4VisAttributes* visAttributesHolder = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));  // blue for the holder
    visAttributesHolder->SetVisibility(true);
    visAttributesHolder->SetForceSolid(true);
    logicDetectorHolder->SetVisAttributes(visAttributesHolder);
  #endif

  }

void ICESPICEDetectorConstruction::PIPS100Detector() {
  // Assuming that the detector window and housing are positioned relative to the detector's dimensions.

  DetectorActiveArea = 49.9*mm2; // Active area of the detector
  DetectorThickness = 100.*micrometer; // Thickness of the detector
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
                                      Silicon,
                                      "Detector");

  solidDetectorWindow = new G4Tubs("DetectorWindow",
                                    0,  // Inner radius
                                    std::sqrt(DetectorActiveArea / 3.14),  // Outer radius
                                    DetectorWindowThickness / 2.,  // Half-height
                                    0.*deg,  // Start angle
                                    360.*deg);  // Spanning angle

  logicDetectorWindow = new G4LogicalVolume(solidDetectorWindow,
                                      Silicon,
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
  auto detectorHousing = CADMesh::TessellatedMesh::FromPLY("./cad_files/pips100/detector_housing.PLY");
  solidDetectorHousing = detectorHousing->GetSolid();
  logicDetectorHousing = new G4LogicalVolume(solidDetectorHousing,
                                            StainlessSteel,
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
                                              Acetal,
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
    visAttributesHolder->SetVisibility(true);
    visAttributesHolder->SetForceSolid(true);
    logicDetectorHolder->SetVisAttributes(visAttributesHolder);
  #endif

  }

void ICESPICEDetectorConstruction::Bi207SourceBacking() {
  // Create the source backing
  SourceBackingThickness = 0.25*25.4*mm;;
  G4double SourceBackingDiameter = 1.0*25.4*mm; // 1 inch;

  // create a cylinder that is 1 inch in diameter and 0.5 in thick

  solidSourceBacking=  new G4Tubs("SourceBacking",
                                  0,  // Inner radius
                                  SourceBackingDiameter/2.0,  // Outer radius
                                  SourceBackingThickness/2.0,  // Half-height
                                  0.*deg,  // Start angle
                                  360.*deg);  // Spanning angle

  logicSourceBacking = new G4LogicalVolume(solidSourceBacking, StainlessSteel, "SourceBacking");

  // // add an ring with a hole in the center with a diameter of 0.75 in and and outer diameter of 1 in
  G4double ringThickness = 1.0*mm; // 0.125 inches
  G4double ringInnerDiameter = 0.6*25.4*mm; // 0.75 inches
  G4double ringOuterDiameter = 1.0*25.4*mm; // 1 inch

  G4Tubs* solidRing = new G4Tubs("SourceRing",
                                  ringInnerDiameter/2.0,  // Inner radius
                                  ringOuterDiameter/2.0,  // Outer radius
                                  ringThickness/2.0,  // Half-height
                                  0.*deg,  // Start angle
                                  360.*deg);  // Spanning angle

  G4LogicalVolume* logicRing = new G4LogicalVolume(solidRing, StainlessSteel, "SourceRing");

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
  G4VisAttributes* simpleSourceBackingVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  simpleSourceBackingVisAtt->SetVisibility(true);
  simpleSourceBackingVisAtt->SetForceSolid(true);
  logicSourceBacking->SetVisAttributes(simpleSourceBackingVisAtt);

  G4VisAttributes* simpleRingVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); 
  simpleRingVisAtt->SetVisibility(true);
  simpleRingVisAtt->SetForceSolid(true);
  logicRing->SetVisAttributes(simpleRingVisAtt);

}

void ICESPICEDetectorConstruction::ICESPICE_5N42_1x1x1_8in() {
    auto mesh = CADMesh::TessellatedMesh::FromOBJ("./cad_files/5N42_1x1x1_8in.obj");

    // attenuator: tantalum 5 slot attenuator
    auto solidAttenuator = mesh->GetSolid("attenuator_5_slot");
    auto logicAttenuator = new G4LogicalVolume(solidAttenuator, Tantalum, "Attenuator");
    G4VPhysicalVolume* physiAttenuator = new G4PVPlacement(0,			             //no rotation
        G4ThreeVector(0, 0, 0), //at (0,0,0)
                              "Attenuator",		             //its name
                              logicAttenuator,		             //its logical volume
                              physiWorld,			             //its mother  volume
                              false,			                     //no boolean operation
                              0);			                     //copy number

    G4VisAttributes* simpleAttenuatorVisAtt= new G4VisAttributes(G4Colour(1,0,0)); //red
    simpleAttenuatorVisAtt->SetVisibility(true);
    simpleAttenuatorVisAtt->SetForceSolid(true);
    logicAttenuator->SetVisAttributes(simpleAttenuatorVisAtt);   

    // magnets: 5 1"x1"x1/8" N42 magnets
    std::vector<std::string> magnetNames = {
        "1x1x1_8in_square_magnet",
        "1x1x1_8in_square_magnet001",
        "1x1x1_8in_square_magnet002",
        "1x1x1_8in_square_magnet003",
        "1x1x1_8in_square_magnet004"
    };

    for (int i = 0; i < magnetNames.size(); i++) {
        auto solidMagnet = mesh->GetSolid(magnetNames[i]);
        auto logicMagnet = new G4LogicalVolume(solidMagnet, NdFeB, "Magnet" + std::to_string(i));
        G4VPhysicalVolume* physiMagnet = new G4PVPlacement(0,			             //no rotation
            G4ThreeVector(0, 0, 0), //at (0,0,0)
                                "Magnet" + std::to_string(i),		             //its name
                                logicMagnet,		             //its logical volume
                                physiWorld,			             //its mother  volume
                                false,			                     //no boolean operation
                                0);			                     //copy number

        G4VisAttributes* simpleMagnetVisAtt= new G4VisAttributes(G4Colour(0,1,0)); //green    
        simpleMagnetVisAtt->SetVisibility(true);
        simpleMagnetVisAtt->SetForceSolid(true);
        logicMagnet->SetVisAttributes(simpleMagnetVisAtt);
    }

    /*

    // magnet holders: Aluminum  
    std::vector<std::string> magnetHolderNames = {
        "magnet_holder",
        "magnet_holder001",
        "magnet_holder002",
        "magnet_holder003",
        "magnet_holder004"
    };

    for (int i = 0; i < magnetHolderNames.size(); i++) {
        auto solidMagnetHolder = mesh->GetSolid(magnetHolderNames[i]);
        auto logicMagnetHolder = new G4LogicalVolume(solidMagnetHolder, Aluminum, "MagnetHolder" + std::to_string(i));
        G4VPhysicalVolume* physiMagnetHolder = new G4PVPlacement(0,			             //no rotation
            G4ThreeVector(0, 0, 0), //at (0,0,0)
                                "MagnetHolder" + std::to_string(i),		             //its name
                                logicMagnetHolder,		             //its logical volume
                                physiWorld,			             //its mother  volume
                                false,			                     //no boolean operation
                                0);			                     //copy number

        G4VisAttributes* simpleMagnetHolderVisAtt= new G4VisAttributes(G4Colour(0,0,1)); //blue
        simpleMagnetHolderVisAtt->SetVisibility(true);
        simpleMagnetHolderVisAtt->SetForceSolid(true);
        logicMagnetHolder->SetVisAttributes(simpleMagnetHolderVisAtt);
    }

    // magnet holder caps: Aluminum 
    std::vector<std::string> magnetHolderCapNames = {
        "magnet_holder_cap",
        "magnet_holder_cap001",
        "magnet_holder_cap002",
        "magnet_holder_cap003",
        "magnet_holder_cap004",
        "magnet_holder_cap005",
        "magnet_holder_cap006",
        "magnet_holder_cap007",
        "magnet_holder_cap008",
        "magnet_holder_cap009"
    };

    for (int i = 0; i < magnetHolderCapNames.size(); i++) {
        auto solidMagnetHolderCap = mesh->GetSolid(magnetHolderCapNames[i]);
        auto logicMagnetHolderCap = new G4LogicalVolume(solidMagnetHolderCap, Aluminum, "MagnetHolderCap" + std::to_string(i));
        G4VPhysicalVolume* physiMagnetHolderCap = new G4PVPlacement(0,			             //no rotation
            G4ThreeVector(0, 0, 0), //at (0,0,0)
                                "MagnetHolderCap" + std::to_string(i),		             //its name
                                logicMagnetHolderCap,		             //its logical volume
                                physiWorld,			             //its mother  volume
                                false,			                     //no boolean operation
                                0);			                     //copy number

        G4VisAttributes* simpleMagnetHolderCapVisAtt= new G4VisAttributes(G4Colour(0,0,1)); //blue
        simpleMagnetHolderCapVisAtt->SetVisibility(true);
        simpleMagnetHolderCapVisAtt->SetForceSolid(true);
        logicMagnetHolderCap->SetVisAttributes(simpleMagnetHolderCapVisAtt);
    }

    // 5- magnet_brace: Aluminum
    std::vector<std::string> magnetBraceNames = {
        "magnet_brace",
        "magnet_brace001",
        "magnet_brace002",
        "magnet_brace003",
        "magnet_brace004"
    };

    for (int i = 0; i < magnetBraceNames.size(); i++) {
        auto solidMagnetBrace = mesh->GetSolid(magnetBraceNames[i]);
        auto logicMagnetBrace = new G4LogicalVolume(solidMagnetBrace, Aluminum, "MagnetBrace" + std::to_string(i));
        G4VPhysicalVolume* physiMagnetBrace = new G4PVPlacement(0,			             //no rotation
            G4ThreeVector(0, 0, 0), //at (0,0,0)
                                "MagnetBrace" + std::to_string(i),		             //its name
                                logicMagnetBrace,		             //its logical volume
                                physiWorld,			             //its mother  volume
                                false,			                     //no boolean operation
                                0);			                     //copy number

        G4VisAttributes* simpleMagnetBraceVisAtt= new G4VisAttributes(G4Colour(0,0,1)); //blue
        simpleMagnetBraceVisAtt->SetVisibility(true);
        simpleMagnetBraceVisAtt->SetForceSolid(true);
        logicMagnetBrace->SetVisAttributes(simpleMagnetBraceVisAtt);
    }

    // 5 bolts: Stainless Steel
    std::vector<std::string> hardware = {
      "socket_head_cap_screw_ai_HX-SHCS_0.19-32x0.5625x0.5625-N",
      "socket_head_cap_screw_ai_HX-SHCS_0.19-32x0.5625x0.5625-N001",
      "socket_head_cap_screw_ai_HX-SHCS_0.19-32x0.5625x0.5625-N002",
      "socket_head_cap_screw_ai_HX-SHCS_0.19-32x0.5625x0.5625-N003",
      "socket_head_cap_screw_ai_HX-SHCS_0.19-32x0.5625x0.5625-N004"
    };

    for (int i = 0; i < hardware.size(); i++) {
        auto solidHardware = mesh->GetSolid(hardware[i]);
        auto logicHardware = new G4LogicalVolume(solidHardware, StainlessSteel, "Hardware" + std::to_string(i));
        G4VPhysicalVolume* physiHardware = new G4PVPlacement(0,			             //no rotation
            G4ThreeVector(0, 0, 0), //at (0,0,0)
                                "Hardware" + std::to_string(i),		             //its name
                                logicHardware,		             //its logical volume
                                physiWorld,			             //its mother  volume
                                false,			                     //no boolean operation
                                0);			                     //copy number

        G4VisAttributes* simpleHardwareVisAtt= new G4VisAttributes(G4Colour(1,1,0)); //yellow
        simpleHardwareVisAtt->SetVisibility(true);
        simpleHardwareVisAtt->SetForceSolid(true);
        logicHardware->SetVisAttributes(simpleHardwareVisAtt);
    }

    // 10: bolts and 10: nuts: Nylon

    std::vector<std::string> hardwareNylon = {
      "socket_head_cap_screw_ai_HX-SHCS_0.112-48x0.375x0.375-N",
      "socket_head_cap_screw_ai_HX-SHCS_0.112-48x0.375x0.375-N001",
      "socket_head_cap_screw_ai_HX-SHCS_0.112-48x0.375x0.375-N002",
      "socket_head_cap_screw_ai_HX-SHCS_0.112-48x0.375x0.375-N003",
      "socket_head_cap_screw_ai_HX-SHCS_0.112-48x0.375x0.375-N004",
      "socket_head_cap_screw_ai_HX-SHCS_0.112-48x0.375x0.375-N005",
      "socket_head_cap_screw_ai_HX-SHCS_0.112-48x0.375x0.375-N006",
      "socket_head_cap_screw_ai_HX-SHCS_0.112-48x0.375x0.375-N007",
      "socket_head_cap_screw_ai_HX-SHCS_0.112-48x0.375x0.375-N008",
      "socket_head_cap_screw_ai_HX-SHCS_0.112-48x0.375x0.375-N009",

      "machine_screw_nut_hex_ai_MSHXNUT_0.112-48-S-N",
      "machine_screw_nut_hex_ai_MSHXNUT_0.112-48-S-N001",
      "machine_screw_nut_hex_ai_MSHXNUT_0.112-48-S-N002",
      "machine_screw_nut_hex_ai_MSHXNUT_0.112-48-S-N003",
      "machine_screw_nut_hex_ai_MSHXNUT_0.112-48-S-N004",
      "machine_screw_nut_hex_ai_MSHXNUT_0.112-48-S-N005",
      "machine_screw_nut_hex_ai_MSHXNUT_0.112-48-S-N006",
      "machine_screw_nut_hex_ai_MSHXNUT_0.112-48-S-N007",
      "machine_screw_nut_hex_ai_MSHXNUT_0.112-48-S-N008",
      "machine_screw_nut_hex_ai_MSHXNUT_0.112-48-S-N009"
    };

    for (int i = 0; i < hardwareNylon.size(); i++) { 
        auto solidHardwareNylon = mesh->GetSolid(hardwareNylon[i]);
        auto logicHardwareNylon = new G4LogicalVolume(solidHardwareNylon, Nylon, "HardwareNylon" + std::to_string(i));
        G4VPhysicalVolume* physiHardwareNylon = new G4PVPlacement(0,			             //no rotation
            G4ThreeVector(0, 0, 0), //at (0,0,0)
                                "HardwareNylon" + std::to_string(i),		             //its name
                                logicHardwareNylon,		             //its logical volume
                                physiWorld,			             //its mother  volume
                                false,			                     //no boolean operation
                                0);			                     //copy number

        G4VisAttributes* simpleHardwareNylonVisAtt= new G4VisAttributes(G4Colour(1,1,0)); //yellow
        simpleHardwareNylonVisAtt->SetVisibility(true);
        simpleHardwareNylonVisAtt->SetForceSolid(true);
        logicHardwareNylon->SetVisAttributes(simpleHardwareNylonVisAtt);
    }

  */

                                                                                                   
}

void ICESPICEDetectorConstruction::ICESPICE_3N42_1x1x1_16in() {
  // Add the attenuator at the origin
    auto attenuator = CADMesh::TessellatedMesh::FromPLY("./cad_files/attenuator_3_slots_1_16in.PLY");
    auto solidAttenuator = attenuator->GetSolid();
    auto logicAttenuator = new G4LogicalVolume(solidAttenuator, Tantalum, "Attenuator");
    // rotate 180 degrees to match the CAD file
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateY(0*deg);

    G4VPhysicalVolume* physiAttenuator = new G4PVPlacement(rot,			             //no rotation
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

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    auto magnet = CADMesh::TessellatedMesh::FromPLY("./cad_files/1x1x1_16in_square_magnet.PLY");
    auto solidMagnet = magnet->GetSolid();
    auto logicMagnet = new G4LogicalVolume(solidMagnet, NdFeB, "Magnet");

    // Calculate the placement and rotation for each magnet
    G4double placementRadius = 3.38*mm;  // Adjusting for the corner to be at 3.5mm
    G4int numMagnets = 3;
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
}

void ICESPICEDetectorConstruction::ICESPICE_5N42_1x1x1_16in() {
  // Add the attenuator at the origin
    auto attenuator = CADMesh::TessellatedMesh::FromPLY("./cad_files/attenuator_5_slots_1_16in.PLY");
    auto solidAttenuator = attenuator->GetSolid();
    auto logicAttenuator = new G4LogicalVolume(solidAttenuator, Tantalum, "Attenuator");
    // rotate 180 degrees to match the CAD file
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateY(0*deg);

    G4VPhysicalVolume* physiAttenuator = new G4PVPlacement(rot,			             //no rotation
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

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    auto magnet = CADMesh::TessellatedMesh::FromPLY("./cad_files/1x1x1_16in_square_magnet.PLY");
    auto solidMagnet = magnet->GetSolid();
    auto logicMagnet = new G4LogicalVolume(solidMagnet, NdFeB, "Magnet");

    // Calculate the placement and rotation for each magnet
    G4double placementRadius = 3.38*mm;  // Adjusting for the corner to be at 3.5mm
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
}

void ICESPICEDetectorConstruction::ICESPICE_6N42_1x1x1_16in() {
  // Add the attenuator at the origin
    auto attenuator = CADMesh::TessellatedMesh::FromPLY("./cad_files/attenuator_6_slots_1_16in.PLY");
    auto solidAttenuator = attenuator->GetSolid();
    auto logicAttenuator = new G4LogicalVolume(solidAttenuator, Tantalum, "Attenuator");
    // rotate 180 degrees to match the CAD file
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateY(0*deg);

    G4VPhysicalVolume* physiAttenuator = new G4PVPlacement(rot,			             //no rotation
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

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
    auto magnet = CADMesh::TessellatedMesh::FromPLY("./cad_files/1x1x1_16in_square_magnet.PLY");
    auto solidMagnet = magnet->GetSolid();
    auto logicMagnet = new G4LogicalVolume(solidMagnet, NdFeB, "Magnet");

    // Calculate the placement and rotation for each magnet
    G4double placementRadius = 3.5*mm;  // Adjusting for the corner to be at 3.5mm
    G4int numMagnets = 6;
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
}