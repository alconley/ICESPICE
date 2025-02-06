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
#define MAG 0

#define ICESPICE_5N42_1x1x1_8in_FLAG 0 
#define ICESPICE_3N42_1x1x1_16in_FLAG 0 
#define ICESPICE_5N42_1x1x1_16in_FLAG 0 
#define ICESPICE_6N42_1x1x1_16in_FLAG 0 

#define PIPS1000 1
#define PIPS500 0
#define PIPS300 0
#define PIPS100 0

#define DETECTORHOLDER 1 // AC: Volume for detector holder

#define BI207SOURCEBACKING 0

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ICESPICEDetectorConstruction::ICESPICEDetectorConstruction()
  :physiWorld(NULL), logicWorld(NULL), solidWorld(NULL),

    physiDetector(NULL), logicDetector(NULL), solidDetector(NULL),
    physiDetectorWindow(NULL), logicDetectorWindow(NULL), solidDetectorWindow(NULL),
    physiDetectorHousing(NULL), logicDetectorHousing(NULL), solidDetectorHousing(NULL),
    physiDetectorHolder(NULL), logicDetectorHolder(NULL), solidDetectorHolder(NULL),

    Aluminum(NULL), Silicon(NULL), StainlessSteel(NULL), Tantalum(NULL), Nylon(NULL),
    NdFeB(NULL), Acetal(NULL), Vacuum(NULL), SiO2(NULL), Bi(NULL),

    f207BiSourceEnabled(false), physi207BiSource(NULL), logic207BiSource(NULL), solid207BiSource(NULL), visAttributes207BiSource(NULL),

    fMessenger(0)  
{
  fField.Put(0);
  WorldSizeXY=WorldSizeZ=0;
  DetectorPosition=-30.*mm; // AC
  Source207BiPosition=70.*mm; // AC
  Source207BiThickness=1000.0*nanometer; // AC
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
  Bi = nist->FindOrBuildMaterial("G4_Bi");

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

  // // Define Galactic Vacuum for the world volume
  // G4Material* gal_vac = nist->FindOrBuildMaterial("G4_Galactic");
  // Vacuum = gal_vac;

  // define vacuum to be 10^-2 Torr
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* vacuum = new G4Material(name="vacuum", 1.33e-11 * g/cm3, ncomponents=1,
                                      kStateGas, 293.15 * kelvin, 1.333 * pascal);
  vacuum->AddMaterial(air, fractionmass=1.);
  Vacuum = vacuum;

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

if (f207BiSourceEnabled) {
    FSU207BiSource();
}

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

    // 207Bi Source Enable
    G4GenericMessenger::Command& toggleSourceCommand =
    fMessenger->DeclareMethod("FSU207BiSourceEnable",
                               &ICESPICEDetectorConstruction::ToggleFSU207BiSource,
                               "Enable or disable the Bi-207 source.");
    toggleSourceCommand.SetParameterName("enable", true);
    toggleSourceCommand.SetDefaultValue("false");

    // 207Bi Source Position
    G4GenericMessenger::Command& sourcePosition
      = fMessenger->DeclareMethodWithUnit("FSU207BiSourcePosition","mm",
                                  &ICESPICEDetectorConstruction::Set207BiSourcePosition, 
                                  "Set the position of the source");
    sourcePosition.SetParameterName("position", true);
    sourcePosition.SetDefaultValue("70.");

    // 207Bi Source Thickness
    G4GenericMessenger::Command& sourceThicknessCommand =
    fMessenger->DeclareMethodWithUnit("FSU207BiSourceThickness", "nanometer",
                                       &ICESPICEDetectorConstruction::UpdateFSU207BiSourceThickness,
                                       "Set the thickness of the Bi-207 source.");
    sourceThicknessCommand.SetParameterName("thickness", true);
    sourceThicknessCommand.SetRange("thickness > 0");
    sourceThicknessCommand.SetDefaultValue("1000.0");

}

void ICESPICEDetectorConstruction::SetDetectorPosition(G4double val) {
    DetectorPosition = val;
    physiDetector->SetTranslation(G4ThreeVector(0, 0, DetectorPosition-DetectorThickness/2.));
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
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
  // logicDetector = new G4LogicalVolume(solidDetector,
  //                                     Silicon,
  //                                     "Detector");
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

    // /*

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

  // */

                                                                                                   
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

void ICESPICEDetectorConstruction::FSU207BiSource() {
    // Define the geometry, material, and placement of the 207Bi source
    G4double sourceRadius = 2.5 * mm;    // Radius of the source
    G4double sourceBackingRadius = sourceRadius + 1.0 * mm; // Radius of the source backing
    G4double sourceBackingThickness = 10.0 * mm;            // Thickness of the source backing
    G4double holderRadius = sourceBackingRadius + 3.0 * mm; // Radius of the holder
    G4double holderLength = sourceBackingThickness + 0.5 * mm; // Length of the holder

    // Create a cylinder for the 207Bi source
    solid207BiSource = new G4Tubs("FSU207BiSource",
                                  0.0, sourceRadius,
                                  Source207BiThickness / 2,
                                  0.0 * deg, 360.0 * deg);

    logic207BiSource = new G4LogicalVolume(solid207BiSource,
                                           Bi, // Material: Bismuth
                                           "FSU207BiSource");

    // Position the source at the specified SourcePosition
    G4ThreeVector sourcePositionVector(0, 0, Source207BiPosition + Source207BiThickness / 2);

    physi207BiSource = new G4PVPlacement(nullptr,                // No rotation
                                         sourcePositionVector,   // Placement position
                                         logic207BiSource,       // Logical volume
                                         "FSU207BiSource",               // Name
                                         logicWorld,             // Parent volume
                                         false,                  // No boolean operation
                                         0);                     // Copy number

    // Set visualization attributes for the source
    visAttributes207BiSource = new G4VisAttributes(G4Colour(1.0, 0.5, 1.0)); // Light purple
    visAttributes207BiSource->SetVisibility(true);
    visAttributes207BiSource->SetForceSolid(true);
    logic207BiSource->SetVisAttributes(visAttributes207BiSource);

    // Add source backing behind the source
    G4Tubs* solidSourceBacking = new G4Tubs("SourceBacking",
                                            0.0, sourceBackingRadius,
                                            sourceBackingThickness / 2,
                                            0.0 * deg, 360.0 * deg);

    G4LogicalVolume* logicSourceBacking = new G4LogicalVolume(solidSourceBacking,
                                                              StainlessSteel, // Material for backing
                                                              "SourceBacking");

    // Position the backing directly behind the source
    G4ThreeVector backingPosition(0, 0, Source207BiThickness / 2 + sourceBackingThickness / 2);
    // G4ThreeVector backingPosition(0, 0, Source207BiThickness / 2 - sourceBackingThickness / 2);


    new G4PVPlacement(nullptr,              // No rotation
                      backingPosition,      // Placement position relative to the source
                      logicSourceBacking,   // Logical volume
                      "SourceBacking",      // Name
                      logic207BiSource,     // Parent volume (source)
                      false,                // No boolean operation
                      0);                   // Copy number

    // Set visualization attributes for the source backing
    G4VisAttributes* visAttributesSourceBacking = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3)); // Gray
    visAttributesSourceBacking->SetVisibility(true);
    visAttributesSourceBacking->SetForceSolid(true);
    logicSourceBacking->SetVisAttributes(visAttributesSourceBacking);

    // Add the holder housing around the source backing
    G4Tubs* solidHolder = new G4Tubs("Holder",
                                     sourceBackingRadius, holderRadius,
                                     holderLength / 2,
                                     0.0 * deg, 360.0 * deg);

    G4LogicalVolume* logicHolder = new G4LogicalVolume(solidHolder,
                                                       StainlessSteel, // Material for the holder
                                                       "Holder");

    // Position the holder to enclose the backing
    G4ThreeVector holderPosition(0, 0, Source207BiThickness / 2 + sourceBackingThickness / 2 - 0.25 * mm);

    new G4PVPlacement(nullptr,            // No rotation
                      holderPosition,     // Placement position relative to the source
                      logicHolder,        // Logical volume
                      "Holder",           // Name
                      logic207BiSource,   // Parent volume (source)
                      false,              // No boolean operation
                      0);                 // Copy number

    // Set visualization attributes for the holder
    G4VisAttributes* visAttributesHolder = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3)); // Dark gray
    visAttributesHolder->SetVisibility(true);
    visAttributesHolder->SetForceSolid(true);
    logicHolder->SetVisAttributes(visAttributesHolder);

    // Apply GPS commands for particle source
    G4UImanager* UI = G4UImanager::GetUIpointer();
    G4String command;

    // Set the source center at SourcePosition
    // Adjust the GPS center to match the source position
    std::ostringstream positionCommand;
    positionCommand << "/gps/pos/centre 0 0 " << (Source207BiPosition + Source207BiThickness / 2) / mm << " mm";
    UI->ApplyCommand(positionCommand.str());

    // Set isotropic angular distribution
    command = "/gps/ang/type iso";
    UI->ApplyCommand(command);

    // Set volume and shape of the source
    command = "/gps/pos/type Volume";
    UI->ApplyCommand(command);
    command = "/gps/pos/shape Cylinder";
    UI->ApplyCommand(command);
    command = "/gps/pos/radius 2.5 mm";
    UI->ApplyCommand(command);

    // Set the height of the cylinder based on the source thickness
    std::ostringstream halfZCommand;
    halfZCommand << "/gps/pos/halfz " << (Source207BiThickness) / nm << " nanometer";
    UI->ApplyCommand(halfZCommand.str());

    command = "/gps/pos/confine FSU207BiSource";
    UI->ApplyCommand(command);

    // Set particle type and ion definition
    command = "/gps/particle ion";
    UI->ApplyCommand(command);

    // Define the ion parameters (example: Z=83 for Bi, A=207, E=0)
    G4int Z = 83;  // Atomic number of Bi
    G4int A = 207; // Mass number
    G4double E = 0; // Excitation energy
    std::ostringstream ionCommand;
    ionCommand << "/gps/ion " << Z << " " << A << " 0 " << E;
    UI->ApplyCommand(ionCommand.str());

    // Set monoenergetic energy distribution
    command = "/gps/ene/type Mono";
    UI->ApplyCommand(command);
    command = "/gps/ene/mono 0 eV";
    UI->ApplyCommand(command);

    G4cout << "207Bi source geometry and GPS configuration updated successfully." << G4endl;
}

void ICESPICEDetectorConstruction::UpdateFSU207BiSourceThickness(G4double newThickness) {
    G4cout << "Updating 207Bi source thickness to: " << newThickness / nanometer << " nanometer" << G4endl;

    // Check if the thickness is actually changing
    if (Source207BiThickness == newThickness) {
        G4cout << "Thickness is already set to the given value. No update required." << G4endl;
        return;
    }

    // Update the thickness variable
    Source207BiThickness = newThickness;

    // Completely remove the old source geometry
    if (physi207BiSource) {
        G4cout << "Removing existing 207Bi source." << G4endl;
        delete physi207BiSource;
        physi207BiSource = nullptr;
    }

    if (logic207BiSource) {
        delete logic207BiSource;
        logic207BiSource = nullptr;
    }

    if (solid207BiSource) {
        delete solid207BiSource;
        solid207BiSource = nullptr;
    }

    if (visAttributes207BiSource) {
        delete visAttributes207BiSource;
        visAttributes207BiSource = nullptr;
    }

    // Recreate the 207Bi source with the updated thickness
    if (f207BiSourceEnabled) {
        G4cout << "Recreating 207Bi source with updated thickness." << G4endl;
        FSU207BiSource();
    }

    // Notify the RunManager about the geometry change
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void ICESPICEDetectorConstruction::ToggleFSU207BiSource(G4bool enable) {
    G4cout << "ToggleFSU207BiSource called. Enable: " << enable << G4endl;

    // Check if the state is changing
    if (f207BiSourceEnabled == enable) {
        G4cout << "No change in state. Returning." << G4endl;
        return;
    }

    f207BiSourceEnabled = enable;

    // If enabling the source, recreate it
    if (f207BiSourceEnabled) {
        G4cout << "Creating 207Bi source." << G4endl;
        FSU207BiSource();
    } 
    // If disabling the source, remove it
    else {
        G4cout << "Removing 207Bi source." << G4endl;

        // Remove the physical volume
        if (physi207BiSource) {
            delete physi207BiSource;
            physi207BiSource = nullptr;
        }

        // Remove the logical volume
        if (logic207BiSource) {
            delete logic207BiSource;
            logic207BiSource = nullptr;
        }

        // Remove the solid volume
        if (solid207BiSource) {
            delete solid207BiSource;
            solid207BiSource = nullptr;
        }

        // Remove the visualization attributes
        if (visAttributes207BiSource) {
            delete visAttributes207BiSource;
            visAttributes207BiSource = nullptr;
        }
    }

    // Notify the RunManager about geometry changes
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// void ICESPICEDetectorConstruction::Set207BiSourcePosition(G4double val) {
//     Source207BiPosition = val;

//     if (f207BiSourceEnabled) {
//         // Recreate the source with the new position
//         ToggleFSU207BiSource(false);
//         ToggleFSU207BiSource(true);
//     }
// }

void ICESPICEDetectorConstruction::Set207BiSourcePosition(G4double val) {
    Source207BiPosition = val;

    if (physi207BiSource) {
        physi207BiSource->SetTranslation(G4ThreeVector(0, 0, Source207BiPosition + Source207BiThickness / 2));
        G4cout << "Updated source position to: " << Source207BiPosition / mm << " mm" << G4endl;

        // Force geometry and visualization updates
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        G4RunManager::GetRunManager()->ReinitializeGeometry();
    } else {
        G4cout << "Warning: Source physical volume not initialized!" << G4endl;
    }
}