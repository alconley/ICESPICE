
#ifndef ICESPICEDetectorConstruction_h
#define ICESPICEDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Cache.hh"
#include "G4MagneticField.hh"

#include "G4Tubs.hh"  // Ensure this header is included for cylindrical volumes

class G4Box;
class G4Trd;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class ICESPICETabulatedField3D;

class G4GenericMessenger;

class G4Tubs; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ICESPICEDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  ICESPICEDetectorConstruction();
  ~ICESPICEDetectorConstruction();

public:  
     
  G4VPhysicalVolume* Construct();
  void ConstructSDandField();
     
public:

  void PrintDetectorParameters();
                    
  G4double GetWorldSizeXY()  const       {return WorldSizeXY;};
  G4double GetWorldSizeZ() const          {return WorldSizeZ;}; 
  void PIPS100Detector();
  void PIPS300Detector(); 
  void PIPS500Detector(); 
  void PIPS1000Detector();
  void ICESPICE_5N42_1x1x1_8in();
  void ICESPICE_3N42_1x1x1_16in();
  void ICESPICE_5N42_1x1x1_16in();
  void ICESPICE_6N42_1x1x1_16in();
  void Bi207SourceBacking();

  void SetDetectorPosition(G4double val); 
  G4double GetDetectorPosition() const {return DetectorPosition;}; 

  void SetSourcePosition(G4double val);

  const G4VPhysicalVolume* GetWorld() const          {return physiWorld;};           
  const G4VPhysicalVolume* GetMeasureVolume() const { return physiDetector; } 
  const G4VPhysicalVolume* GetSiliconPV() const { return physiDetector; } 
  const G4VPhysicalVolume* GetDetectorWindowPV() const { return physiDetectorWindow; }
  const G4VPhysicalVolume* GetDetectorHousingPV() const { return physiDetectorHousing; }
  const G4VPhysicalVolume* GetDetectorHolderPV() const { return physiDetectorHolder; }


private:

  G4double           SSD;
  G4double           zOffset;
  
  G4double           WorldSizeXY;
  G4double           WorldSizeZ;
  G4VPhysicalVolume* physiWorld;
  G4LogicalVolume*   logicWorld;  
  G4Box*             solidWorld;

  G4double           DetectorPosition;   
  G4double           SourcePosition;

  G4VPhysicalVolume* physiDetector; 
  G4LogicalVolume*   logicDetector;
  // G4VSolid*          solidDetector; // uncomment if using CADMesh
  G4Tubs*            solidDetector; 
  G4double           DetectorActiveArea;
  G4double           DetectorThickness;

  
  G4VPhysicalVolume* physiDetectorWindow;
  G4LogicalVolume*   logicDetectorWindow;
  G4Tubs*            solidDetectorWindow; 
  G4double           DetectorWindowThickness;

  G4VPhysicalVolume* physiDetectorHousing;
  G4LogicalVolume*   logicDetectorHousing;
  G4VSolid*          solidDetectorHousing;   

  G4VPhysicalVolume* physiDetectorHolder;
  G4LogicalVolume*   logicDetectorHolder;
  G4VSolid*          solidDetectorHolder;

  G4VPhysicalVolume* physiSourceBacking;
  G4LogicalVolume*   logicSourceBacking;
  G4Tubs*            solidSourceBacking;
  G4double           SourceBackingThickness;

  G4Material*        Aluminum;
  G4Material*        Silicon;
  G4Material*        StainlessSteel;
  G4Material*        Tantalum;
  G4Material*        Nylon;
  G4Material*        NdFeB;
  G4Material*        Acetal;
  G4Material*        Vacuum;  

  G4Cache<G4MagneticField*> fField;  //pointer to the thread-local fields

  G4GenericMessenger* fMessenger;  // Messenger for dynamic configuration

private:

  void DefineMaterials();
  void DefineCommands();
  G4VPhysicalVolume* ConstructCalorimeter();     
};


#endif


