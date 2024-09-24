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
//  Modified by: Alex Conley
//    *****************************************
//    *                                       *
//    *    ICESPICEDetectorConstruction.hh     *
//    *                                       *
//    *****************************************
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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

  G4Material* GetWorldMaterial()         {return WorldMaterial;};
  G4Material* GetAttenuatorMaterial()           {return AttenuatorMaterial;}; 
  G4Material* GetMagnetMaterial()           {return MagnetMaterial;}; 

  void PIPS100Detector();
  void PIPS300Detector(); 
  void PIPS500Detector(); 
  void PIPS1000Detector();
  void ICESPICE();
  void PIPSTransmissionDetector(G4double activeArea, G4double thickness); 

  void SetDetectorPosition(G4double val); 
  G4double GetDetectorPosition() const {return DetectorPosition;}; 

  const G4VPhysicalVolume* GetWorld() const          {return physiWorld;};           
  const G4VPhysicalVolume* GetMeasureVolume() const { return physiDetector; } 
  const G4VPhysicalVolume* GetSiliconPV() const { return physiDetector; } 
  const G4VPhysicalVolume* GetAttenuatorPV() const { return physiAttenuator; }
  const G4VPhysicalVolume* GetDetectorWindowPV() const { return physiDetectorWindow; }
  const G4VPhysicalVolume* GetDetectorHousingPV() const { return physiDetectorHousing; }


private:

  G4double           SSD;
  G4double           zOffset;
  
  G4double           WorldSizeXY;
  G4double           WorldSizeZ;
  G4VPhysicalVolume* physiWorld;
  G4LogicalVolume*   logicWorld;  
  G4Box*             solidWorld;
  G4Material*        WorldMaterial;

  G4Material*        MagnetMaterial; 

  G4double           DetectorPosition;   

  G4VPhysicalVolume* physiDetector; 
  G4LogicalVolume*   logicDetector;
  // G4VSolid*          solidDetector; // uncomment if using CADMesh
  G4Tubs*            solidDetector; 
  G4Material*        DetectorMaterial;
  G4double           DetectorActiveArea;
  G4double           DetectorThickness;

  G4VPhysicalVolume* physiAttenuator;
  G4LogicalVolume*   logicAttenuator;
  G4VSolid*          solidAttenuator;
  G4Material*        AttenuatorMaterial; 
  
  G4VPhysicalVolume* physiDetectorWindow;
  G4LogicalVolume*   logicDetectorWindow;
  G4Tubs*            solidDetectorWindow; 
  G4Material*        DetectorWindowMaterial; 
  G4double           DetectorWindowThickness;

  G4VPhysicalVolume* physiDetectorHousing;
  G4LogicalVolume*   logicDetectorHousing;
  G4VSolid*          solidDetectorHousing;   
  G4Material*        DetectorHousingMaterial; 

  G4VPhysicalVolume* physiTransmissionDetector; 
  G4LogicalVolume*   logicTransmissionDetector;
  G4Tubs*          solidTransmissionDetector;

  G4Cache<G4MagneticField*> fField;  //pointer to the thread-local fields

  G4GenericMessenger* fMessenger;  // Messenger for dynamic configuration

private:

  void DefineMaterials();
  void DefineCommands();
  G4VPhysicalVolume* ConstructCalorimeter();     
};


#endif


