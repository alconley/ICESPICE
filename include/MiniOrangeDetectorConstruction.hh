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
//    *    MiniOrangeDetectorConstruction.hh     *
//    *                                       *
//    *****************************************
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef MiniOrangeDetectorConstruction_h
#define MiniOrangeDetectorConstruction_h 1

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
class MiniOrangeTabulatedField3D;

class G4GenericMessenger;

class G4Tubs; // AC

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MiniOrangeDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  MiniOrangeDetectorConstruction();
  ~MiniOrangeDetectorConstruction();

public:  
     
  G4VPhysicalVolume* Construct();
  void ConstructSDandField();
     
public:

  void PrintDetectorParameters();
                    
  G4double GetWorldSizeXY()  const       {return WorldSizeXY;};
  G4double GetWorldSizeZ() const          {return WorldSizeZ;}; 

  G4double GetAttenuatorVolumeRadius() const {return AttenuatorVolumeRadius;}; // AC
  G4double GetAttenuatorVolumeHeight() const {return AttenuatorVolumeHeight;}; // AC

  G4double GetMagnetWidth()         {return MagnetWidth;}; // AC
  G4double GetMagnetHeight()         {return MagnetHeight;}; // AC
  G4double GetMagnetLength()         {return MagnetLength;}; // AC

  G4Material* GetWorldMaterial()         {return WorldMaterial;};
  G4Material* GetAttenuatorMaterial()           {return AttenuatorMaterial;}; // AC
  G4Material* GetMagnetMaterial()           {return MagnetMaterial;}; // AC

  void UpdateDetectorComponents(); // AC
  void SetDetectorPosition(G4double val); // AC
  G4double GetDetectorPosition() const {return DetectorPosition;}; // AC
  G4double GetDetectorActiveArea() const {return DetectorActiveArea;}; // AC

  void SetDetectorThickness(G4double val); // AC
  G4double GetDetectorThickness() const {return DetectorThickness;}; // AC
  G4double GetDetectorWindowThickness() const {return DetectorWindowThickness;}; // AC
  G4double GetDetectorHousingThickness() const {return DetectorHousingThickness;}; // AC
  G4double GetDetectorHousingOuterDiameter() const {return DetectorHousingOuterDiameter;}; // AC

  const G4VPhysicalVolume* GetWorld() const          {return physiWorld;};           
  const G4VPhysicalVolume* GetAttenuator() const      {return physiAttenuator;}; // AC
  const G4VPhysicalVolume* GetMagnet() const      {return physiMagnet;}; // AC
  const G4VPhysicalVolume* GetMeasureVolume() const { return physiDetector; } // AC
  const G4VPhysicalVolume* GetSiliconPV() const { return physiDetector; } // AC

private:
  
  G4double           WorldSizeXY;
  G4double           WorldSizeZ;

  G4double           MeasureVolumeRadius; // AC
  G4double           MeasureVolumeHeight; // AC
  G4double           MeasureVolumePosition; // AC

  G4double           AttenuatorVolumeRadius; // AC
  G4double           AttenuatorVolumeHeight; // AC
  G4double           AttenuatorVolumePosition; // AC

  G4double           MagnetWidth; // AC
  G4double           MagnetHeight; // AC
  G4double           MagnetLength; // AC

  G4double           DetectorPosition; // AC  
  G4double           DetectorActiveArea; // AC
  G4double           DetectorThickness; // AC
  G4double           DetectorWindowThickness; // AC
  G4double           DetectorHousingThickness; // AC
  G4double           DetectorHousingOuterDiameter; // AC

  G4double           SSD;
  G4double           zOffset;

  G4VPhysicalVolume* physiWorld;
  G4LogicalVolume*   logicWorld;  
  G4Box*             solidWorld;
  
  G4VPhysicalVolume* physiMagnet;
  G4LogicalVolume*   logicMagnet;
  G4Box*             solidMagnet;

  G4VPhysicalVolume* physiAttenuator; // AC
  G4LogicalVolume*   logicAttenuator; // AC
  G4Tubs*            solidAttenuator; // AC

  G4VPhysicalVolume* physiDetector; // AC
  G4LogicalVolume*   logicDetector; // AC
  G4Tubs*            solidDetector; // AC

  G4VPhysicalVolume* physiDetectorWindow; // AC
  G4LogicalVolume*   logicDetectorWindow; // AC
  G4Tubs*            solidDetectorWindow; // AC

  G4VPhysicalVolume* physiDetectorHousing; // AC
  G4LogicalVolume*   logicDetectorHousing; // AC
  G4Tubs*            solidDetectorHousing; // AC

  G4Material*        WorldMaterial;
  G4Material*        AttenuatorMaterial; // AC
  G4Material*        MagnetMaterial; // AC
  G4Material*        DetectorMaterial; // AC
  G4Material*        DetectorHousingMaterial; // AC

  G4Cache<G4MagneticField*> fField;  //pointer to the thread-local fields

  G4GenericMessenger* fMessenger;  // Messenger for dynamic configuration

private:

  void DefineMaterials();
  void DefineCommands();
  G4VPhysicalVolume* ConstructCalorimeter();     
};


#endif


