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
//
/// \file nuSolIncidentParticles/nuSolIncidentParticles/src/DetectorConstruction.cc
/// \brief Implementation of the nuSolIncidentParticles::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

namespace nuSolIncidentParticles
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_Fe");

  // Elements needed
   fGd = nistManager->FindOrBuildMaterial("G4_Gd");
   fGa = nistManager->FindOrBuildMaterial("G4_Ga");
   fCe = nistManager->FindOrBuildMaterial("G4_Ce");
   fAl = nistManager->FindOrBuildMaterial("G4_Al");
   fC = nistManager->FindOrBuildMaterial("G4_C");
   fH = nistManager->FindOrBuildMaterial("G4_H");
   fO = nistManager->FindOrBuildMaterial("G4_O");

  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density


  // Veto Material
  G4int ncomponents;
  fEljin_200 = new G4Material("Eljin", density = 1.02*g/cm3, ncomponents = 2);
  fEljin_200->AddMaterial(fC, 91.55* perCent);
  fEljin_200->AddMaterial(fH, 8.45* perCent);
  
  //GAGG(Ce) material
  fGAGG = new G4Material("GAGG", density = 6.63*g/cm3, ncomponents = 5);
  fGAGG->AddMaterial(fGd, 44.2  * perCent);
  fGAGG->AddMaterial(fAl, 5.06  * perCent);
  fGAGG->AddMaterial(fGa, 19.6  * perCent);
  fGAGG->AddMaterial(fO , 18.0  * perCent);
  fGAGG->AddMaterial(fCe, 13.14 * perCent);
  
  // Heat shield Material
  fCarbonAerogel = new G4Material("CarbonAerogel", density = 0.045*g/cm3, ncomponents = 1);
  fCarbonAerogel->AddMaterial(fC, 100* perCent);

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // My Parameters
  G4double GAGGRadius = 10*cm;
  G4double GAGGHalfZ = 0.5*m + 5*cm; // approximately 200 kg of GAGG + 5 cm for electronics

  G4double VetoRadius = GAGGRadius + 2*cm; // 2 cm larger in every direction
  G4double VetoHalfZ = GAGGHalfZ + 2*cm;

  // 9.36 cm is ~ CSDA range of 300 MeV p in Fe https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html
  // 1.16 is equivalent to 9.36 at the max angle coming directly from the sun at 8 solar radii
  // 1 cm for from behind
  G4double shieldRadius = VetoRadius + 1.16*cm;
  G4double shieldHalfZ = VetoHalfZ + (9.36 * cm + 1 * cm)/2.0;

  // dimensions match Parker
  G4double heatShieldHalfZ = 4.5*2.54*cm/2.0;
  G4double heatShieldRadius = 2.1*m;
  
  // Geometry parameters
  G4int nofLayers = 10;
  G4double absoThickness = 10.*mm;
  G4double gapThickness =  5.*mm;
  G4double calorSizeXY  = 10.*cm;

  auto layerThickness = absoThickness + gapThickness;
  auto calorThickness = nofLayers * layerThickness;
  auto worldSizeXY = 10*m;//1.2 * calorSizeXY;
  auto worldSizeZ  = 10*m;//1.2 * calorThickness;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  auto radShieldMaterial = G4Material::GetMaterial("G4_Fe");
  auto gapMaterial = G4Material::GetMaterial("liquidArgon");

  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name

  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                         // at (0,0,0)
    worldLV,                                 // its logical volume
    "World",                                 // its name
    nullptr,                                 // its mother  volume
    false,                                   // no boolean operation
    0,                                       // copy number
    fCheckOverlaps);                         // checking overlaps

  // Heat Shield
  G4Tubs *heatShield = new G4Tubs("heatShield", // name
				  0, // inner radius
				  heatShieldRadius, // outer radius
				  heatShieldHalfZ, // half z
				  0, // start angle
				  2*M_PI); // end angle
  
  G4LogicalVolume *heatShieldLog = new  G4LogicalVolume(
                 heatShield,     // its solid
                 fCarbonAerogel,  // its material
                 "heatShieldLog");   // its name
  
  fHeatShieldPV = new G4PVPlacement(nullptr,  // no rotation
		    G4ThreeVector(0, 0, -shieldHalfZ - gapThickness - heatShieldHalfZ),          // out from the shield.
		    heatShieldLog,                  // its logical volume
		    "heatShield",            // its name
		    worldLV,                  // its mother  volume
		    false,                    // no boolean operation
		    0,                        // copy number
		    fCheckOverlaps);          // checking overlaps

  
  // radiation Shield
  G4Tubs *radiationShield = new G4Tubs("radiationShield", // name
				  0, // inner radius
				  shieldRadius, // outer radius
				  shieldHalfZ, // half z
				  0, // start angle
				  2*M_PI); // end angle
  
  G4LogicalVolume *radiationShieldLog = new  G4LogicalVolume(
                 radiationShield,     // its solid
                 radShieldMaterial,  // its material
                 "radiationShieldLog");   // its name
  
  fRadShieldPV = new G4PVPlacement(nullptr,  // no rotation
		    G4ThreeVector(0, 0, 0),          // at 0,0,0
		    radiationShieldLog,                  // its logical volume
		    "radiationShield",            // its name
		    worldLV,                  // its mother  volume
		    false,                    // no boolean operation
		    0,                        // copy number
		    fCheckOverlaps);          // checking overlaps

  // Veto
  G4Tubs *Veto = new G4Tubs("Veto", // name
				  0, // inner radius
				  VetoRadius, // outer radius
				  VetoHalfZ, // half z
				  0, // start angle
				  2*M_PI); // end angle
  
  G4LogicalVolume *VetoLog = new  G4LogicalVolume(
                 Veto,     // its solid
                 fEljin_200,  // its material
                 "VetoLog");   // its name
  
  fVetoPV = new G4PVPlacement(nullptr,  // no rotation
	            G4ThreeVector(0, 0, -shieldHalfZ + VetoHalfZ + 9.36*cm),          // 15 cm from front of shield
		    VetoLog,                  // its logical volume
		    "Veto",            // its name
		    radiationShieldLog,                  // its mother  volume
		    false,                    // no boolean operation
		    0,                        // copy number
		    fCheckOverlaps);          // checking overlaps
  
  // GAGG
  G4Tubs *GAGG = new G4Tubs("GAGG", // name
				  0, // inner radius
				  GAGGRadius, // outer radius
				  GAGGHalfZ, // half z
				  0, // start angle
				  2*M_PI); // end angle
  
  G4LogicalVolume *GAGGLog = new  G4LogicalVolume(
                 GAGG,     // its solid
                 fEljin_200,  // its material
                 "GAGGLog");   // its name
  
  fGAGGPV = new G4PVPlacement(nullptr,  // no rotation
	            G4ThreeVector(0, 0, 0),          // 15 cm from front of shield
		    GAGGLog,                  // its logical volume
		    "GAGG",            // its name
		    VetoLog,                  // its mother  volume
		    false,                    // no boolean operation
		    0,                        // copy number
		    fCheckOverlaps);          // checking overlaps
  
  
  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  G4VisAttributes *heatShieldVis = new G4VisAttributes(G4Colour(1.0,1.0,1.0, 1.0));
  heatShieldVis->SetVisibility(true);
  heatShieldVis->SetForceSolid(true);
  heatShieldLog->SetVisAttributes(heatShieldVis);

  G4VisAttributes *radiationShieldVis = new G4VisAttributes(G4Colour(0.7,0.7,0.7, 0.5));
  radiationShieldVis->SetVisibility(true);
  radiationShieldVis->SetForceSolid(true);
  radiationShieldLog->SetVisAttributes(radiationShieldVis);

  G4VisAttributes *VetoVis = new G4VisAttributes(G4Colour(0.0,0.0,1, 0.5));
  VetoVis->SetVisibility(true);
  VetoVis->SetForceSolid(true);
  VetoLog->SetVisAttributes(VetoVis);

  G4VisAttributes *GAGGVis = new G4VisAttributes(G4Colour(1.0,0.9,0.0, 0.5));
  GAGGVis->SetVisibility(true);
  GAGGVis->SetForceSolid(true);
  GAGGLog->SetVisAttributes(GAGGVis);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

