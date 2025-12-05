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

#ifndef nuSolIncidentParticlesDetectorConstruction_h
#define nuSolIncidentParticlesDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;
class G4Material;

namespace nuSolIncidentParticles
{

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined
/// via G4GlobalMagFieldMessenger class.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

  public:
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // get methods
    //
    const G4VPhysicalVolume* GetAbsorberPV() const;
    const G4VPhysicalVolume* GetGapPV() const;
    const G4VPhysicalVolume* GetGAGGPV() const;
    const G4VPhysicalVolume* GetVetoPV() const;
    const G4VPhysicalVolume* GetRadShieldPV() const;
    const G4VPhysicalVolume* GetHeatShieldPV() const;

  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

  
  // Materials & Elements
  
  G4Material* fEljin_200;
  G4Material* fGd;
  G4Material* fGa;
  G4Material* fCe;
  G4Material* fAl;
  G4Material* fC;
  G4Material* fH;
  G4Material* fO;
  G4Material* fGAGG;
  G4Material* fCarbonAerogel;
  
    // data members
    //
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
                                      // magnetic field messenger

    G4VPhysicalVolume* fAbsorberPV = nullptr; // the absorber physical volume
    G4VPhysicalVolume* fGapPV = nullptr;      // the gap physical volume
  
    G4VPhysicalVolume* fGAGGPV = nullptr; // the GAG physical volume
    G4VPhysicalVolume* fVetoPV = nullptr;      // the veto physical volume
    G4VPhysicalVolume* fHeatShieldPV = nullptr; // the heat shield physical volume
    G4VPhysicalVolume* fRadShieldPV = nullptr;      // the radiation shield physical volume

    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps
};

// inline functions

inline const G4VPhysicalVolume* DetectorConstruction::GetAbsorberPV() const {
  return fAbsorberPV;
}

inline const G4VPhysicalVolume* DetectorConstruction::GetGapPV() const  {
  return fGapPV;
}

inline const G4VPhysicalVolume* DetectorConstruction::GetGAGGPV() const {
  return fGAGGPV;
}

inline const G4VPhysicalVolume* DetectorConstruction::GetVetoPV() const  {
  return fVetoPV;
}

inline const G4VPhysicalVolume* DetectorConstruction::GetRadShieldPV() const  {
  return fRadShieldPV;
}

inline const G4VPhysicalVolume* DetectorConstruction::GetHeatShieldPV() const  {
  return fHeatShieldPV;
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

