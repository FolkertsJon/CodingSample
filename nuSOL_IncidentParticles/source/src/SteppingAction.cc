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
/// \file nuSolIncidentParticles/nuSolIncidentParticles/src/SteppingAction.cc
/// \brief Implementation of the nuSolIncidentParticles::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

using namespace nuSolIncidentParticles;

namespace nuSolIncidentParticles
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(const DetectorConstruction* detConstruction,
                               EventAction* eventAction)
  : fDetConstruction(detConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step

  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  auto postVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }

  if ( volume == fDetConstruction->GetAbsorberPV() ) {
    fEventAction->AddAbs(edep,stepLength);
  }

  if ( volume == fDetConstruction->GetGapPV() ) {
    fEventAction->AddGap(edep,stepLength);
  }

  if ( volume == fDetConstruction->GetHeatShieldPV()
       || volume == fDetConstruction->GetRadShieldPV()) {
    fEventAction->SetTouchedShield();
  }

  if ( volume == fDetConstruction->GetGAGGPV() ) {
    fEventAction->AddGAGG(edep,stepLength);
    if ( postVolume == fDetConstruction->GetVetoPV() ) {
      // get analysis manager
      auto analysisManager = G4AnalysisManager::Instance();

      G4ThreeVector pos = step->GetTrack()->GetPosition();
      G4ThreeVector mom = step->GetTrack()->GetMomentum();
      
      // fill ntuple
      analysisManager->FillNtupleIColumn(1, 0, step->GetTrack()->GetDefinition()->GetPDGEncoding());
      analysisManager->FillNtupleIColumn(1, 1, fEventAction->GetEventNumber());
      analysisManager->FillNtupleDColumn(1, 2, pos.x());
      analysisManager->FillNtupleDColumn(1, 3, pos.y());
      analysisManager->FillNtupleDColumn(1, 4, pos.z());
      analysisManager->FillNtupleDColumn(1, 5, mom.x());
      analysisManager->FillNtupleDColumn(1, 6, mom.y());
      analysisManager->FillNtupleDColumn(1, 7, mom.z());
      analysisManager->AddNtupleRow(1);
      // add the particle saving here
    }
  }

  if ( volume == fDetConstruction->GetVetoPV() ) {
    fEventAction->AddVeto(edep,stepLength);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
