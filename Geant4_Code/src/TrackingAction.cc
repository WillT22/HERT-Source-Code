//

#include "TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "HistoManager.hh"
#include "G4Track.hh"
#include "G4Trajectory.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "EventAction.hh"

G4double TrackingAction::beamParticle[4] = {-99.,-99.,-99.,-99.};

TrackingAction::TrackingAction()
{;}

TrackingAction::~TrackingAction()
{;}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
    fpTrackingManager->SetStoreTrajectory(EventAction::particleTrack);
    fpTrackingManager->SetTrajectory(new G4Trajectory(aTrack));

    G4ParticleDefinition* pdef = aTrack->GetDefinition();

    //  G4cout << " ==== TrackingAction::PreUserTrackingAction ====  "
    //        <<  aTrack->GetTrackID()   << " particle " 
    //        << pdef->GetParticleName() << " K E " 
    //        << aTrack->GetKineticEnergy() << G4endl;  

    if (aTrack->GetTrackID() == 1) {
        G4ThreeVector vv = aTrack->GetPosition();
        G4double kE = aTrack->GetKineticEnergy();
        G4ThreeVector direction = aTrack->GetMomentumDirection();
        G4double weight = aTrack->GetWeight();
        G4int pid = aTrack->GetDynamicParticle()->GetPDGcode();

        for (int i = 0; i < 3; i++)beamParticle[i] = vv(i);
        beamParticle[3] = kE;
        //     G4cout << " ==== TrackingAction::PreUserTrackingAction ====   pos "
        //            <<  " x " << vv(0)   <<  " y " << vv(1)  <<  " z " << vv(2)
        //            << " KE " << kE << G4endl;
        G4double x = vv.x(), y = vv.y(), z = vv.z();
        G4double theta = direction.theta(), phi = direction.phi();
        if (phi < 0.) phi += twopi;
        G4double cost = std::cos(theta);

        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

        G4double r = vv.mag();
        G4double dr = analysisManager->GetH1Width(2);
        G4double dv = 2 * twopi * r * r * dr;

        // fill histograms
        analysisManager->FillH1(1, kE);
        if (dv > 0.) analysisManager->FillH1(2, r, 1. / dv);
        analysisManager->FillH1(3, cost);
        analysisManager->FillH1(4, phi);
        analysisManager->FillH2(1, x, y);
        analysisManager->FillH2(2, y, z);
        analysisManager->FillH2(3, z, x);
        analysisManager->FillH2(4, phi, cost);
        analysisManager->FillH2(5, phi, theta);

        // fill ntuple
        analysisManager->FillNtupleIColumn(0, pid);
        analysisManager->FillNtupleDColumn(1, kE);
        analysisManager->FillNtupleDColumn(2, x);
        analysisManager->FillNtupleDColumn(3, y);
        analysisManager->FillNtupleDColumn(4, z);
        analysisManager->FillNtupleDColumn(5, theta);
        analysisManager->FillNtupleDColumn(6, phi);
        analysisManager->FillNtupleDColumn(7, weight);
        analysisManager->AddNtupleRow();
    }

    G4TouchableHistory* theTouchable
        = (G4TouchableHistory*)(aTrack->GetTouchable());
    
    G4VPhysicalVolume* calPhys = theTouchable->GetVolume(0);

    //  G4cout << calPhys->GetName() << " " << aTrack->GetTrackID() << " " 
    //         << pdef->GetParticleName()   << G4endl;
}



