
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
//#include "G4ParticleGun.hh" // not needed with G4GeneralParticleSource
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
    fParticleGun(0)
{
    fParticleGun = new G4GeneralParticleSource();
}
//    G4int n_particle = 1;
//    particleGun = new G4ParticleGun(n_particle);
//
//    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
//    G4String particleName;
//    particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="e-"));
//    //particleGun->SetParticleDefinition(G4Geantino::GeantinoDefinition())
//    particleGun->SetParticleEnergy(2.0*MeV);
//    particleGun->SetParticlePosition(G4ThreeVector(0.0*cm, 0.0*cm, -4.0*cm));


	PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4int i = anEvent->GetEventID() % 3;
//	G4ThreeVector v(0.0,0.0,1.0);
//	switch(i) {
//	case 0:
//		break;
//	case 1:
//		v.setY(0.0);
//		break;
//	case 2:
//		v.setZ(0.0);
//		break;
//}
//	particleGun->SetParticleMomentumDirection(v);
fParticleGun->GeneratePrimaryVertex(anEvent);

}


