#include "PhysicsList.hh"
#include "G4ParticleTypes.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4SystemOfUnits.hh"

// --- Added Hadronic and HP Includes ---
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4NeutronTrackingCut.hh"


PhysicsList::PhysicsList(): G4VModularPhysicsList()
{
    // Particles
    particleList = new G4DecayPhysics("decays");

    // EM physics (Enforcing Option 4)
    emPhysicsList = new G4EmStandardPhysics_option4();
    
    // Hadronic physics (QGSP_BIC_HP and required companions)
    hadronPhys.push_back(new G4RadioactiveDecayPhysics());
    hadronPhys.push_back(new G4IonBinaryCascadePhysics());
    hadronPhys.push_back(new G4EmExtraPhysics());
    hadronPhys.push_back(new G4HadronElasticPhysicsHP());
    hadronPhys.push_back(new G4StoppingPhysics());
    hadronPhys.push_back(new G4HadronPhysicsQGSP_BIC_HP());
    hadronPhys.push_back(new G4NeutronTrackingCut());
}

PhysicsList::~PhysicsList()
{
    delete particleList;
    delete emPhysicsList;
    
    // Clean up hadronic physics list memory
    for(size_t i=0; i<hadronPhys.size(); i++) {
        delete hadronPhys[i];
    }
    hadronPhys.clear();
}

void PhysicsList::ConstructParticle()
{
    // In this method, static member functions should be called
    // for all particles which you want to use.
    particleList->ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
    // Define transportation process
    AddTransportation();
    
    // Construct EM and Particle processes
    emPhysicsList->ConstructProcess();
    particleList->ConstructProcess();
    
    // Construct Hadronic processes
    for(size_t i=0; i<hadronPhys.size(); i++) {
        hadronPhys[i]->ConstructProcess();
    }
}

void PhysicsList::SetCuts()
{
    // Suppress error messages even in case e/gamma/proton do not exist            
    G4int temp = GetVerboseLevel();                                                
    SetVerboseLevel(0);                                                           
    
    SetCutValue(0.1*mm, "electron");

    // Retrieve verbose level
    SetVerboseLevel(temp);  
}

void PhysicsList::AddPhysicsList(const G4String& name)
{
   if (name == "G4standard_exp") {
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option2();

  } else if (name == "G4standard_fast") {
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option1();
  }  
}

void PhysicsList::SetStandardList(G4bool flagHP, G4bool glauber)
{
    // Left empty as per your original file, initialization is now handled in the constructor
}