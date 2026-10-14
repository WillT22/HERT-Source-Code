/*
written by R. Ian Crocker on 3 Mar 2008
Revised to support QGSP_BIC_HP hadronic vector arrays
*/

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include <vector> // Required for std::vector

class G4VPhysicsConstructor;

class PhysicsList: public G4VModularPhysicsList
{
  public:
    PhysicsList();
    virtual ~PhysicsList();

    // The AddPhysicsList method usually needs to be public if 
    // called from a macro or PhysicsListMessenger
    void AddPhysicsList(const G4String& name);
    void SetStandardList(G4bool flagHP = false, G4bool glauber = false);
    void List();

  protected:
    // Construct particles and physics processes
    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
    virtual void SetCuts() override;
   
  private:
    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;

    G4VPhysicsConstructor* emPhysicsList;
    G4VPhysicsConstructor* particleList;
    
    // Vector array to hold QGSP_BIC_HP and its required companions
    std::vector<G4VPhysicsConstructor*> hadronPhys;
    
    G4bool dump;  
};

#endif