
#ifndef SetSensDet_h
#define SetSensDet_h 1

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;
#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "REPTILEhit.hh"



class SetSensDet : public G4VSensitiveDetector
{
  public:
      SetSensDet(G4String);
      virtual ~SetSensDet();

      virtual void Initialize(G4HCofThisEvent*);
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();
      virtual void DrawAll();
      virtual void PrintAll();
      inline G4String GetThisName() const{return thisName;} 

  private:
      REPTILEhitsCollection*  SiDiscColl;      
      REPTILEhitsCollection*  shieldColl; 
      REPTILEhitsCollection*  BeDiscColl;
      G4String thisName; 
      G4int sensID;
};

#endif

