#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    void printFront(const G4Event*);
    void resetArrays();
    G4int GetDetNumber(G4String thisDetName);
  public:  
    static G4String bkSnam[17];
    static G4String bkGnam[7]; 
    static G4String ftSnam[4]; 
    static G4String ftGnam[4]; 
    static G4String shlnam[22];
    static G4String detName[21];
    static G4String BeName;

  private:
    G4int SiSensID[17];
    G4int backSensID[7];        
    G4int backGardID[7];        
    G4int frontSensID[4]; 
    G4int frontGardID[4]; 
    G4int shieldID[12]; 
    G4int BeDiscID; 

    G4int    evtNum;
    G4int    numTrk;
    G4double engyDep[21];
    G4double trkLeng[21];    
    G4int    numGamma[21];
    G4int    numElec[21];    
    G4int    numPosr[21];
    G4double minKeGamm[21]; 
    G4double minKeElec[21]; 
    G4double minKePosr[21];
    G4double maxKeGamm[21]; 
    G4double maxKeElec[21]; 
    G4double maxKePosr[21];      
 };

#endif
