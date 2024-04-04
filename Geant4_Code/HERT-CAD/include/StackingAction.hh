#ifndef StackingAction_H
#define StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4Track;

class StackingAction : public G4UserStackingAction
{
  public:
    StackingAction();
    virtual ~StackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void PrepareNewEvent();
    G4int GetDetNumber(G4String thisDetName);    
  private:

    static G4int nGamma[21];
    static G4int nElectron[21];
    static G4int nPositron[21];
    static G4double eMinGamma[21];
    static G4double eMinElectron[21];
    static G4double eMinPositron[21];
    static G4double eMaxGamma[21];
    static G4double eMaxElectron[21];
    static G4double eMaxPositron[21];
    static G4String detName[9];
//    static G4String shldAlName[6];
//    static G4String shldWName[6];

  public:
    static G4int GetNGamma(G4int i) { return nGamma[i]; }
    static G4int GetNElectron(G4int i) { return nElectron[i]; }
    static G4int GetNPositron(G4int i) { return nPositron[i]; }
    static G4double GetEMinGamma(G4int i) { return eMinGamma[i]; }
    static G4double GetEMinElectron(G4int i) { return eMinElectron[i]; }
    static G4double GetEMinPositron(G4int i) { return eMinPositron[i]; }
    static G4double GetEMaxGamma(G4int i) { return eMaxGamma[i]; }
    static G4double GetEMaxElectron(G4int i) { return eMaxElectron[i]; }
    static G4double GetEMaxPositron(G4int i) { return eMaxPositron[i]; }    
};

#endif

