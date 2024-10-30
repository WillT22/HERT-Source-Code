//

#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"

class TrackingAction : public G4UserTrackingAction
{
  public:
    TrackingAction();
    virtual ~TrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
    
  private:  
    
    static G4double beamParticle[4];
    
  public:
    static G4double GetPrimaryParticle(G4int i) { return beamParticle[i]; }
    //static void     SetPrimaryParticle(G4int i, G4double val){ primaryParticle[i] = val; }
     
};

#endif
