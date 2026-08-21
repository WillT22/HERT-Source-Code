#ifndef REPThit_h_
#define REPThit_h_ 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class REPTILEhit : public G4VHit 
{
 public:
   REPTILEhit();
   virtual ~REPTILEhit();
   REPTILEhit(const REPTILEhit&);
   const REPTILEhit& operator=(const REPTILEhit&);
   int operator==(const REPTILEhit&) const;

   inline void* operator new(size_t);
   inline void  operator delete(void*);

   virtual void Draw();
   virtual void Print();
      
 public:
   inline void AddEnergy(G4double de)
   { Edep += de; }
   
   inline void AddStep(G4double dl){
     TrackLength += dl; 
     nStep++;
   }
   
   inline G4double GetEdep() const{ return Edep; }
   
   inline G4double GetTrak() const{ return TrackLength; }
   
   inline G4int GetNStep() const{ return nStep; }
   
    
 private:
   G4double Edep;
   G4double TrackLength;
   G4int    nStep;
   G4ThreeVector pos;
   
};

typedef G4THitsCollection<REPTILEhit> REPTILEhitsCollection;

extern G4Allocator<REPTILEhit> REPTILEhitAllocator;

inline void* REPTILEhit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) REPTILEhitAllocator.MallocSingle();
  return aHit;
}

inline void REPTILEhit::operator delete(void* aHit)
{
  REPTILEhitAllocator.FreeSingle((REPTILEhit*) aHit);
}


#endif /*REPTILEbkstkKHIT_HH_*/
