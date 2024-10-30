#ifndef Hit_h_
#define Hit_h_ 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class Hit : public G4VHit 
{
 public:
   Hit();
   virtual ~Hit();
   Hit(const Hit&);
   const Hit& operator=(const Hit&);
   int operator==(const Hit&) const;

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

typedef G4THitsCollection<Hit> HitsCollection;

extern G4Allocator<Hit> HitAllocator;

inline void* Hit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) HitAllocator.MallocSingle();
  return aHit;
}

inline void Hit::operator delete(void* aHit)
{
  HitAllocator.FreeSingle((Hit*) aHit);
}


#endif 