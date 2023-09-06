#include "Hit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<Hit> HitAllocator;

Hit::Hit()
:Edep(0.),TrackLength(0.),nStep(0)
{;}

Hit::~Hit()
{;}

Hit::Hit(const Hit& right):G4VHit()
{
  Edep = right.Edep;
  TrackLength = right.TrackLength;
  nStep= right.nStep;
}

const Hit& Hit::operator=(const Hit& right)
{
  Edep = right.Edep;
  TrackLength = right.TrackLength;
  nStep= right.nStep;
  return *this;
}

int Hit::operator==(const Hit& right) const
{
  return (this==&right);
}

void Hit::Draw()
{
 
 }

void Hit::Print()
{;}

