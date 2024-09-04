#include "REPTILEhit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<REPTILEhit> REPTILEhitAllocator;

REPTILEhit::REPTILEhit()
:Edep(0.),TrackLength(0.),nStep(0)
{;}

REPTILEhit::~REPTILEhit()
{;}

REPTILEhit::REPTILEhit(const REPTILEhit& right):G4VHit()
{
  Edep = right.Edep;
  TrackLength = right.TrackLength;
  nStep= right.nStep;
}

const REPTILEhit& REPTILEhit::operator=(const REPTILEhit& right)
{
  Edep = right.Edep;
  TrackLength = right.TrackLength;
  nStep= right.nStep;
  return *this;
}

int REPTILEhit::operator==(const REPTILEhit& right) const
{
  return (this==&right);
}

void REPTILEhit::Draw()
{
/**	
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager){
    G4Circle circle(pos);
    circle.SetScreenSize(1.4);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
    }
 */  
 }

void REPTILEhit::Print()
{;}

