//
//
//
//
//
//

#include "SetSensDet.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4TouchableHandle.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4THitsCollection.hh"

SetSensDet::SetSensDet(G4String name):G4VSensitiveDetector(name),sensID(-1){
  

  thisName = name;
  G4String collNam="SiDiscColl";
  
  if( thisName.find("Al") < 10  || thisName.find("W") < 10) collNam="shieldColl";  
  if( thisName.find("Be") < 10 ) collNam="BeDiscColl";
  
  collectionName.insert(collNam);

  G4cout <<  " **********name "<< thisName << G4endl;    
}
//--------------------------------------------------------------------------
SetSensDet::~SetSensDet()
{;}
//--------------------------------------------------------------------------

void SetSensDet::Initialize(G4HCofThisEvent* HCE){

//        G4cout << "   <REPTsetSensDet Initialize> collectionName " << collectionName[0] 
//               << "   sensID " << sensID << " name "<< GetThisName() << G4endl;
       
      if(  collectionName[0] == "SiDiscColl" ){
        SiDiscColl = new REPTILEhitsCollection(SensitiveDetectorName,
                                                          collectionName[0]); 
        SiDiscColl->insert(new REPTILEhit);
        if(sensID<0){ sensID = GetCollectionID(0); }
        HCE->AddHitsCollection(sensID, SiDiscColl); 
      }
      if(  collectionName[0] == "shieldColl" ){
        shieldColl = new REPTILEhitsCollection(SensitiveDetectorName,
                                                          collectionName[0]); 
        shieldColl->insert(new REPTILEhit);
        if(sensID<0){ sensID = GetCollectionID(0); }
        HCE->AddHitsCollection(sensID, shieldColl);
      }  
      if(  collectionName[0] == "BeDiscColl" ){
        BeDiscColl = new REPTILEhitsCollection(SensitiveDetectorName,
                                                          collectionName[0]); 
        BeDiscColl->insert(new REPTILEhit);
        if(sensID<0){ sensID = GetCollectionID(0); }
        HCE->AddHitsCollection(sensID, BeDiscColl);
      }  

        
//        G4cout << "   <REPTsetSensDet Initialize> collectionName " << collectionName[0] 
//               << "   sensID " << sensID << " name "<< GetThisName() << G4endl;


}
//--------------------------------------------------------------------------
G4bool SetSensDet::ProcessHits(G4Step* aStep,G4TouchableHistory*){
	
  G4double edep = aStep->GetTotalEnergyDeposit(); 
  G4double stepl = aStep->GetStepLength();
  G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();

  
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

  G4VPhysicalVolume* thisVol = theTouchable->GetVolume(0);
  
//  G4cout << " ==== REPTsetSensDet ==== <ProcessHits> vol name " 
//         << thisVol->GetName() << " <edep> MeV "<< edep/MeV << G4endl;

  G4int rpNumber = theTouchable->GetReplicaNumber();  
  G4String thisDetName = thisVol->GetName();
  
//  G4cout << " ==== REPTsetSensDet ==== <ProcessHits> vol name " 
//         << thisVol->GetName() << " thisDetName " 
//         << thisDetName        << " ReplicaNumber "
//         << rpNumber << G4endl;
      
  if( thisDetName.find("Al") < 10 || thisDetName.find("W") < 10    ){
  	//G4cout <<  " start shieldColl rpNumber "  << G4endl;   	   
               (*shieldColl)[rpNumber]->AddEnergy(edep);           
               if( particle==G4Electron::ElectronDefinition() ||
                   particle==G4Positron::PositronDefinition()    )
                                    (*shieldColl)[rpNumber]->AddStep(stepl);                              
  } else if( thisDetName.find("Be") < 10 ){
  	//G4cout <<  " start BediscColl rpNumber "  << G4endl;   	   
               (*BeDiscColl)[rpNumber]->AddEnergy(edep);           
               if( particle==G4Electron::ElectronDefinition() ||
                   particle==G4Positron::PositronDefinition()    )
                                    (*BeDiscColl)[rpNumber]->AddStep(stepl);                              
  		
  }else {
  	//G4cout << " start SiDiscColl " << G4endl;
               (*SiDiscColl)[rpNumber]->AddEnergy(edep);
               if( particle==G4Electron::ElectronDefinition() ||
                   particle==G4Positron::PositronDefinition()    )
                                    (*SiDiscColl)[rpNumber]->AddStep(stepl);
  	}
   
//  G4cout << " ==== REPTsetSensDet ==== <ProcessHits>   END " 
//         << G4endl;
 
	return true;

}
//--------------------------------------------------------------------------
void SetSensDet::EndOfEvent(G4HCofThisEvent*)
{;}
//--------------------------------------------------------------------------
void SetSensDet::clear()
{;} 
//--------------------------------------------------------------------------
void SetSensDet::DrawAll()
{;} 
//--------------------------------------------------------------------------
void SetSensDet::PrintAll()
{;} 
