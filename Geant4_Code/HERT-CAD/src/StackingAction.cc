#include "StackingAction.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

G4int StackingAction::nGamma[21]    = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
G4int StackingAction::nElectron[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
G4int StackingAction::nPositron[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


G4double StackingAction::eMinGamma[21] = {
	DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,
	DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,
	DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX			};
	
G4double StackingAction::eMinElectron[21] = {
	DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,
	DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,
	DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX		};

G4double StackingAction::eMinPositron[21] = {
	DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,
	DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,
	DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX			};
	
G4double StackingAction::eMaxGamma[21]    = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
G4double StackingAction::eMaxElectron[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
G4double StackingAction::eMaxPositron[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};		

G4String StackingAction::detName[9] = {
    "detector_1","detector_2",
    "detector_3","detector_4",
    "detector_5","detector_6",
    "detector_7","detector_8", "detector_9" };

//G4String StackingAction::shldAlName[9] = {
//	          "AlTube","AlBackPlate","AlWasher1","AlWasher2","AlWasher3","AlWasher4",
// 	                         "AlWasher5","AlWasher6","AlWasher7"};
//
//G4String StackingAction::shldWName[10] = {
//		          "WTube","Wann","WDisc","Washer1","Washer2","Washer3","Washer4",
// 	                         "Washer5","Washer6","Washer7"};
//
/**
 * Note the names for Al and W  are names of "Physical Volumes" 
 * and are therefore different from  "Sensitive det names"    
 */
StackingAction::StackingAction()
{;}

StackingAction::~StackingAction()
{;}

G4ClassificationOfNewTrack 
                           StackingAction::ClassifyNewTrack(const G4Track * aTrack){
                           	
    
    G4TouchableHistory* theTouchable
   = (G4TouchableHistory*)(aTrack->GetTouchable());

//-------
   
   if(theTouchable){
//  	G4cout << " <StackingAction::ClassifyNewTrack> HistoryDepth " << 
//  	          theTouchable->GetHistoryDepth() << G4endl;
  	 G4ThreeVector vv = aTrack->GetPosition();
  	 G4double kE = aTrack->GetKineticEnergy();
 
     G4ParticleDefinition* pdef = aTrack->GetDefinition() ;
  	 
//     G4cout << " <StackingAction::ClassifyNewTrack> physical Volume " << 
//   	        theTouchable->GetVolume(0)->GetName()  << " ID "
//   	        <<  aTrack->GetTrackID()  << " particle " 
//            <<  pdef->GetParticleName()<< " KE " << kE
//   	        << G4endl;
   }

//-------

  if(theTouchable&&theTouchable->GetHistoryDepth()>=0) {
//  	G4cout << " <StackingAction::ClassifyNewTrack> HistoryDepth " << 
//  	          theTouchable->GetHistoryDepth() << G4endl;
    G4VPhysicalVolume* calPhys = theTouchable->GetVolume(0);
    if(calPhys){
     G4int iDet = GetDetNumber(calPhys->GetName()) ;     
//     G4cout << " <StackingAction::ClassifyNewTrack> physical Volume " << 
//   	            calPhys->GetName() << " "<< iDet << G4endl;
    
     G4ParticleDefinition * particleType = aTrack->GetDefinition();
     G4double kEngy = aTrack->GetKineticEnergy();
    
     if(iDet >= 0){
       if(particleType==G4Gamma::GammaDefinition()){
          nGamma[iDet]++;
          if(eMinGamma[iDet]>kEngy) eMinGamma[iDet]=kEngy; 
          if(eMaxGamma[iDet]<kEngy) eMaxGamma[iDet]=kEngy; 
       }
       else if(particleType==G4Electron::ElectronDefinition()){
          nElectron[iDet]++;
          if(eMinElectron[iDet]>kEngy) eMinElectron[iDet]=kEngy; 
          if(eMaxElectron[iDet]<kEngy) eMaxElectron[iDet]=kEngy; 
       }
       else if(particleType==G4Positron::PositronDefinition()){
         nPositron[iDet]++;
         if(eMinPositron[iDet]>kEngy) eMinPositron[iDet]=kEngy; 
         if(eMaxPositron[iDet]>kEngy) eMaxPositron[iDet]=kEngy;         
       }
     }
  }

  } 
   return fUrgent;
}

void StackingAction::PrepareNewEvent(){ 

  for(int i=0;i<21;i++){
    nGamma[i] = 0; 
    nElectron[i] = 0;
    nPositron[i] = 0;
    eMinGamma[i] = DBL_MAX;
    eMinElectron[i] = DBL_MAX;
    eMinPositron[i] = DBL_MAX;
    eMaxGamma[i] = 0.0;
    eMaxElectron[i] = 0.0;
    eMaxPositron[i] = 0.0;    
  }
     
}
G4int StackingAction::GetDetNumber(G4String thisDetName){
	  G4int ik = -99;
	  for(int i=0;i<18;i++){
	  	if(thisDetName.find(detName[i]) < 10){
	     ik = i;		
	  	 break;	
	  	}
	  }
      if(thisDetName.find("Al") < 10)ik = 18;		
	  if(thisDetName.find("W" ) < 10)ik = 19;	
	  if(thisDetName.find("Be") < 10)ik = 20;	
	  
//	  G4cout << " <StackingAction> det/shld name  " << thisDetName << " ik "<< ik << G4endl;
     
	  return ik;	  		  
}    
