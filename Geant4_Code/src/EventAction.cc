/**
     Store/print event data
*/
#include "EventAction.hh"
#include "Hit.hh"
#include "StackingAction.hh"
#include "TrackingAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4RunManager.hh"
#include "globals.hh"


G4String  EventAction::bkSnam[]  ={"detector_1","detector_2",
                                "detector_3","detector_4",
                                "detector_5","detector_6",
                                "detector_7","detector_8", "detector_9"};

G4String EventAction::bkGnam[] = {"R3Gard","R4Gard","R5Gard","R6Gard",
                                      "R7Gard","R8Gard","R9Gard"};
                                      
G4String EventAction::ftSnam[] = {"R1Sens","R2Sens"};   

G4String EventAction::ftGnam[] = {"R1Gard","R2Gard"};   

G4String EventAction::shlnam[] = {"solid_alHousing","solid_TaColl",
	                                  "solid_TaTooth","solid_WFr1","solid_WFr2","solid_WFr3","solid_WFr4","solid_WFr5","solid_WFr6","solid_WFr7",
                                       "solid_AlFrShim","solid_WFrIn1","solid_WFrIn2","solid_WFrIn3","solid_AlChShim1","solid_AlChShim2","solid_AlChShim3","solid_AlChShim4",
	                                  "solid_BackW","solid_AlBkShim1","solid_AlBkShim2",
	                                  "solid_AlBkPlate"};

G4String EventAction::detName[]= {"detector_1","detector_2",
    "detector_3","detector_4","detector_5","detector_6","detector_7","detector_8", "detector_9", "Al","W","Be"};

G4String EventAction::BeName = "solid_BeWin";

std::fstream f;

//.................................................... 	  
 EventAction::EventAction(){
    for(int i=0; i<9 ; i++){
   	          SiSensID[i] = -1;
    }   
 }
//.................................................... 
 EventAction::~EventAction(){
 }

 void EventAction::BeginOfEventAction(const G4Event* evt){
//  G4cout << " <EventAction> BeginOfEventAction " << G4endl;
  //G4RunManager::GetRunManager()->rndmSaveThisEvent(); 

  // Store incident energy
  G4PrimaryParticle* primary = evt->GetPrimaryVertex()->GetPrimary();
  G4double incidentEnergy = primary->GetKineticEnergy();
  // Write incident energy
  f.open("EnergyDepositResult.txt", std::ios::app | std::ios::out);
  f << std::setw(9) << std::setprecision(6) << incidentEnergy/MeV;

  G4String collNam = "/SiDiscColl" ;   
  G4String collNm1 = "/shieldColl" ;  
  G4String collNm2 = "/BeDiscColl" ;                       	    
  for(int i=0; i<9; i++){
    SiSensID[i] = G4SDManager::GetSDMpointer()->GetCollectionID(bkSnam[i]+collNam);

  // 	G4cout << " SiSensID[i] " << SiSensID[i] << 
    //	 "bkSnam[i]+collNam "<< bkSnam[i]+collNam << G4endl;     
      }	 
    	 
//    	G4cout << " backGardID[i] " << backGardID[i] << 
//    	 "bkGnam[i]+collNm1 "<< bkGnam[i]+collNam << G4endl;     
          
   //for(int i=0; i<10 ; i++){ 
   	//shieldID[i] = G4SDManager::GetSDMpointer()->GetCollectionID(shlnam[i]+collNm1);
   	//G4cout << " shieldID[i] " << shieldID[i] << 
	//	"shlnam[i]+collNm1 "<< shlnam[i]+collNm1 << G4endl;
    //}
	
    
 //    BeDiscID = G4SDManager::GetSDMpointer()->GetCollectionID(BeName+collNm2);	 
     resetArrays();
  }
//**************************************************************
  void EventAction::EndOfEventAction(const G4Event* evt){

   evtNum = evt->GetEventID();
  
   // Store particle tracking
   G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
   numTrk = 0;
   if (trajectoryContainer) {
   	numTrk = trajectoryContainer->entries();
   }

//   G4cout << " <EventAction> printFront " << G4endl;         
   printFront(evt);
}
//*********************************************************
  void EventAction::printFront(const G4Event* evt){

// get the HCofThisEvent: it contains all the hits collections
//that have been defined (one hit collection may be associated to
//each detector).
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  HitsCollection* bkC = 0;
  G4String name;

//   G4cout << " <EventAction> printFront 01 " << G4endl;         
  
// k is number of time to  print the output data for each detector
// i is the detector number. RIC 4.5.08 I think...?
  f << "          ";
  for(size_t i=0; i<9 ; i++){
   G4double totE  = 0.;
   G4double totL  = 0.;
   G4double Etw = 0.;
   G4int    nStep = 0;
//   if (k==0 && HCE){
    	 bkC = (HitsCollection*)(HCE->GetHC(SiSensID[i]));
    	 name =  bkSnam[i];
//    }	 
//    if (k==1 && HCE) {
//    	bkC = (HitsCollection*)(HCE->GetHC(frontGardID[i]));
//    	 name =  ftGnam[i];
//    }	     	
    if (bkC){
       G4int nHit = bkC->entries();
//       G4cout << bkSnam[i] << " engy "  << (*bkC)[0]->GetEdep() << G4endl;
//       Etw = (*bkC)[0];
//       G4cout << (*bkC)[0];
       for (G4int ii=0;ii<nHit;ii++){
        totE += (*bkC)[ii]->GetEdep(); 
        totL += (*bkC)[ii]->GetTrak();
        nStep += (*bkC)[ii]->GetNStep();
      }
/*      G4cout << name 
           << " Edep & trk len (MeV,mm): " 
           << std::setw(10) << std::setprecision(6)            
           << totE/MeV  
           << std::setw(10) << std::setprecision(6) 
           << totL/mm << G4endl;   
*/
       f << std::setw(10) << std::setprecision(6)
           //    << Etw
           //    << std::setw(10) << std::setprecision(6)
           << totE / MeV;
   //    << std::setw(10) << std::setprecision(6)
   //    << totL/mm 
   //    << std::setw(10) << std::setprecision(6)
   //    << nStep
		   
    G4int detIndex =  GetDetNumber(name); 
    engyDep[detIndex] = totE/MeV;   
//    trkLeng[detIndex] = totL/mm;
    }
    
    }
   f << std::endl;
   f.close();
/*    f.open("EnergyThroughWindow.txt", std::ios::app | std::ios::out);
    f << " EthrWin (MeV)): \n";
    for (int i = 2; i = 2; i++)
    {
        G4double Ethrwin = 0.;

        bkC = (HitsCollection*)(HCE->GetHC("BeFrontDisc"));
        G4int nHit = bkC->entries();
        f << std::setw(10) << std::setprecision(6)
            <<nHit
    }
    f.close();
 */ }

   //------------------------------------------------------------------ 
void  EventAction::resetArrays(){ 
	
    for(int k=0; k<21 ; k++){ 
       engyDep[k] = 0.0;
       trkLeng[k] = 0.0;    
       numGamma[k] = 0;
       numElec[k] = 0;    
       numPosr[k] = 0;
       minKeGamm[k] = 0.0; 
       minKeElec[k] = 0.0; 
       minKePosr[k] = 0.0;
       maxKeGamm[k] = 0.0; 
       maxKeElec[k] = 0.0; 
       maxKePosr[k] = 0.0;       
	   } 
}    
//------------------------------------------------------------------ 
G4int EventAction::GetDetNumber(G4String thisDetName){
	  G4int ik = -99;
	  //G4cout << " REPTeventAction::GetDetNumber <enter> " << thisDetName << G4endl;
	  
	  for(int i=0;i<20;i++){ 	
	  	if(thisDetName.find(detName[i]) < 10){
	     ik = i;	
	  	 break;	
	  	}
	  }
	  if(thisDetName.find("Al") < 10)ik = 18;		
      if(thisDetName.find("W") < 10 )ik = 19;
      if(thisDetName.find("Be") < 10)ik = 20;		
	  
//    G4cout << " REPTeventAction::GetDetNumber <return> " 
//           << thisDetName << " "<< ik << G4endl;
	  
	  return ik;	  		  
}    
//------------------------------------------------------------------
