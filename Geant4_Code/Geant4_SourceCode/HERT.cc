/* this is the first attempt to write a code for Geant4
Ian Crocker, 2 Mar 2008*/
//Revised by H. Zhao, 08/12/2020

#include "G4RunManager.hh"   //| 
#include "G4UImanager.hh"    //| General part of GEANT4 simulation, no need to change
#include "G4UIterminal.hh"   //|
#include "G4VisExecutive.hh" //| added 06/17/2020 to use /vis/
#include "G4UIExecutive.hh"  //| added 06/17/2020 to use interactive mode

#include "DetectorConstruction.hh"   //|
#include "DetectorMessenger.hh"   //|
#include "PhysicsList.hh"            //|
#include "PrimaryGeneratorAction.hh" //| head files for reptile
                                     //|
// 20 March 08  -->                  //|
#include "RunAction.hh"              //|
#include "TrackingAction.hh"         //|
#include "EventAction.hh"            //|
#include "StackingAction.hh"         //|
#include "SetSensDet.hh"             //|
// <-- 20 March 08  


int main(int argc,char** argv)
{

  // Detect interactive mode (if no arguments) and define UI session - added 06/17/2020
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  G4VUserDetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);

  G4VUserPhysicsList* physics = new PhysicsList;
  physics->SetVerboseLevel(1);
  runManager->SetUserInitialization(physics);

  // set mandatory user action class
  G4VUserPrimaryGeneratorAction* gen_action = new PrimaryGeneratorAction;
  runManager->SetUserAction(gen_action);

  // 20 March 08  -->	
  runManager->SetUserAction(new RunAction);
  runManager->SetUserAction(new EventAction);
  runManager->SetUserAction(new TrackingAction);
  runManager->SetUserAction(new StackingAction);
  // <-- 20 March 08  

  // Initialize G4 kernel
  runManager->Initialize();

  //added 06/17/2020 to use /vis/
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  if(argc>1) // execute an argument macro file if exist
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);    
  }
  else{ // start interactive session
    G4UIsession* session = new G4UIterminal();
    session->SessionStart();
    delete session;
  }
 
	// Job termination
	//
	// Free the store: user actions, physics_list and detector_description are
	//                 owned and deleted by the run manager, so they should not
	//                 be deleted in the main() program !
	//
	
        delete visManager; //|added 06/17/2020 to use /vis/
        delete runManager;
        
	return 0;
}
