// By LYK (2/14/2019)
// Goal: To allow the update of geometry of the window (both material and depth) using messenger (or commands)
// Steps:
// 1. 

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include <sstream>
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//DetectorMessenger::DetectorMessenger()
DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det)
{
  
    fDirectory = new G4UIdirectory("/REPTileGeant_UserMessenger/");
  fDirectory->SetGuidance("UI commands specific to this example.");

  fDetDirectory = new G4UIdirectory("/REPTileGeant_UserMessenger/det/");
  fDetDirectory->SetGuidance("Detector construction control");

  fWindMatCmd = new G4UIcmdWithAString("/REPTile2/det/setWindowMaterial",this);
  fWindMatCmd->SetGuidance("Select Material of the Window.: Al or Be");
  fWindMatCmd->SetParameterName("materialChoice",false);
  fWindMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fWindDepCmd = new G4UIcmdWithADoubleAndUnit("/REPTile2/det/setWindowDepth",this);
  fWindDepCmd->SetGuidance("Define the window depth: (in mm)");
  fWindDepCmd->SetParameterName("windowLengthChoice",false);
  fWindDepCmd->SetRange("windowLengthChoice>0.");
  fWindDepCmd->SetUnitCategory("Length");
  fWindDepCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fWindMatCmd;
  delete fWindDepCmd;
  delete fDirectory;
  delete fDetDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fWindMatCmd )
   { fDetector->SetWindowMaterial(newValue);}

  if( command == fWindDepCmd ) {
    fDetector->SetWindowDepth(fWindDepCmd->GetNewDoubleValue(newValue));
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
