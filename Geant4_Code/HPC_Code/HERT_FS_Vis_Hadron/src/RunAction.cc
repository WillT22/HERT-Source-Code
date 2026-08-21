//
//
//

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

std::fstream f_run;

RunAction::RunAction()
: G4UserRunAction(), fHistoManager(0)
{
    fHistoManager = new HistoManager();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
    delete fHistoManager;
}



void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    if ( analysisManager->IsActive() ) {
        analysisManager->OpenFile();
    }

    std::remove("EnergyDepositResult.txt");
    f_run.open("EnergyDepositResult.txt", std::ios::app | std::ios::out);
    f_run << "Einc(MeV) " << "Edep(MeV) "
        << "Detector1 " << "Detector2 " << "Detector3 " << "Detector4 " << "Detector5 " << "Detector6 " << "Detector7 " << "Detector8 " << "Detector9" << std::endl;
    f_run.close();
}


void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    if ( analysisManager->IsActive() ) {
        analysisManager->Write();
        analysisManager->CloseFile();
    }

  //compute statistics: mean and rms
  //

  //print
  //
 
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << G4endl;
     
}
