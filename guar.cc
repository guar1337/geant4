#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "QBBC.hh"
#include "Randomize.hh"
#include "G4Timer.hh"

int main(int argc,char** argv)
{
	G4Timer timer;
	timer.Start();
	//root output file
	TFile *outFile=NULL;	
	TTree *outTree=NULL;
	outFile = new TFile("gurney_bis.root","RECREATE");
	outTree=new TTree("simevents","MC events");

	// Detect interactive mode (if no arguments) and define UI session
	//
	G4UIExecutive* ui = 0;
	if ( argc == 1 ) {
	ui = new G4UIExecutive(argc, argv);
	}
	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	
	// Construct the default run manager
	//
	G4RunManager *runManager = new G4RunManager;
	// Set mandatory initialization classes
	runManager->SetUserInitialization(new DetectorConstruction());
	runManager->SetUserInitialization(new QBBC);	
	runManager->SetUserAction(new PrimaryGeneratorAction());
	runManager->SetUserAction(new EventAction(outTree));	
	runManager->Initialize();
		
	// Initialize visualization
	//
	G4VisManager *visManager = new G4VisExecutive;
	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();
	G4UImanager *UImanager = G4UImanager::GetUIpointer();

	// Process macro or start UI session
	//
	if ( ! ui ) 
	{ 
		// batch mode
		G4String command = "/control/execute init_vis.mac";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command+fileName);
	}
	else
	{ 
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		delete ui;
	}

	// Job termination
	// Free the store: user actions, physics_list and detector_description are
	// owned and deleted by the run manager, so they should not be deleted 
	// in the main() program !

	delete visManager;
	delete runManager;
	timer.Stop();
	printf("Zajelo to %f sekund\n",timer.GetRealElapsed());
	outFile->cd();
	outTree->Write();
	outFile->Close();
	
return 0;
}