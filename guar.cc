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
	//G4Timer timer;
	//timer.Start();
	//root output file
	TFile *outFile=NULL;	
	TTree *outTree=NULL;
	std::string geoID = std::to_string(cs::runNo/10);
	std::string targetMass = std::to_string(cs::tarMass);
	TString fName = "/home/zalewski/dataTmp/MC/v3/mc_geo" + geoID + "_" + targetMass + "H_dt.root";

	outFile = new TFile(fName.Data(),"RECREATE");
	outTree=new TTree("calibrated","MC events");

	
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
	G4UIExecutive* ui = 0;
	G4VisManager *visManager = new G4VisExecutive;
	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();
	G4UImanager *UImanager = G4UImanager::GetUIpointer();
	printf("\n\nAAA\n\n%s\t%d\n\nAAA\n\n",argv[1], argc);
	
	if ( argc == 2)
	{	
		printf("argc = %d\n", argc);
		//strcmp gives 0 for equal strings
		if (strcmp(argv[1],"init_vis.mac")==0)
		{
			
			ui = new G4UIExecutive(argc, argv);
			UImanager->ApplyCommand("/control/execute init_vis.mac");
			ui->SessionStart();
			delete ui;
		}

		else
		{
			G4String cmd = "/control/execute ";
			G4String macroName = argv[1];
			UImanager->ApplyCommand(cmd + macroName);
		}

	}

	// Job termination
	// Free the store: user actions, physics_list and detector_description are
	// owned and deleted by the run manager, so they should not be deleted 
	// in the main() program !

	delete visManager;
	delete runManager;
	//timer.Stop();
	//printf("Zajelo to %f sekund\n",timer.GetRealElapsed());
	outFile->cd();
	outTree->Write();
	outFile->Close();
	
return 0;
}