//file:///home/zalewski/aku/geant4/src/EventAction.cc
#include "EventAction.hh"
#include "g4root.hh"
#include "Randomize.hh"
#include <iomanip>
#include "cesiumHit.hh"
#include "siliconHit.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"

#include <sys/ioctl.h> // For ioctl, TIOCGWINSZ
#include <unistd.h> // For STDOUT_FILENO

EventAction::EventAction():
	G4UserEventAction()
{
	tree=NULL;
	lv2H=NULL;
	lv6He=NULL;
	lvBeam=NULL;
}
EventAction::EventAction(TTree *T):
	G4UserEventAction(),
	fsiliconHCID(-1),
	fcesiumHCID(-1)
{
	tree=T;
	//Work on getting TLorentVector to the outTree
	//Success? Getting TLorentzVector out of the tree is possible, but I still get warnings from TTree
	v2H = new TVector3();
	v6He = new TVector3();
	lv2H = new TLorentzVector();
	lv5He = new TLorentzVector();
	lv6He = new TLorentzVector();
	lvBeam = new TLorentzVector();
	lv2H_CM = new TLorentzVector();
	lv6He_CM = new TLorentzVector();

	struct winsize size;
	ioctl(STDOUT_FILENO,TIOCGWINSZ,&size);
	consoleWidth = size.ws_col-8;
	
	tree->Bronch("flvBeam.",		"TLorentzVector",		&lvBeam);
	tree->Bronch("lv2H.",		"TLorentzVector", 	&lv2H);
	tree->Bronch("lv5He.",		"TLorentzVector", 	&lv5He);
	tree->Bronch("lv6He.",		"TLorentzVector", 	&lv6He);
	tree->Bronch("lv2H_CM.",	"TLorentzVector", 	&lv2H_CM);
	tree->Bronch("lv6He_CM.",	"TLorentzVector", 	&lv6He_CM);
	
	//Deuterium part 
	tree->Branch("cal_CsI_L", 	cal_CsI_L, 	"cal_CsI_L[16]/D");
	tree->Branch("cal_SQX_L", 	cal_SQX_L, 	"cal_SQX_L[32]/D");
	tree->Branch("cal_SQY_L",	cal_SQY_L, 	"cal_SQY_L[16]/D");
	tree->Branch("fsqlang",		&sqlang, 	"fsqlang/D");
	tree->Branch("fsqrang5He",	&sqrang5He, 	"fsqrang5He/D");
	tree->Branch("fsqlde",		&sqlde,		"fsqlde/D");
	tree->Branch("fsqletot",	&sqletot,	"fsqletot/D");

	//Helium part

	tree->Branch("cal_CsI_R",		cal_CsI_R,		"cal_CsI_R[16]/D");
	tree->Branch("cal_SQY_R",		cal_SQY_R,		"cal_SQY_R[16]/D");
	tree->Branch("cal_SQX_R",		cal_SQX_R,		"cal_SQX_R[32]/D");
	tree->Branch("fsqrang",		&sqrang,		"fsqrang/D");
	tree->Branch("fsqrde",		&sqrde,		"fsqrde/D");
	tree->Branch("fsqretot",	&sqretot,	"fsqretot/D");

	//BEAM
	tree->Branch("fbeamT",	&beamT,	"fbeamT/D");
	tree->Branch("tof",		&tof,	"tof/D");
	tree->Branch("fThetaCM",		&thetaCM,	"fThetaCM/D");
	tree->Branch("fevx",		&evx,		"fevx/D");
	tree->Branch("fevy",		&evy,		"fevy/D");
	tree->Branch("fevz",		&evz,		"fevz/D");

	tree->Branch("fX6He",	&X6He,	"fX6He/D");
	tree->Branch("fY6He",	&Y6He,	"fY6He/D");
	tree->Branch("fZ6He",	&Z6He,	"fZ6He/D");
	tree->Branch("fSQX_R_strip",	&fSQX_R_strip,	"fSQX_R_strip/I");


	tree->Branch("fX2H",		&X2H,		"fX2H/D");
	tree->Branch("fY2H",		&Y2H,		"fY2H/D");
	tree->Branch("fZ2H",		&Z2H,		"fZ2H/D");
	tree->Branch("fSQX_L_strip",	&fSQX_L_strip,	"fSQX_L_strip/I");

	tree->Branch("fEve2H",	&fEve2H,	"fEve2H/B");
	tree->Branch("fEve6He",	&fEve6He,"fEve6He/B");

	tree->Branch("MWPC_1_X", &MWPC_1_X,	 "MWPC_1_X/D");
	tree->Branch("MWPC_2_X", &MWPC_2_X,	 "MWPC_2_X/D");
	tree->Branch("MWPC_1_Y", &MWPC_1_Y,	 "MWPC_1_Y/D");
	tree->Branch("MWPC_2_Y", &MWPC_2_Y,	 "MWPC_2_Y/D");

	tree->Branch("nx1", &nx1,	 "nx1/D");
	tree->Branch("nx2", &nx2,	 "nx2/D");
	tree->Branch("ny1", &ny1,	 "ny1/D");
	tree->Branch("ny2", &ny2,	 "ny2/D");

	//create ghost columns for compatibility with experimental data
	tree->Branch("cal_SQ300", cal_SQ300,	 "cal_SQ300[16]/D");
	tree->Branch("tSQX_L", 	tSQX_L,	 "tSQX_L[16]/s");
	tree->Branch("tSQX_R", 	tSQX_R,	 "tSQX_R[16]/s");
	tree->Branch("geo", 		&geo,	 "geo/I");

}

EventAction::~EventAction()
{
	if(v2H) delete v2H;
	if(v6He) delete v6He;
	if(lv2H) delete lv2H;
	if(lv5He) delete lv5He;
	if(lv6He) delete lv6He;
	if(lv6He_CM) delete lv6He_CM;
	if(lvBeam) delete lvBeam;
	if(lv2H_CM) delete lv2H_CM;
}

void EventAction::BeginOfEventAction(const G4Event *event)
{
	Long_t nEntries = G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed();
/*
	std::cout << "[";
	int pos = consoleWidth * progress;
	for (int i = 0; i < consoleWidth; ++i)
	{
		if (i < pos) std::cout << "#";
		else if (i == pos) std::cout << ">";
		else std::cout << ".";
	}
	std::cout << "] " << int(progress * 100.0) << " %\r";
	std::cout.flush();

	if( event->GetEventID() % ( nEntries / 100 ) == 0) progress += 0.01;
*/	// initialisation per event
	auto sdManager = G4SDManager::GetSDMpointer();
	fsiliconHCID = sdManager->GetCollectionID("siliconColl");
	fcesiumHCID = sdManager->GetCollectionID("cesiumColl");

	//2H scattered
	G4PrimaryParticle *PrimaryParticle_2H = event->GetPrimaryVertex(0)->GetPrimary(0);
	v2H->SetXYZ(PrimaryParticle_2H->GetPx(), PrimaryParticle_2H->GetPy(), PrimaryParticle_2H->GetPz());
	lv2H->SetVectM(*v2H, PrimaryParticle_2H->GetMass());

	//6He scattered
	G4PrimaryParticle *PrimaryParticle_6He = event->GetPrimaryVertex(0)->GetPrimary(1);
	v6He->SetXYZ(PrimaryParticle_6He->GetPx(), PrimaryParticle_6He->GetPy(), PrimaryParticle_6He->GetPz());
	lv6He->SetVectM(*v6He, PrimaryParticle_6He->GetMass());

	//beam, (deuterium and helium in CM)
		
	ParticleInfo *particleInfo=(ParticleInfo*)PrimaryParticle_2H->GetUserInformation();

	tmp_lvBeam = new G4LorentzVector(particleInfo->Get_LV_Beam());	
	lvBeam->SetPxPyPzE(tmp_lvBeam->px(), tmp_lvBeam->py(), tmp_lvBeam->pz(), tmp_lvBeam->e());

	tmp_lv2H_CM = new G4LorentzVector(particleInfo->Get_LV_2H_CM());
	lv2H_CM->SetPxPyPzE(tmp_lv2H_CM->px(), tmp_lv2H_CM->py(), tmp_lv2H_CM->pz(), tmp_lv2H_CM->e());
	
	phiCM = tmp_lv2H_CM->phi();

	tmp_lv5He = new G4LorentzVector(particleInfo->Get_LV_5He());
	lv5He->SetPxPyPzE(tmp_lv5He->px(), tmp_lv5He->py(), tmp_lv5He->pz(), tmp_lv5He->e());

	tmp_lv6He_CM = new G4LorentzVector(particleInfo->Get_LV_6He_CM());
	lv6He_CM->SetPxPyPzE(tmp_lv6He_CM->px(), tmp_lv6He_CM->py(), tmp_lv6He_CM->pz(), tmp_lv6He_CM->e());
	thetaCM = particleInfo->Get_thetaCM();

	TLorentzVector m_lvTar(0.0,0.0,0.0,cs::mass2H);
	TLorentzVector m_lvCM = *lvBeam + m_lvTar;
	corrthetaCM = lv6He_CM->Vect().Angle(m_lvCM.Vect());

	MWPC_1_X = particleInfo->Get_MWPC_1_X();
	MWPC_2_X = particleInfo->Get_MWPC_2_X();
	MWPC_1_Y = particleInfo->Get_MWPC_1_Y();
	MWPC_2_Y = particleInfo->Get_MWPC_2_Y();

	nx1 = particleInfo->Get_nx1();
	nx2 = particleInfo->Get_nx2();
	ny1 = particleInfo->Get_ny1();
	ny2 = particleInfo->Get_ny2();

	sqlang = lvBeam->Angle(*v2H) * TMath::RadToDeg();
	sqrang5He = lvBeam->Angle(lv5He->Vect()) * TMath::RadToDeg();
	sqrang = lvBeam->Angle(*v6He) * TMath::RadToDeg();

	beamT = lvBeam->E()-lvBeam->M();
	tof = (cs::tofBase)/(lvBeam->Beta()*cs::c);

	//double E_IN_CM_deut = IN_CM_deut.e();

	evx = event->GetPrimaryVertex(0)->GetX0()*mm;
	evy = event->GetPrimaryVertex(0)->GetY0()*mm;
	evz = event->GetPrimaryVertex(0)->GetZ0()*mm;

	
	//printf("evX: %f\tevY: %f\tevZ: %f\n\n", evx, evy, evz);
}


void EventAction::EndOfEventAction(const G4Event *event)
{
	//G4cout<<"I am starting new EVENT"<<G4endl;
	G4HCofThisEvent *HCofThisEvent = event->GetHCofThisEvent();

	if(!HCofThisEvent)
	{
		printf("Didn't find HitsCollection!\n");
		return;
	}

	if (!HCofThisEvent) 
	{
		G4ExceptionDescription msg1;
		msg1 << "Bida calkowita" << G4endl; 
		G4Exception("no cos nie poszlo",
		"ej,", JustWarning, msg1);
		return;
	}

	siliconHitsCollection *SiHC = nullptr;
	cesiumHitsCollection *CsIHC = nullptr;
	SiHC = static_cast<siliconHitsCollection*>(HCofThisEvent->GetHC(fsiliconHCID));
	CsIHC = static_cast<cesiumHitsCollection*>(HCofThisEvent->GetHC(fcesiumHCID));

	if (fsiliconHCID<0 || fcesiumHCID<0) 
	{
		G4ExceptionDescription msg2;
		msg2 << "Nie ma Hit collection w EventAction" << G4endl; 
		G4Exception("no cos nie poszlo",
		"ej,", JustWarning, msg2);
		return;
	}

	SipixelNo=0;
	SistripNo=0;
	SidetectorNo=0;
	CsIpixelNo=0;
	CsIstripNo=0;
	CsIdetectorNo=0;
	sqrde=0;
	sqlde=0;
	sqretot=0;
	sqletot=0;
	fSQX_L_strip=0;
	fSQX_R_strip=0;

		//zeroing Silicon detectors
	std::fill(cal_SQX_L,cal_SQX_L+32,0.0);
	std::fill(cal_SQY_L,cal_SQY_L+16,0.0);
	std::fill(cal_SQX_R,cal_SQX_R+32,0.0);
	std::fill(cal_SQY_R,cal_SQY_R+16,0.0);
	std::fill(cal_SQ300,cal_SQ300+16,0.0);
	//zeroing CesiumIodide detectors

	std::fill(cal_CsI_L,cal_CsI_L+16,0.0);
	std::fill(cal_CsI_R,cal_CsI_R+16,0.0);

	fEve2H = false;
	fEve6He = false;

	int Si_n_hit = SiHC->entries();
	if (Si_n_hit>0)
	{
		for ( int iii = 0 ; iii < Si_n_hit; iii++)
		{
			//siliconHit *hit	= (siliconHit*) SiHC->GetHit(iii);
			SidetectorNo= ((*SiHC)[iii])->GetDetectorNo();
			SipixelNo= 	((*SiHC)[iii])->GetPixelNo();
			SistripNo= 	((*SiHC)[iii])->GetStripNo();
			G4ThreeVector positionAccu = ((*SiHC)[iii])->GetPos();

			if (SidetectorNo == 1)
			{
				if (fEve6He == false)
				{
					X6He = ((*SiHC)[iii])->GetPos().x();
					Y6He = ((*SiHC)[iii])->GetPos().y();
					Z6He = ((*SiHC)[iii])->GetPos().z();
					fEve6He = true;
				} 
				
				cal_SQX_R[SipixelNo]+=((*SiHC)[iii])->GetEnergy();
				cal_SQY_R[SistripNo]+=((*SiHC)[iii])->GetEnergy();
				//printf("SipixelNo: %d\tSistripNo: %d\tEnergy: %f\n", SipixelNo, SistripNo, ((*SiHC)[iii])->GetEnergy());
				sqrde+=((*SiHC)[iii])->GetEnergy();
				fSQX_R_strip = SipixelNo;
			}
				
			else if (SidetectorNo == 0)
			{
				if (fEve2H == false)
				{
					X2H = ((*SiHC)[iii])->GetPos().x();
					Y2H = ((*SiHC)[iii])->GetPos().y();
					Z2H = ((*SiHC)[iii])->GetPos().z();
					fEve2H = true;
				} 

				cal_SQX_L[SipixelNo]+=((*SiHC)[iii])->GetEnergy();
				cal_SQY_L[SistripNo]+=((*SiHC)[iii])->GetEnergy();
				sqlde+=((*SiHC)[iii])->GetEnergy();
				fSQX_L_strip = SipixelNo;
			}
		}	
	}

	int CsI_n_hit = CsIHC->entries();
	if (CsI_n_hit>0)
	{
		for ( int iii = 0; iii < CsI_n_hit; iii++)
		{
			CsIdetectorNo= 	((*CsIHC)[iii])->GetDetectorNo();
			CsIpixelNo= 	((*CsIHC)[iii])->GetPixelNo();
			CsIstripNo= 	((*CsIHC)[iii])->GetStripNo();
			//printf("DetNo: %d\tEdep: %f\thitNo: %d\n", CsIdetectorNo,((*CsIHC)[iii])->GetEnergy(), iii);

			if(CsIdetectorNo == 0)
			{
				cal_CsI_L[CsIpixelNo+4*CsIstripNo]+=((*CsIHC)[iii])->GetEnergy();
				sqletot+=((*CsIHC)[iii])->GetEnergy();
				//printf("Energia w de:\t%f\tentries: %d\n", cal_CsI_R[CsIpixelNo+4*CsIstripNo], CsI_n_hit);
			}

			else if(CsIdetectorNo == 1)
			{
				cal_CsI_R[CsIpixelNo+4*CsIstripNo]+=((*CsIHC)[iii])->GetEnergy();
				sqretot+=((*CsIHC)[iii])->GetEnergy();
				
				//if (cal_CsI_R[CsIpixelNo+4*CsIstripNo]<0.1) printf("Energia w he:\t%f\n", cal_CsI_R[CsIpixelNo+4*CsIstripNo]);
			}
		}
	}

	geo = cs::runNo/10;
	if (sqlde>0)
	{
		
	}
	
tree->Fill();

	if(tmp_lv5He) delete tmp_lv5He;
	if(tmp_lv6He_CM) delete tmp_lv6He_CM;
	if(tmp_lvBeam) delete tmp_lvBeam;
	if(tmp_lv2H_CM) delete tmp_lv2H_CM;

}
