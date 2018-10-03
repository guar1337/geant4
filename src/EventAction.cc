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
	lv6He = new TLorentzVector();
	lvBeam = new TLorentzVector();
	lv2H_CM = new TLorentzVector();
	lv6He_CM = new TLorentzVector();
	
	tree->Bronch("lvBeam.",		"TLorentzVector",		&lvBeam);
	tree->Bronch("lv2H.",		"TLorentzVector", 	&lv2H);
	tree->Bronch("lv6He.",		"TLorentzVector", 	&lv6He);
	tree->Bronch("lv2H_CM.",	"TLorentzVector", 	&lv2H_CM);
	tree->Bronch("lv6He_CM.",	"TLorentzVector", 	&lv6He_CM);
	
	//Deuterium part 
	tree->Branch("Tdeut",		&Tdeut,"Tdeut/D");
	tree->Branch("CsIdeut", 	CsIdeut, "CsIdeut[16]/D");
	tree->Branch("SideutX", 	SideutX, "SideutX[16]/D");
	tree->Branch("SideutY",		SideutY, "SideutY[16]/D");
	tree->Branch("deutEDEP",	&deutEDEP, "deutEDEP/D");

	//Helium part

	tree->Branch("The",&The,"The/D");
	tree->Branch("CsIhe", CsIhe, "CsIhe[16]/D");
	tree->Branch("SiheY", SiheY, "SiheY[16]/D");
	tree->Branch("SiheX", SiheX, "SiheX[16]/D");
	tree->Branch("heEDEP", &heEDEP, "heEDEP/D");

	//BEAM
	tree->Branch("mazz",&mazz,"mazz/D");
	tree->Branch("Tbeam",&Tbeam,"Tbeam/D");
	tree->Branch("thetaCM",&thetaCM,"thetaCM/D");
	tree->Branch("phiCM",&phiCM,"phiCM/D");	
	tree->Branch("Xpos",&Xpos,"Xpos/D");
	tree->Branch("Ypos",&Ypos,"Ypos/D");
	tree->Branch("Zpos",&Zpos,"Zpos/D");

	tree->Branch("X6He",&X6He,"X6He/D");
	tree->Branch("Y6He",&Y6He,"Y6He/D");
	tree->Branch("Z6He",&Z6He,"Z6He/D");

	tree->Branch("X2H",&X2H,"X2H/D");
	tree->Branch("Y2H",&Y2H,"Y2H/D");
	tree->Branch("Z2H",&Z2H,"Z2H/D");

	//ReCo
	tree->Branch("reTheta6He", &reTheta6He, "reTheta6He/D");
	tree->Branch("reTheta2H", &reTheta2H, "reTheta2H/D");
	tree->Branch("labAngHe", &labAngHe, "labAngHe/D");
	tree->Branch("labAng2H", &labAng2H, "labAng2H/D");

	//ReCo

}

EventAction::~EventAction()
{
	if(v2H) delete v2H;
	if(v6He) delete v6He;
	if(lv2H) delete lv2H;
	if(lv6He_CM) delete lv6He_CM;
	if(lvBeam) delete lvBeam;
	if(lv2H_CM) delete lv2H_CM;

	if(tmp_lv6He_CM) delete tmp_lv6He_CM;
	if(tmp_lvBeam) delete tmp_lvBeam;
	if(tmp_lv2H_CM) delete tmp_lv2H_CM;
}


void EventAction::BeginOfEventAction(const G4Event *event)
{

	// initialisation per event
	auto sdManager = G4SDManager::GetSDMpointer();
	fsiliconHCID = sdManager->GetCollectionID("sensSilicon/siliconColl");
	fcesiumHCID = sdManager->GetCollectionID("sensCesium/cesiumColl");

	//2H scattered
	G4PrimaryParticle *PrimaryParticle_2H = event->GetPrimaryVertex(0)->GetPrimary(0);
	v2H->SetXYZ(PrimaryParticle_2H->GetPx(), PrimaryParticle_2H->GetPy(), PrimaryParticle_2H->GetPz());
	Tdeut = PrimaryParticle_2H->GetKineticEnergy()/MeV;
	thetaDeut = 180.*(PrimaryParticle_2H->GetMomentumDirection().getTheta())/double(CLHEP::pi);
	phiDeut = 180.*(PrimaryParticle_2H->GetMomentumDirection().getPhi())/double(CLHEP::pi);
	lv2H->SetVectM(*v2H, PrimaryParticle_2H->GetMass()+Tdeut);

	//6He scattered
	G4PrimaryParticle *PrimaryParticle_6He = event->GetPrimaryVertex(0)->GetPrimary(1);
	v6He->SetXYZ(PrimaryParticle_6He->GetPx(), PrimaryParticle_6He->GetPy(), PrimaryParticle_6He->GetPz());
	The = PrimaryParticle_6He->GetKineticEnergy();
	thetaHe = 180.*(PrimaryParticle_6He->GetMomentumDirection().getTheta())/double(CLHEP::pi);
	phiHe = 180.*(PrimaryParticle_6He->GetMomentumDirection().getPhi())/double(CLHEP::pi);
	lv6He->SetVectM(*v6He, PrimaryParticle_6He->GetMass()+The);

	//beam, deuterium and helium in CM
	
	ParticleInfo *particleInfo=(ParticleInfo*)PrimaryParticle_2H->GetUserInformation();
	tmp_lvBeam = new G4LorentzVector(particleInfo->Get_LV_Beam());
	lvBeam->SetPxPyPzE(tmp_lvBeam->px(), tmp_lvBeam->py(), tmp_lvBeam->pz(), tmp_lvBeam->m());
	tmp_lv2H_CM = new G4LorentzVector(particleInfo->Get_LV_Beam());
	lv2H_CM->SetPxPyPzE(tmp_lv2H_CM->px(), tmp_lv2H_CM->py(), tmp_lv2H_CM->pz(), tmp_lv2H_CM->m());
	tmp_lv6He_CM = new G4LorentzVector(particleInfo->Get_LV_Beam());
	lv6He_CM->SetPxPyPzE(tmp_lv6He_CM->px(), tmp_lv6He_CM->py(), tmp_lv6He_CM->pz(), tmp_lv6He_CM->m());

	
	Tbeam =lvBeam->E();

	//double E_IN_CM_deut = IN_CM_deut.e();

	Xpos = event->GetPrimaryVertex(0)->GetX0()/mm;
	Ypos = event->GetPrimaryVertex(0)->GetY0()/mm;
	Zpos = event->GetPrimaryVertex(0)->GetZ0()/1000/mm;

}


void EventAction::EndOfEventAction(const G4Event* event)
{

//G4cout<<"I am starting new EVENT"<<G4endl;
G4HCofThisEvent* HCofThisEvent = event->GetHCofThisEvent();

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
	}	//ending scope on if on HCofThisEvent

siliconHitsCollection *SiHC;
cesiumHitsCollection *CsIHC;
SiHC = static_cast<siliconHitsCollection*>(HCofThisEvent->GetHC(fsiliconHCID));
CsIHC = static_cast<cesiumHitsCollection*>(HCofThisEvent->GetHC(fcesiumHCID));

	if (!SiHC && !CsIHC) 
	{
	G4ExceptionDescription msg2;
	msg2 << "Bida czesciowa" << G4endl; 
	G4Exception("no cos nie poszlo",
	"ej,", JustWarning, msg2);
	return;


	}//ending scope on if on SiHC && CsIHC

	G4int SipixelNo=0;
	G4int SistripNo=0;
	G4int SidetectorNo=0;
	G4int CsIpixelNo=0;
	G4int CsIstripNo=0;
	G4int CsIdetectorNo=0;
	stuck=0;
auto Si_n_hit = SiHC->entries();
if (Si_n_hit>0)
	{
	stuck=((*SiHC)[0])->GetParticle()->GetPDGMass();
	if (stuck>1.0)
	{
	mazz=stuck;
	}
	else
	{
	mazz=0;
	}
	//G4cout<<stuck<<" uuuLALA"<<G4endl;
	for ( int iii = 0 ; iii < Si_n_hit; iii++)
	{
		//siliconHit *hit	= (siliconHit*) SiHC->GetHit(iii);
		SidetectorNo= ((*SiHC)[iii])->GetDetectorNo();


		if (SidetectorNo == 1)
		{	
		SipixelNo= 	((*SiHC)[iii])->GetPixelNo();
		SistripNo= 	((*SiHC)[iii])->GetStripNo();
		SiheX[SipixelNo]+=((*SiHC)[iii])->GetEnergy();
		SiheY[SistripNo]+=((*SiHC)[iii])->GetEnergy();
	G4ThreeVector positionAccu = ((*SiHC)[iii])->GetPos();
	//G4cout<<SipixelNo<<" JESTEM WOW, SUPER WOW "<<positionAccu<<G4endl;
		}
		else if(SidetectorNo == 0)
		{
		SipixelNo= 	((*SiHC)[iii])->GetPixelNo();
		SistripNo= 	((*SiHC)[iii])->GetStripNo();
		SideutX[SipixelNo]+=((*SiHC)[iii])->GetEnergy();
		SideutY[SistripNo]+=((*SiHC)[iii])->GetEnergy();
		}
	}
	
	}

auto CsI_n_hit = CsIHC->entries();
if (CsI_n_hit>0)
	{
	//G4cout<<"CsI_n_hit= "<<CsI_n_hit<<G4endl;
	for ( int iii = 0; iii < CsI_n_hit; iii++)
		{
		CsIdetectorNo= 	((*CsIHC)[iii])->GetDetectorNo();
		if(CsIdetectorNo == 1)
		{
		CsIpixelNo= 	((*CsIHC)[iii])->GetPixelNo();
		CsIstripNo= 	((*CsIHC)[iii])->GetStripNo();
		CsIhe[CsIpixelNo+4*CsIstripNo]+=((*CsIHC)[iii])->GetEnergy();
		}
		else if(CsIdetectorNo == 0)
		{
		CsIpixelNo= 	((*CsIHC)[iii])->GetPixelNo();
		CsIstripNo= 	((*CsIHC)[iii])->GetStripNo();
		CsIdeut[CsIpixelNo+4*CsIstripNo]+=((*CsIHC)[iii])->GetEnergy();
		}
	}
	}
//if((CsI_n_hit+Si_n_hit)>0)
//	{
	tree->Fill();
//	}
		//zeroing Silicon detectors
	SipixelNo=0;
	SistripNo=0;
	SidetectorNo=0;
	std::fill(SideutX,SideutX+16,0);
	std::fill(SideutY,SideutY+16,0);
	std::fill(SiheX,SiheX+16,0);
	std::fill(SiheY,SiheY+16,0);

	//zeroing CesiumIodide detectors
	CsIpixelNo=0;
	CsIstripNo=0;
	CsIdetectorNo=0;
	std::fill(CsIdeut,CsIdeut+16,0.0);
	std::fill(CsIhe,CsIhe+16,0.0);
	
}