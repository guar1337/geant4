#include "EventAction.hh"
#include "cesiumHit.hh"
#include "siliconHit.hh"

#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"	
#include <G4SDManager.hh>

#include "g4root.hh"
#include "Randomize.hh"
#include <iomanip>


EventAction::EventAction():
	G4UserEventAction()
{
	tree=NULL;
	v2H=NULL;
	v6He=NULL;
	vBEAM=NULL;

}
EventAction::EventAction(TTree *T):
	G4UserEventAction(),
	fsiliconHCID(-1),
	fcesiumHCID(-1)
{
	
	tree=T;

	lv2H = new G4LorentzVector();
	lv6He = new G4LorentzVector();
	lvBEAM = new G4LorentzVector();
	lvDeutCM = new G4LorentzVector();
	tree->Bronch("vBEAM.","G4LorentzVector",&vBEAM);
	tree->Bronch("v2H.","G4LorentzVector",&v2H);
	tree->Bronch("v6He.","G4LorentzVector",&v6He);
	tree->Bronch("vDeutCM.","G4LorentzVector",&vDeutCM);

	//Deuterium part 
	//tree->Branch("thetaDeut",&thetaDeut,"thetaDeut/D");
	//tree->Branch("phiDeut",&phiDeut,"phiDeut/D");
	tree->Branch("Tdeut",&Tdeut,"Tdeut/D");
	tree->Branch("CsIdeut", CsIdeut, "CsIdeut[16]/D");
	tree->Branch("SideutX", SideutX, "SideutX[16]/D");
	tree->Branch("SideutY", SideutY, "SideutY[16]/D");
	tree->Branch("deutEDEP", &deutEDEP, "deutEDEP/D");

	//Helium part
	//tree->Branch("thetaHe",&thetaHe,"thetaHe/D");
	//tree->Branch("phiHe",&phiHe,"phiHe/D");
	tree->Branch("The",&The,"The/D");
	tree->Branch("CsIhe", CsIhe, "CsIhe[16]/D");
	tree->Branch("SiheY", SiheY, "SiheY[16]/D");
	tree->Branch("SiheX", SiheX, "SiheX[16]/D");
	tree->Branch("heEDEP", &heEDEP, "heEDEP/D");

	//BEAM
	tree->Branch("mazz",&mazz,"mazz/D");
	tree->Branch("Tbeam",&Tbeam,"Tbeam/D");
	tree->Branch("thetaCM",&thetaCM,"thetaCM/D");
	//tree->Branch("phiCM",&phiCM,"phiCM/D");	
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
	
	
	tree->Branch("exp1",&exp1,"exp1/d");

	//ReCo

}

EventAction::~EventAction()
{
	if(v2H) delete v2H;
	if(v6He) delete v6He;
	if(vBEAM) delete vBEAM;
	if(vDeutCM) delete vDeutCM;
}


void EventAction::BeginOfEventAction(const G4Event *event)
{
/*
	// initialisation per event
	auto sdManager = G4SDManager::GetSDMpointer();
	fsiliconHCID = sdManager->GetCollectionID("sensSilicon/siliconColl");
	fcesiumHCID = sdManager->GetCollectionID("sensCesium/cesiumColl");
*/
	//Beam
	ParticleInfo* particleInfo=(ParticleInfo*)p->GetUserInformation();
	G4LorentzVector V_BEAM(particleInfo->Get_LV_Beam());
	G4LorentzVector CM_2H(particleInfo->Get_LV_DeutCM());

	G4PrimaryParticle *p = event->GetPrimaryVertex(0)->GetPrimary(0);
	v2H = p->GetMomentumDirection();
	Tdeut = p->GetKineticEnergy()/MeV;
	thetaDeut = 180.*(p->GetMomentumDirection().getTheta())/double(CLHEP::pi);
	phiDeut = 180.*(p->GetMomentumDirection().getPhi())/double(CLHEP::pi);
	lv2H->setVectM(*v2H, p->GetMass());
	
	
	
	vDeutCM = p->GetMomentumDirection();
	Tbeam =V_BEAM.e();
	vBEAM->SetPxPyPzE(px,py,pz,Tbeam);

	double E_IN_CM_deut = IN_CM_deut.e();
	lvDeutCM->SetPxPyPzE(px,py,pz,E_IN_CM_deut);

			
	p= event->GetPrimaryVertex(0)->GetPrimary(1);	//scattered helium
	v6He = p->GetMomentumDirection();
	The = p->GetKineticEnergy();
	thetaHe = 180.*(p->GetMomentumDirection().getTheta())/double(CLHEP::pi);
	phiHe = 180.*(p->GetMomentumDirection().getPhi())/double(CLHEP::pi);
	lv6He->setVectM(*v6He, p->GetMass());
	mass6He = 5605.53497;

//G4cout<<E_IN_CM_deut-mass2H<<" CM= "<<Tbeam-mass6He<<" LAB= "<<G4endl;
	Xpos = event->GetPrimaryVertex(0)->GetX0()/mm;
	Ypos = event->GetPrimaryVertex(0)->GetY0()/mm;
	Zpos = event->GetPrimaryVertex(0)->GetZ0()/1000/mm;

}


void EventAction::EndOfEventAction(const G4Event* event)
{

//G4cout<<"I am starting new EVENT"<<G4endl;
G4HCofThisEvent* hce = event->GetHCofThisEvent();

if(!hce)
{
	G4cout<<"dupa"<<G4endl;
	return;
}

if (!hce) 
{
	G4ExceptionDescription msg1;
	msg1 << "Bida calkowita" << G4endl; 
	G4Exception("no cos nie poszlo",
	"ej,", JustWarning, msg1);
	return;
	}	//ending scope on if on hce
siliconHitsCollection * SiHC = NULL;
cesiumHitsCollection * CsIHC = NULL;
SiHC = static_cast<siliconHitsCollection*>(hce->GetHC(fsiliconHCID));
CsIHC = static_cast<cesiumHitsCollection*>(hce->GetHC(fcesiumHCID));

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
if((CsI_n_hit+Si_n_hit)>0)
	{
	tree->Fill();
	}
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