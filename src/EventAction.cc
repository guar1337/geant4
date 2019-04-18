#include "EventAction.hh"
<<<<<<< HEAD
#include "g4root.hh"
#include "Randomize.hh"
#include <iomanip>
#include "cesiumHit.hh"
#include "siliconHit.hh"
=======
#include "cesiumHit.hh"
#include "siliconHit.hh"

>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4HCofThisEvent.hh"
<<<<<<< HEAD
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
=======
#include "G4VHitsCollection.hh"	
#include <G4SDManager.hh>

#include "g4root.hh"
#include "Randomize.hh"
#include <iomanip>
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217


EventAction::EventAction():
	G4UserEventAction()
{
	tree=NULL;
<<<<<<< HEAD
	lv2H=NULL;
	lv6He=NULL;
	lvBeam=NULL;
=======
	v2H=NULL;
	v6He=NULL;
	vBEAM=NULL;

>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
}
EventAction::EventAction(TTree *T):
	G4UserEventAction(),
	fsiliconHCID(-1),
	fcesiumHCID(-1)
{
<<<<<<< HEAD
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
	tree->Branch("sqlesum",		&sqlesum,	"sqlesum/D");
	tree->Branch("CsIdeut", 	CsIdeut, 	"CsIdeut[16]/D");
	tree->Branch("SideutX", 	SideutX, 	"SideutX[16]/D");
	tree->Branch("SideutY",		SideutY, 	"SideutY[16]/D");
	tree->Branch("sqlang",		&sqlang, 	"sqlang/D");
	tree->Branch("sqlde",		&sqlde,		"sqlde/D");
	tree->Branch("sqletot",		&sqletot,	"sqletot/D");
	tree->Branch("sqlphi",		&sqlphi,		"sqlphi/D");
	tree->Branch("sqltheta",	&sqltheta,	"sqltheta/D");


	//Helium part

	tree->Branch("sqresum",		&sqresum,	"sqresum/D");
	tree->Branch("CsIhe",		CsIhe,		"CsIhe[16]/D");
	tree->Branch("SiheY",		SiheY,		"SiheY[16]/D");
	tree->Branch("SiheX",		SiheX,		"SiheX[16]/D");
	tree->Branch("sqrang",		&sqrang,		"sqrang/D");
	tree->Branch("sqrde",		&sqrde,		"sqrde/D");
	tree->Branch("sqretot",		&sqretot,	"sqretot/D");
	tree->Branch("sqrphi",		&sqrphi,		"sqrphi/D");
	tree->Branch("sqrtheta",	&sqrtheta,	"sqrtheta/D");
	tree->Branch("t_sqrang",		&t_sqrang, 	"t_sqrang/D");

	//BEAM
	tree->Branch("beamT",&beamT,"TbebeamTam/D");
	tree->Branch("thetaCM",&thetaCM,"thetaCM/D");
	tree->Branch("phiCM",&phiCM,"phiCM/D");
=======
	
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
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
	tree->Branch("Xpos",&Xpos,"Xpos/D");
	tree->Branch("Ypos",&Ypos,"Ypos/D");
	tree->Branch("Zpos",&Zpos,"Zpos/D");

	tree->Branch("X6He",&X6He,"X6He/D");
	tree->Branch("Y6He",&Y6He,"Y6He/D");
	tree->Branch("Z6He",&Z6He,"Z6He/D");

	tree->Branch("X2H",&X2H,"X2H/D");
	tree->Branch("Y2H",&Y2H,"Y2H/D");
	tree->Branch("Z2H",&Z2H,"Z2H/D");
<<<<<<< HEAD
=======

	//ReCo
	tree->Branch("reTheta6He", &reTheta6He, "reTheta6He/D");
	tree->Branch("reTheta2H", &reTheta2H, "reTheta2H/D");
	tree->Branch("labAngHe", &labAngHe, "labAngHe/D");
	tree->Branch("labAng2H", &labAng2H, "labAng2H/D");
	
	
	tree->Branch("exp1",&exp1,"exp1/d");

	//ReCo

>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
}

EventAction::~EventAction()
{
	if(v2H) delete v2H;
	if(v6He) delete v6He;
<<<<<<< HEAD
	if(lv2H) delete lv2H;
	if(lv6He_CM) delete lv6He_CM;
	if(lvBeam) delete lvBeam;
	if(lv2H_CM) delete lv2H_CM;

	if(tmp_lv6He_CM) delete tmp_lv6He_CM;
	if(tmp_lvBeam) delete tmp_lvBeam;
	if(tmp_lv2H_CM) delete tmp_lv2H_CM;
=======
	if(vBEAM) delete vBEAM;
	if(vDeutCM) delete vDeutCM;
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
}


void EventAction::BeginOfEventAction(const G4Event *event)
{
<<<<<<< HEAD

	// initialisation per event
	auto sdManager = G4SDManager::GetSDMpointer();
	fsiliconHCID = sdManager->GetCollectionID("siliconColl");
	fcesiumHCID = sdManager->GetCollectionID("cesiumColl");

	//2H scattered
	G4PrimaryParticle *PrimaryParticle_2H = event->GetPrimaryVertex(0)->GetPrimary(0);
	v2H->SetXYZ(PrimaryParticle_2H->GetPx(), PrimaryParticle_2H->GetPy(), PrimaryParticle_2H->GetPz());
	sqlesum = PrimaryParticle_2H->GetKineticEnergy()/MeV;
	sqltheta = 180.*(PrimaryParticle_2H->GetMomentumDirection().getTheta())/double(CLHEP::pi);
	sqlphi = 180.*(PrimaryParticle_2H->GetMomentumDirection().getPhi())/double(CLHEP::pi);
	lv2H->SetVectM(*v2H, PrimaryParticle_2H->GetMass()+sqlesum);

	//6He scattered
	G4PrimaryParticle *PrimaryParticle_6He = event->GetPrimaryVertex(0)->GetPrimary(1);
	v6He->SetXYZ(PrimaryParticle_6He->GetPx(), PrimaryParticle_6He->GetPy(), PrimaryParticle_6He->GetPz());
	sqresum = PrimaryParticle_6He->GetKineticEnergy();
	sqrtheta = 180.*(PrimaryParticle_6He->GetMomentumDirection().getTheta())/double(CLHEP::pi);
	sqrphi = 180.*(PrimaryParticle_6He->GetMomentumDirection().getPhi())/double(CLHEP::pi);
	lv6He->SetVectM(*v6He, PrimaryParticle_6He->GetMass()+sqresum);

	//beam, (deuterium and helium in CM)
	
	ParticleInfo *particleInfo=(ParticleInfo*)PrimaryParticle_2H->GetUserInformation();

	tmp_lvBeam = new G4LorentzVector(particleInfo->Get_LV_Beam());	
	lvBeam->SetPxPyPzE(tmp_lvBeam->px(), tmp_lvBeam->py(), tmp_lvBeam->pz(), tmp_lvBeam->e());

	tmp_lv2H_CM = new G4LorentzVector(particleInfo->Get_LV_2H_CM());
	lv2H_CM->SetPxPyPzE(tmp_lv2H_CM->px(), tmp_lv2H_CM->py(), tmp_lv2H_CM->pz(), tmp_lv2H_CM->e());

	tmp_lv6He_CM = new G4LorentzVector(particleInfo->Get_LV_6He_CM());
	lv6He_CM->SetPxPyPzE(tmp_lv6He_CM->px(), tmp_lv6He_CM->py(), tmp_lv6He_CM->pz(), tmp_lv6He_CM->e());

	sqlang = 180.0 * (lvBeam->Angle(*v2H))/double(CLHEP::pi);
	sqrang = 180.0 * (lvBeam->Angle(*v6He))/double(CLHEP::pi);

	beamT =lvBeam->E()-lvBeam->M();

	//double E_IN_CM_deut = IN_CM_deut.e();

=======
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
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
	Xpos = event->GetPrimaryVertex(0)->GetX0()/mm;
	Ypos = event->GetPrimaryVertex(0)->GetY0()/mm;
	Zpos = event->GetPrimaryVertex(0)->GetZ0()/1000/mm;

<<<<<<< HEAD
	t_sqrang = -180.0*atan(-sin(2*lvBeam->Angle(*v2H))/(cos(2*lvBeam->Angle(*v2H))+6.01888589/2.01410177))/double(CLHEP::pi);

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
=======
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
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
{
	G4ExceptionDescription msg1;
	msg1 << "Bida calkowita" << G4endl; 
	G4Exception("no cos nie poszlo",
	"ej,", JustWarning, msg1);
	return;
<<<<<<< HEAD
	}	//ending scope on if on HCofThisEvent

siliconHitsCollection *SiHC = nullptr;
cesiumHitsCollection *CsIHC = nullptr;
SiHC = static_cast<siliconHitsCollection*>(HCofThisEvent->GetHC(fsiliconHCID));
CsIHC = static_cast<cesiumHitsCollection*>(HCofThisEvent->GetHC(fcesiumHCID));

if (fsiliconHCID<0 || fcesiumHCID<0) 
{
=======
	}	//ending scope on if on hce
siliconHitsCollection * SiHC = NULL;
cesiumHitsCollection * CsIHC = NULL;
SiHC = static_cast<siliconHitsCollection*>(hce->GetHC(fsiliconHCID));
CsIHC = static_cast<cesiumHitsCollection*>(hce->GetHC(fcesiumHCID));

	if (!SiHC && !CsIHC) 
	{
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
	G4ExceptionDescription msg2;
	msg2 << "Bida czesciowa" << G4endl; 
	G4Exception("no cos nie poszlo",
	"ej,", JustWarning, msg2);
	return;
<<<<<<< HEAD
}//ending scope on if on SiHC && CsIHC

G4int SipixelNo=0;
G4int SistripNo=0;
G4int SidetectorNo=0;
G4int CsIpixelNo=0;
G4int CsIstripNo=0;
G4int CsIdetectorNo=0;



auto Si_n_hit = SiHC->entries();
if (Si_n_hit>0)
{
=======


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
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
	for ( int iii = 0 ; iii < Si_n_hit; iii++)
	{
		//siliconHit *hit	= (siliconHit*) SiHC->GetHit(iii);
		SidetectorNo= ((*SiHC)[iii])->GetDetectorNo();
<<<<<<< HEAD
		SipixelNo= 	((*SiHC)[iii])->GetPixelNo();
		SistripNo= 	((*SiHC)[iii])->GetStripNo();
		G4ThreeVector positionAccu = ((*SiHC)[iii])->GetPos();

		if (SidetectorNo == 1)
		{	
			SiheX[SipixelNo]+=((*SiHC)[iii])->GetEnergy();
			SiheY[SistripNo]+=((*SiHC)[iii])->GetEnergy();			
		}
			
		else if(SidetectorNo == 0)
		{
			SideutX[SipixelNo]+=((*SiHC)[iii])->GetEnergy();
			SideutY[SistripNo]+=((*SiHC)[iii])->GetEnergy();
		}
	}	
}
=======


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
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217

auto CsI_n_hit = CsIHC->entries();
if (CsI_n_hit>0)
	{
<<<<<<< HEAD
	for ( int iii = 0; iii < CsI_n_hit; iii++)
	{
		CsIdetectorNo= 	((*CsIHC)[iii])->GetDetectorNo();
		CsIpixelNo= 	((*CsIHC)[iii])->GetPixelNo();
		CsIstripNo= 	((*CsIHC)[iii])->GetStripNo();

		if(CsIdetectorNo == 1)
		{
			CsIhe[CsIpixelNo+4*CsIstripNo]+=((*CsIHC)[iii])->GetEnergy();
		}
		else if(CsIdetectorNo == 0)
		{
			CsIdeut[CsIpixelNo+4*CsIstripNo]+=((*CsIHC)[iii])->GetEnergy();
		}
	}
}
if(Si_n_hit>0)
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
=======
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
	
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
}