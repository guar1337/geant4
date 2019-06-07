//file:///home/guar/aku/geant4/src/EventAction.cc
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
	tree->Branch("CsIdeut", 	CsIdeut, 	"CsIdeut[16]/D");
	tree->Branch("SideutX", 	SideutX, 	"SideutX[16]/D");
	tree->Branch("SideutY",		SideutY, 	"SideutY[16]/D");
	tree->Branch("sqlang",		&sqlang, 	"sqlang/D");
	tree->Branch("fsqlang",		&fsqlang, 	"fsqlang/D");
	tree->Branch("sqlde",		&sqlde,		"sqlde/D");
	tree->Branch("sqletot",		&sqletot,	"sqletot/D");
	tree->Branch("sqlphi",		&sqlphi,		"sqlphi/D");
	tree->Branch("sqltheta",	&sqltheta,	"sqltheta/D");


	//Helium part

	tree->Branch("CsIhe",		CsIhe,		"CsIhe[16]/D");
	tree->Branch("SiheY",		SiheY,		"SiheY[16]/D");
	tree->Branch("SiheX",		SiheX,		"SiheX[16]/D");
	tree->Branch("sqrang",		&sqrang,		"sqrang/D");
	tree->Branch("sqrde",		&sqrde,		"sqrde/D");
	tree->Branch("sqretot",		&sqretot,	"sqretot/D");
	tree->Branch("sqrphi",		&sqrphi,		"sqrphi/D");
	tree->Branch("sqrtheta",	&sqrtheta,	"sqrtheta/D");
	tree->Branch("fsqrang",		&fsqrang, 	"fsqrang/D");

	//BEAM
	tree->Branch("beamT",	&beamT,"beamT/D");
	tree->Branch("thetaCM",	&thetaCM,"thetaCM/D");
	tree->Branch("phiCM",	&phiCM,"phiCM/D");
	tree->Branch("evx",		&evx,	"evx/D");
	tree->Branch("evy",		&evy,	"evy/D");
	tree->Branch("evz",		&evz,	"evz/D");

	tree->Branch("X6He",		&X6He,"X6He/D");
	tree->Branch("Y6He",		&Y6He,"Y6He/D");
	tree->Branch("Z6He",		&Z6He,"Z6He/D");

	tree->Branch("X2H",		&X2H,"X2H/D");
	tree->Branch("Y2H",		&Y2H,"Y2H/D");
	tree->Branch("Z2H",		&Z2H,"Z2H/D");

	tree->Branch("fEve2H",		&fEve2H,"fEve2H/B");
	tree->Branch("fEve6He",		&fEve6He,"fEve6He/B");
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

	fsqlang = 180.0 * (lvBeam->Angle(*v2H))/double(CLHEP::pi);
	fsqrang = 180.0 * (lvBeam->Angle(*v6He))/double(CLHEP::pi);

	beamT =lvBeam->E()-lvBeam->M();

	//double E_IN_CM_deut = IN_CM_deut.e();

	evx = event->GetPrimaryVertex(0)->GetX0()*mm;
	evy = event->GetPrimaryVertex(0)->GetY0()*mm;
	evz = event->GetPrimaryVertex(0)->GetZ0()*mm;
	
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

	//zeroing Silicon detectors
std::fill(SideutX,SideutX+16,0.0);
std::fill(SideutY,SideutY+16,0.0);
std::fill(SiheX,SiheX+16,0.0);
std::fill(SiheY,SiheY+16,0.0);
//zeroing CesiumIodide detectors

std::fill(CsIdeut,CsIdeut+16,0.0);
std::fill(CsIhe,CsIhe+16,0.0);

fEve2H = false;
fEve6He = false;

auto Si_n_hit = SiHC->entries();
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
			
			SiheX[SipixelNo]+=((*SiHC)[iii])->GetEnergy();
			SiheY[SistripNo]+=((*SiHC)[iii])->GetEnergy();
			sqrde+=((*SiHC)[iii])->GetEnergy();
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

			SideutX[SipixelNo]+=((*SiHC)[iii])->GetEnergy();
			SideutY[SistripNo]+=((*SiHC)[iii])->GetEnergy();
			sqlde+=((*SiHC)[iii])->GetEnergy();
		}
	}	
}

auto CsI_n_hit = CsIHC->entries();
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
			CsIdeut[CsIpixelNo+4*CsIstripNo]+=((*CsIHC)[iii])->GetEnergy();
			sqletot+=((*CsIHC)[iii])->GetEnergy();
			//printf("Energia w de:\t%f\tentries: %d\n", CsIhe[CsIpixelNo+4*CsIstripNo], CsI_n_hit);
		}

		else if(CsIdetectorNo == 1)
		{
			CsIhe[CsIpixelNo+4*CsIstripNo]+=((*CsIHC)[iii])->GetEnergy();
			sqretot+=((*CsIHC)[iii])->GetEnergy();
			
			//if (CsIhe[CsIpixelNo+4*CsIstripNo]<0.1) printf("Energia w he:\t%f\n", CsIhe[CsIpixelNo+4*CsIstripNo]);
		}
	}
}


	tree->Fill();


}
