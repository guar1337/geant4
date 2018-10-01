#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "ParticleInfo.hh"


#include <fstream>
#include <iostream>
#include <iomanip> 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......	

	PrimaryGeneratorAction::PrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction()
	{
	// default particle kinematic
	particletable = G4ParticleTable::GetParticleTable();
	iontable = G4IonTable::GetIonTable();
	defProt=particletable->FindParticle("proton");
	defNeut=particletable->FindParticle("neutron");
	defAngel=particletable->FindParticle("geantino");

	ELC= new G4EmCalculator();
	G4NistManager* man = G4NistManager::Instance();
	man->SetVerbose(0);
	Deut_target = man->FindOrBuildMaterial("G4_POLYETHYLENE");
	
	beam_spot_radius=7.5*mm;
	tar_thick=20*um;
	Ex6He=1797*keV;
	}
	
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	PrimaryGeneratorAction::~PrimaryGeneratorAction()
	{
	}
	
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	void
	PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
	{

	def6He = iontable->GetIon(2,6);
	def4He = iontable->GetIon(2,4);
	def2H = iontable->GetIon(1,2);

	mass6He=def6He->GetPDGMass();
	mass4He=def4He->GetPDGMass();
	mass2H=def2H->GetPDGMass();
	massNeut=defNeut->GetPDGMass();
	massSum =mass6He + mass2H;	

	double beam_T = 150*MeV;
	double beam_mean_T = CLHEP::RandGauss::shoot(beam_T,beam_T*0.025/2.355);
	beam_T=beam_mean_T;
	
	//generate vertex position
	G4double relative_Z_position;
	//beam is cut by annular DSSD - sizes are given in mm
	Vertex_X = CLHEP::RandGauss::shoot(0,(beam_spot_radius)/2.355);
	Vertex_Y = CLHEP::RandGauss::shoot(0,(beam_spot_radius)/2.355);
	relative_Z_position = CLHEP::RandFlat::shoot(-1.4142*tar_thick, 1.4142*tar_thick);
	Vertex_Z = relative_Z_position*um-Vertex_X*mm;	
	G4ThreeVector VertexPosition(Vertex_X,Vertex_Y,Vertex_Z);

	//Eloss estimation - CHECK
	E_tar_loss = get_E(beam_T, (tar_thick/2)*um+Vertex_Z*um, Deut_target);	
	beam_T= E_tar_loss;
	const double mom6He = sqrt(beam_T*beam_T+2*mass6He*beam_T);
	G4LorentzVector V_BEAM(0,0,mom6He,beam_T+mass6He);
	G4LorentzVector INI_V_2H(0,0,0,mass2H);



//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	START OF elastic PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//ELASTIC SCATTERING VECTORS
	G4LorentzVector EL_V_6He(0,0,mom6He,beam_T+mass6He);
	G4LorentzVector EL_V_2H(0,0,0,mass2H);
	G4LorentzVector EL_V_CM = EL_V_6He+EL_V_2H;
	G4ThreeVector boostVect = EL_V_CM.boostVector();

	EL_V_6He.boost(-boostVect);
	EL_V_2H.boost(-boostVect);
	//CM PART
	double theta_CM = acos(CLHEP::RandFlat::shoot(-1.0,1.0));
	double phi_CM = CLHEP::RandFlat::shoot(0.0,2*double(CLHEP::pi));	
	EL_V_6He.rotate(0,	theta_CM,	0	);
	EL_V_6He.rotate(0,	0,	phi_CM	);
	EL_V_2H.rotate(0,	-theta_CM,	0	);
	EL_V_2H.rotate(0,	0,	-phi_CM	);
	G4LorentzVector LV_deut_CM_EL(EL_V_2H);
	//G4cout<<" ZONK "<<theta_CM<<G4endl;
	//G4cout<<V_BEAM.e()-mass6He<<" =original"<<""<<G4endl;
	EL_V_6He.boost(boostVect);
	EL_V_2H.boost(boostVect);
	//G4cout<<mass6He<<G4endl;
	//Setting primary particles
	//deuterium
	G4PrimaryParticle *PrimaryParticle_2H=new G4PrimaryParticle(def2H); //2 PP
	PrimaryParticle_2H->SetKineticEnergy(EL_V_2H.e()-mass2H);
	PrimaryParticle_2H->SetMomentumDirection(EL_V_2H.vect().unit());
	
	ParticleInfo *ElasticINFO = new ParticleInfo;
	PrimaryParticle_2H->SetUserInformation(ElasticINFO);
	ElasticINFO->Set_LV_Beam(V_BEAM);
	ElasticINFO->Set_LV_2H_CM(LV_deut_CM_EL);

	//Helium 6
	G4PrimaryParticle * PrimaryParticle_6He=new G4PrimaryParticle(def6He);//3 PP
	PrimaryParticle_6He->SetKineticEnergy(EL_V_6He.e()-mass6He);	
	PrimaryParticle_6He->SetMomentumDirection(EL_V_6He.vect().unit() );

	elasticVertex = new G4PrimaryVertex(VertexPosition,0);
	elasticVertex->SetPrimary(PrimaryParticle_2H);
	elasticVertex->SetPrimary(PrimaryParticle_6He);


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	START OF INelastic PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	

	G4LorentzVector INEL_V_6He=V_BEAM;

	G4LorentzVector INEL_V_2H=INI_V_2H;
	G4LorentzVector INEL_V_CM = INEL_V_6He+INEL_V_2H;

	G4ThreeVector INboostVect = INEL_V_CM.boostVector();
	INEL_V_6He.boost(-INboostVect);
	INEL_V_2H.boost(-INboostVect);
	double pumpkin = INEL_V_2H.e() + INEL_V_6He.e();
	double face	= INEL_V_CM.e();
		//Let's go inelastic now!
	double newEcm	= INEL_V_2H.e() + INEL_V_6He.e() - Ex6He;
	double newE2H	= (1/(2*newEcm))*(newEcm*newEcm-mass6He*mass6He+mass2H*mass2H);
	double newE6He	= (1/(2*newEcm))*(newEcm*newEcm+mass6He*mass6He-mass2H*mass2H)+Ex6He;
	double newMom2H	= sqrt(newE2H*newE2H-mass2H*mass2H);
	double newMom6He = sqrt(newE6He*newE6He-mass6He*mass6He);
	printf("\n\nWhatever\n\n");
	INEL_V_2H.setE(newE2H);
	INEL_V_6He.setE(newE6He+Ex6He);
 
	INEL_V_2H.setRho(newMom2H);
	INEL_V_6He.setRho(newMom6He);

	//CM PART
	INEL_V_6He.rotate(0,	theta_CM,	0	);
	INEL_V_6He.rotate(0,	0,	phi_CM	);
	INEL_V_2H.rotate(0,	-theta_CM,	0	);
	INEL_V_2H.rotate(0,	0,	-phi_CM	);
	
	G4LorentzVector LV_deut_CM_inel(INEL_V_2H);

 INEL_V_6He.boost(INboostVect);
 INEL_V_2H.boost(INboostVect);
//G4cout<<pumpkin-(newE2H+newE6He-mass6He-mass2H)<<"	"<<pumpkin<<G4endl;

		//NEUTRON EMISSION - negligible lifetime of excited 6He allows for immidiate decay
 G4LorentzVector V_of_CM_excited_6He = INEL_V_6He;
 G4ThreeVector DECAYboostVect = V_of_CM_excited_6He.boostVector();
 	//neutrons are sitting and waiting for some action
 double decayE=mass6He+Ex6He-(mass4He+2*massNeut);
 double momRatio = CLHEP::RandFlat::shoot(0.1,1);

		//Setting neutron vectors
		//4He + 2n decay of 6He is classical - differences are negligible
 double neut1Theta_CM = acos(CLHEP::RandFlat::shoot(-1.0,1.0));
 double neut1Phi_CM = CLHEP::RandFlat::shoot(0.0,2*double(CLHEP::pi));
 double neut2Theta_CM = acos(CLHEP::RandFlat::shoot(-1.0,1.0));
 double neut2Phi_CM = CLHEP::RandFlat::shoot(0.0,2*double(CLHEP::pi));
 CLHEP::Hep3Vector V_neut_1 (0,0,1);
 CLHEP::Hep3Vector V_neut_2 (0,0,1);
 V_neut_1.rotate(0,	neut1Theta_CM,	0	);
 V_neut_1.rotate(0,	0,	neut1Phi_CM	);
 V_neut_2.rotate(0,	neut2Theta_CM,	0	);
 V_neut_2.rotate(0,	0,	neut2Phi_CM	);
 CLHEP::Hep3Vector V_4He (0,0,1);
 double V1_cos_V2 = cos(V_neut_1.angle(V_neut_2));
 //G4cout<<"angle w.r.t. another1: "<<V_neut_1.angle(V_neut_2)*180/CLHEP::pi<<G4endl;
 double L_fac = momRatio*momRatio+2*momRatio*V1_cos_V2+1;

 double momNeu1=sqrt(2*decayE*mass4He*massNeut/(mass4He*(1+momRatio*momRatio)+massNeut*L_fac));
 double momNeu2 = momRatio*momNeu1;
 double mom4He	= momNeu1*sqrt(L_fac);
 double eneNeu1 = momNeu1*momNeu1/(2*massNeut);
 double eneNeu2 = momNeu2*momNeu2/(2*massNeut);
 double ene4He	= mom4He*mom4He/(2*mass4He);
 
 G4LorentzVector LV_Neut1(V_neut_1,massNeut+eneNeu1);
 G4LorentzVector LV_Neut2(V_neut_2,massNeut+eneNeu2);
 G4LorentzVector LV_4He=(-LV_Neut1-LV_Neut2);
 LV_4He.setE(ene4He+mass4He);
 //LV_4He.setRho(mom4He);
 
 
//mass6He-(mass4He+2*massNeut)

 LV_4He.boost(DECAYboostVect);
 LV_Neut1.boost(DECAYboostVect);
 LV_Neut2.boost(DECAYboostVect);
 //auto T4He = LV_4He.e()-mass4He;
 //auto TN1 = LV_Neut1.e()-massNeut;
 //auto TN2 = LV_Neut2.e()-massNeut;
 //G4cout<<"Yo mamma: "<<(beam_T-(LV_Neut1.e()-massNeut+LV_Neut2.e()-massNeut+LV_4He.e()-mass4He+INEL_V_2H.e()-mass2H))/MeV<<" neu1: "<<TN1-eneNeu1<<" neu2: "<<TN2-eneNeu2<<" differece: "<<(LV_4He.rho()+LV_Neut2.rho()+LV_Neut1.rho())*eV<<G4endl;

 //G4cout<<(beam_T-(T4He+TN1+TN2+INEL_V_2H.e()-mass2H)+(mass6He-(mass4He+2*massNeut)))/keV<<" =before coll, after coll= "<<G4endl;


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	PRIMARIES PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 G4PrimaryParticle *	INEL_PrimaryParticle_2H=new G4PrimaryParticle(def2H); //0 PP
 INEL_PrimaryParticle_2H->SetKineticEnergy(INEL_V_2H.e()-mass2H);
 INEL_PrimaryParticle_2H->SetMomentumDirection(INEL_V_2H.vect().unit());
 
 ParticleInfo *INelasticINFO = new ParticleInfo;
 INEL_PrimaryParticle_2H->SetUserInformation(INelasticINFO);
 INelasticINFO->Set_LV_Beam(V_BEAM);
 INelasticINFO->Set_LV_2H_CM(LV_deut_CM_inel);

 G4PrimaryParticle * INEL_pp4He=new G4PrimaryParticle(def4He);
 INEL_pp4He->SetKineticEnergy(LV_4He.e()-mass4He);
 INEL_pp4He->SetMomentumDirection(LV_4He.vect().unit());
 
 G4PrimaryParticle * INEL_ppNeutron_1=new G4PrimaryParticle(defNeut);
 INEL_ppNeutron_1->SetKineticEnergy(LV_Neut1.e()-massNeut);
 INEL_ppNeutron_1->SetMomentumDirection(LV_Neut1.vect().unit());

 G4PrimaryParticle * INEL_ppNeutron_2=new G4PrimaryParticle(defNeut);
 INEL_ppNeutron_2->SetKineticEnergy(LV_Neut2.e()-massNeut);
 INEL_ppNeutron_2->SetMomentumDirection(LV_Neut2.vect().unit());
 
 inelasticVertex = new G4PrimaryVertex(VertexPosition,0);
	
 inelasticVertex->SetPrimary(INEL_PrimaryParticle_2H);				//0
 inelasticVertex->SetPrimary(INEL_pp4He);				//1
 inelasticVertex->SetPrimary(INEL_ppNeutron_1);		//2
 inelasticVertex->SetPrimary(INEL_ppNeutron_2);		//3


int EL_or_INEL = rand() % 10;
if (EL_or_INEL == 0)
{
	anEvent->AddPrimaryVertex(inelasticVertex);
}
else
{
	anEvent->AddPrimaryVertex(elasticVertex);
}


anEvent->AddPrimaryVertex(elasticVertex);
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	END OF INELASTIC PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	}
	
	double 
	PrimaryGeneratorAction::get_E(double E, double r, G4Material* mat)
{
	double R0=ELC->GetRange(E,def6He,mat);
	double R;
	R = R0 - r;
	double E1=ELC->GetKinEnergy(R,def6He,mat);
	return E1;
}

double
PrimaryGeneratorAction::get_R(double E, G4Material* mat)
{
	double R0=ELC->GetRange(E,def6He,mat);
 
	return R0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

