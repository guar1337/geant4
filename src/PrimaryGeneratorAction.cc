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
	excitedStateEnergy_6He=1797*keV;
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

	G4double beam_T = 150*MeV;
	G4double beam_mean_T = CLHEP::RandGauss::shoot(beam_T,beam_T*0.025/2.355);
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
	G4LorentzVector lvBeam(0,0,mom6He,beam_T+mass6He);
	G4LorentzVector lvTarget(0,0,0,mass2H);



//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	START OF elastic PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//ELASTIC SCATTERING VECTORS
	G4LorentzVector lv6He_EL(0,0,mom6He,beam_T+mass6He);
	G4LorentzVector lv2H_EL(0,0,0,mass2H);
	G4LorentzVector lvCM_EL = lv6He_EL+lv2H_EL;
	G4ThreeVector boostVect_EL = lvCM_EL.boostVector();

	lv6He_EL.boost(-boostVect_EL);
	lv2H_EL.boost(-boostVect_EL);
	//CM PART
	double theta_CM = acos(CLHEP::RandFlat::shoot(-1.0,1.0));
	double phi_CM = CLHEP::RandFlat::shoot(0.0,2*double(CLHEP::pi));	
	lv6He_EL.rotate(0,	theta_CM,	0	);
	lv6He_EL.rotate(0,	0,	phi_CM	);
	lv2H_EL.rotate(0,	-theta_CM,	0	);
	lv2H_EL.rotate(0,	0,	-phi_CM	);
	G4LorentzVector lv2H_CM_EL(lv2H_EL);
	G4LorentzVector lv6He_CM_EL(lv6He_EL);
	lv6He_EL.boost(boostVect_EL);
	lv2H_EL.boost(boostVect_EL);
	//Setting primary particles
	//deuterium
	G4PrimaryParticle *PrimaryParticle_2H_EL = new G4PrimaryParticle(def2H); //2 PP
	PrimaryParticle_2H_EL->SetKineticEnergy(lv2H_EL.e()-mass2H);
	PrimaryParticle_2H_EL->SetMomentumDirection(lv2H_EL.vect().unit());
	//I put all particles, that are not in PrimaryVertex into ParticleInfo of 2H
	//That is: lv2H_CM, lv6He_CM and lvBeam
	ParticleInfo *ElasticINFO = new ParticleInfo;
	PrimaryParticle_2H_EL->SetUserInformation(ElasticINFO);
	ElasticINFO->Set_LV_Beam(lvBeam);
	ElasticINFO->Set_LV_2H_CM(lv2H_CM_EL);
	ElasticINFO->Set_LV_6He_CM(lv6He_CM_EL);

	//Helium 6
	G4PrimaryParticle *PrimaryParticle_6He_EL = new G4PrimaryParticle(def6He);//3 PP
	PrimaryParticle_6He_EL->SetKineticEnergy(lv6He_EL.e()-mass6He);	
	PrimaryParticle_6He_EL->SetMomentumDirection(lv6He_EL.vect().unit() );

	elasticVertex = new G4PrimaryVertex(VertexPosition,0);
	elasticVertex->SetPrimary(PrimaryParticle_2H_EL);
	elasticVertex->SetPrimary(PrimaryParticle_6He_EL);


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	START OF INELastic PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
G4LorentzVector lv6He_IN=lvBeam;
G4LorentzVector lv2H_IN=lvTarget;
G4LorentzVector lvCM_IN = lv6He_IN+lv2H_IN;

G4ThreeVector boostVect_IN = lvCM_IN.boostVector();
lv6He_IN.boost(-boostVect_IN);
lv2H_IN.boost(-boostVect_IN);
//Let's go inelastic now!
G4double newEcm	= lv2H_IN.e() + lv6He_IN.e() - excitedStateEnergy_6He;
G4double newE2H	= (1/(2*newEcm))*(newEcm*newEcm-mass6He*mass6He+mass2H*mass2H);
G4double newE6He	= (1/(2*newEcm))*(newEcm*newEcm+mass6He*mass6He-mass2H*mass2H)+excitedStateEnergy_6He;
G4double newMom2H	= sqrt(newE2H*newE2H-mass2H*mass2H);
G4double newMom6He = sqrt(newE6He*newE6He-mass6He*mass6He);

lv2H_IN.setE(newE2H);
lv6He_IN.setE(newE6He+excitedStateEnergy_6He);
 
lv2H_IN.setRho(newMom2H);
lv6He_IN.setRho(newMom6He);

//CM PART
lv6He_IN.rotate(0,	theta_CM,	0);
lv6He_IN.rotate(0,	0,				phi_CM	);
lv2H_IN.rotate(0,		-theta_CM,	0	);
lv2H_IN.rotate(0,		0,				-phi_CM	);
	
G4LorentzVector lv2H_CM_IN(lv2H_IN);
G4LorentzVector lv6He_CM_IN(lv6He_IN);

lv6He_IN.boost(boostVect_IN);
lv2H_IN.boost(boostVect_IN);


//NEUTRON EMISSION - negligible lifetime of excited 6He allows for immidiate decay
G4LorentzVector lv6He_exc_IN = lv6He_IN;
G4ThreeVector boostVect_6He = lv6He_exc_IN.boostVector();
//neutrons are sitting and waiting for some action
G4double decayE = mass6He + excitedStateEnergy_6He -(mass4He + 2*massNeut);
G4double momRatio = CLHEP::RandFlat::shoot(0.1,1);

//Setting neutron vectors
//4He + 2n decay of 6He is classical - differences are negligible
G4double neut1Theta_CM = acos(CLHEP::RandFlat::shoot(-1.0,1.0));
G4double neut1Phi_CM = CLHEP::RandFlat::shoot(0.0,2*G4double(CLHEP::pi));
G4double neut2Theta_CM = acos(CLHEP::RandFlat::shoot(-1.0,1.0));
G4double neut2Phi_CM = CLHEP::RandFlat::shoot(0.0,2*G4double(CLHEP::pi));

G4ThreeVector vneu1 (0,0,1);
G4ThreeVector vneu2 (0,0,1);
vneu1.rotate(0,	neut1Theta_CM,	0	);
vneu1.rotate(0,	0,	neut1Phi_CM	);
vneu2.rotate(0,	neut2Theta_CM,	0	);
vneu2.rotate(0,	0,	neut2Phi_CM	);

G4ThreeVector v4He (0,0,1);
G4double cos_angle_between_neutrons = cos(vneu1.angle(vneu2));
G4double L_fac = momRatio * momRatio + 2*momRatio * cos_angle_between_neutrons + 1;
G4double momNeu1=sqrt(2*decayE*mass4He*massNeut/(mass4He*(1+momRatio*momRatio)+massNeut*L_fac));
G4double momNeu2 = momRatio*momNeu1;
G4double mom4He	= momNeu1*sqrt(L_fac);
G4double eneNeu1 = momNeu1*momNeu1/(2*massNeut);
G4double eneNeu2 = momNeu2*momNeu2/(2*massNeut);
G4double ene4He	= mom4He*mom4He/(2*mass4He);
 
G4LorentzVector lvNeut1(vneu1,massNeut+eneNeu1);
G4LorentzVector lvNeut2(vneu2,massNeut+eneNeu2);
G4LorentzVector lv4He_IN=(-lvNeut1-lvNeut2);
lv4He_IN.setE(ene4He+mass4He);

lv4He_IN.boost(boostVect_6He);
lvNeut1.boost(boostVect_6He);
lvNeut2.boost(boostVect_6He);

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	PRIMARIES PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
G4PrimaryParticle *	PrimaryParticle_2H_IN=new G4PrimaryParticle(def2H); //0 PP
PrimaryParticle_2H_IN->SetKineticEnergy(lv2H_IN.e()-mass2H);
PrimaryParticle_2H_IN->SetMomentumDirection(lv2H_IN.vect().unit());
 
ParticleInfo *INelasticINFO = new ParticleInfo;
PrimaryParticle_2H_IN->SetUserInformation(INelasticINFO);
INelasticINFO->Set_LV_Beam(lvBeam);
INelasticINFO->Set_LV_2H_CM(lv2H_CM_IN);
INelasticINFO->Set_LV_6He_CM(lv6He_CM_IN);

G4PrimaryParticle * PrimaryParticle_4He_IN=new G4PrimaryParticle(def4He);
PrimaryParticle_4He_IN->SetKineticEnergy(lv4He_IN.e()-mass4He);
PrimaryParticle_4He_IN->SetMomentumDirection(lv4He_IN.vect().unit());
/*
G4PrimaryParticle * PrimaryParticle_Neutron1_IN=new G4PrimaryParticle(defNeut);
PrimaryParticle_Neutron1_IN->SetKineticEnergy(lvNeut1.e()-massNeut);
PrimaryParticle_Neutron1_IN->SetMomentumDirection(lvNeut1.vect().unit());

G4PrimaryParticle * PrimaryParticle_Neutron2_IN=new G4PrimaryParticle(defNeut);
PrimaryParticle_Neutron2_IN->SetKineticEnergy(lvNeut2.e()-massNeut);
PrimaryParticle_Neutron2_IN->SetMomentumDirection(lvNeut2.vect().unit());
*/
inelasticVertex = new G4PrimaryVertex(VertexPosition,0);
	
inelasticVertex->SetPrimary(PrimaryParticle_2H_IN);				//0
inelasticVertex->SetPrimary(PrimaryParticle_4He_IN);				//1
//inelasticVertex->SetPrimary(PrimaryParticle_Neutron1_IN);		//2
//inelasticVertex->SetPrimary(PrimaryParticle_Neutron2_IN);		//3

inelasticVertex->SetWeight(0.1);
elasticVertex->SetWeight(1.0);
anEvent->AddPrimaryVertex(elasticVertex);
anEvent->AddPrimaryVertex(inelasticVertex);




//anEvent->AddPrimaryVertex(elasticVertex);
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	END OF INELASTIC PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

}
	
	G4double 
	PrimaryGeneratorAction::get_E(G4double E, G4double r, G4Material* mat)
{
	G4double R0=ELC->GetRange(E,def6He,mat);
	G4double R;
	R = R0 - r;
	G4double E1=ELC->GetKinEnergy(R,def6He,mat);
	return E1;
}

G4double
PrimaryGeneratorAction::get_R(G4double E, G4Material* mat)
{
	G4double R0=ELC->GetRange(E,def6He,mat);
 
	return R0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

