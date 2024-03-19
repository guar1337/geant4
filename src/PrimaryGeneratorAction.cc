#include "PrimaryGeneratorAction.hh"
#include "TFile.h"


#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"



#include <fstream>
#include <iostream>
#include <iomanip> 

	

	PrimaryGeneratorAction::PrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction()
	{

		geometryID = cs::runNo/10;
		printf("\ncs::runNo used is: %d inferred geometry is: %d\n", cs::runNo, geometryID);
		G4String beamFilePath;
		if (geometryID==1)
		{
			tar_thick = (60.0)*um;
			tar_angle = 45.0*TMath::DegToRad();
			beamFilePath = "/home/zalewski/dataTmp/beamSource/beamSource1.root";
			tar_pos_Z = 10.0;
		}

		if (geometryID==2)
		{
			tar_thick = (2.0*60.0)*um;
			tar_angle = 6.0*TMath::DegToRad();
			beamFilePath = "/home/zalewski/dataTmp/beamSource/beamSource2.root";
			tar_pos_Z = 10.0;
		}

		if (geometryID==3)
		{
			tar_thick = (2.0*60.0)*um;
			tar_angle = 0.0*TMath::DegToRad();
			beamFilePath = "/home/zalewski/dataTmp/beamSource/beamSource3.root";
			tar_pos_Z = 10.0 - 2.0;
		}


		
		in_lvBeam = new TLorentzVector();
		TFile *inBeamF = new TFile(beamFilePath.data(),"READ");
		if (!inBeamF->IsOpen()) G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction","beam sourceF not found", FatalException, beamFilePath.data());
		inBeamTree = (TTree*)inBeamF->Get("beamSource");
		inBeamTree->SetMakeClass(0);
		inBeamTree->SetBranchAddress("lvBeam", &in_lvBeam);

		inBeamTree->SetBranchAddress("MWPC_1_X", &in_MWPC_1_X);
		inBeamTree->SetBranchAddress("MWPC_2_X", &in_MWPC_2_X);
		inBeamTree->SetBranchAddress("MWPC_1_Y", &in_MWPC_1_Y);
		inBeamTree->SetBranchAddress("MWPC_2_Y", &in_MWPC_2_Y);

		
		// default particle kinematic
		particletable = G4ParticleTable::GetParticleTable();
		iontable = G4IonTable::GetIonTable();

		defNeut=particletable->FindParticle("neutron");

		ELC= new G4EmCalculator();
		G4NistManager* man = G4NistManager::Instance();
		//man->SetVerbose(0);
		silicon_material = man->FindOrBuildMaterial("G4_Si");

		//ElasticINFO = new ParticleInfo;
		//INelasticINFO = new ParticleInfo;


		//PrimaryParticle_2H_EL = new G4PrimaryParticle();
		//PrimaryParticle_6He_EL = new G4PrimaryParticle();
		//PrimaryParticle_3H_DT = new G4PrimaryParticle();
		//PrimaryParticle_5He_DT = new G4PrimaryParticle();
		


		excitedStateEnergy_6He=1797*keV;
		MWPC_equivalent_of_Si = 660*um;

	}

		
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete ELC;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	inputEventNo = anEvent->GetEventID() - entryShift;
	
	if (geometryID==1 && inputEventNo == 5942259)
	{
		inputTreeLoopCounter++;
		entryShift = inputTreeLoopCounter * 5942259;
		//printf("Loop number %d for 1st geo\nTree entry numer: %d\n",inputTreeLoopCounter, anEvent->GetEventID());
	}
	
	if (geometryID==2 && inputEventNo == 15750492)
	{
		inputTreeLoopCounter++;
		entryShift = inputTreeLoopCounter * 15750492;
		//printf("Loop number %d for 2nd geo\nTree entry numer: %d\n",inputTreeLoopCounter, anEvent->GetEventID());
	}

	if (geometryID==3 && inputEventNo == 35570688)
	{
		inputTreeLoopCounter++;
		entryShift = inputTreeLoopCounter * 35570688;
		//printf("Loop number: %d for 3rd geo\nTree entry numer: %d\n",inputTreeLoopCounter, anEvent->GetEventID());
	}
	
	
	def1H = iontable->GetIon(1,1);
	def2H = iontable->GetIon(1,2);
	def3H = iontable->GetIon(1,3);
	def4He = iontable->GetIon(2,4);
	def5He = iontable->GetIon(2,5);
	def6He = iontable->GetIon(2,6);	
	defTar = iontable->GetIon(1,cs::tarMass);

	massNeut = defNeut->GetPDGMass();
	mass1H = def1H->GetPDGMass();
	mass2H = def2H->GetPDGMass();
	mass3H = def3H->GetPDGMass();
	mass4He=def4He->GetPDGMass();
	mass5He=def5He->GetPDGMass();
	mass6He=def6He->GetPDGMass();
	massTar=defTar->GetPDGMass();
	massNeut=defNeut->GetPDGMass();

	if (flag)
	{
		printf("\n\n\nneut: %f\t1H: %f\t2H: %f\t3H: %f\t4He: %f\t5He: %f\t6He: %f\n\n\n", massNeut, mass1H, mass2H, mass3H, mass4He, mass5He, mass6He);
		flag=false;
	}
	
	
	inBeamTree->GetEntry(inputEventNo);
	//generate vertex position
	MWPC_1_X = in_MWPC_1_X + CLHEP::RandFlat::shoot(0.0,1.25)-0.625 - 1.0;
	MWPC_1_Y = in_MWPC_1_Y + CLHEP::RandFlat::shoot(0.0,1.25)-0.625 - 2.1375;

	MWPC_2_X = in_MWPC_2_X + CLHEP::RandFlat::shoot(0.0,1.25)-0.625 + 0.2;
	MWPC_2_Y = in_MWPC_2_Y + CLHEP::RandFlat::shoot(0.0,1.25)-0.625 - 1.125;
	
	MWPC_1_Z = -816.0;
	MWPC_2_Z = -270.0;

	dX=MWPC_2_X-MWPC_1_X;
	dY=MWPC_2_Y-MWPC_1_Y;
	dZ=MWPC_2_Z-MWPC_1_Z;

	
	double ene_beam = in_lvBeam->E();
	double mom_beam = sqrt(ene_beam*ene_beam - cs::mass6He*cs::mass6He);
	G4ThreeVector vBeam(dX, dY, dZ);
	vBeam.setMag(mom_beam);
	G4LorentzVector lvBeam(vBeam, in_lvBeam->E());

	double tmpTarPos = tar_pos_Z + CLHEP::RandFlat::shoot(-tar_thick/2.0, tar_thick/2.0);
	Tcoef=(cos(tar_angle)*tmpTarPos-sin(tar_angle)*MWPC_1_X - cos(tar_angle)*MWPC_1_Z)/(sin(tar_angle)*dX+cos(tar_angle)*dZ);

	evX = MWPC_1_X + dX*Tcoef;
	evY = MWPC_1_Y + dY*Tcoef;
	evZ = MWPC_1_Z + dZ*Tcoef;

	//printf("evX: %f\tevY: %f\tevZ: %f\n", evX, evY, evZ);
	
	G4ThreeVector VertexPosition(evX,evY,evZ);
	G4ThreeVector minusvectBeam(-dX,-dY,-dZ);
	G4ThreeVector vectBeam(dX, dY, dZ);

	
	//printf("evX: %f\tevY: %f\tevZ: %f\n", evX, evY, evZ);

	//Eloss estimation - CHECK
	//I can only use material which are in my setup, but how to get deut gas?
	//That's why I will use silicon equivalent of deuterium gas
	
	const G4LorentzVector lvTarget(0.0, 0.0, 0.0, massTar);

	const double theta_CM = acos(CLHEP::RandFlat::shoot(-1.0, 1.0));
	const double phi_CM = CLHEP::RandFlat::shoot(0.0, 2.0*CLHEP::pi);
	double realTheta;


//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	START OF elastic PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/*
		G4LorentzVector lv6He_EL(lvBeam);
		G4LorentzVector lv2H_EL(lvTarget);
		G4LorentzVector lvCM_EL = lv6He_EL+lv2H_EL;
		G4ThreeVector boostVect_EL = lvCM_EL.boostVector();

		lv6He_EL.boost(-boostVect_EL);
		lv2H_EL.boost(-boostVect_EL);
		//CM PART

		lv6He_EL.setTheta(theta_CM);
		lv6He_EL.setPhi(phi_CM);
		realTheta = lv6He_EL.vect().angle(lvBeam.vect());

		lv2H_EL.setVect(-lv6He_EL.vect());

		G4LorentzVector lv2H_CM_EL(lv2H_EL);
		G4LorentzVector lv6He_CM_EL(lv6He_EL);

		lv6He_EL.boost(boostVect_EL);
		lv2H_EL.boost(boostVect_EL);


		//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		//xxxxxxxxxxxxxx	ELASTIC PRIMARIES PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		//deuterium
		G4PrimaryParticle *PrimaryParticle_2H_EL = new G4PrimaryParticle(defTar);
		Double_t kineticEnergyH2 = lv2H_EL.e()-lv2H_EL.m();
		if (kineticEnergyH2<0.1 && cs::runNo==30 && cs::tarMass==2)
		{
			//printf("No. of event: %d\tenergy: %f\n", inputEventNo, kineticEnergyH2);
			kineticEnergyH2 = 100.0;
		}
		PrimaryParticle_2H_EL->SetKineticEnergy(kineticEnergyH2);

		PrimaryParticle_2H_EL->SetMomentumDirection(lv2H_EL.vect().unit());
		//I put all particles, that are not in PrimaryVertex into ParticleInfo of 2H
		//That is: lv2H_CM, lv6He_CM and lvBeam
		ParticleInfo *ElasticINFO = new ParticleInfo;
		PrimaryParticle_2H_EL->SetUserInformation(ElasticINFO);
		ElasticINFO->Set_LV_Beam(lvBeam);
		ElasticINFO->Set_LV_2H_CM(lv2H_CM_EL);
		ElasticINFO->Set_LV_6He_CM(lv6He_CM_EL);

		ElasticINFO->Set_MWPC_1_X(in_MWPC_1_X);
		ElasticINFO->Set_MWPC_2_X(in_MWPC_2_X);
		ElasticINFO->Set_MWPC_1_Y(in_MWPC_1_Y);
		ElasticINFO->Set_MWPC_2_Y(in_MWPC_2_Y);

		ElasticINFO->Set_nx1(in_nx1);
		ElasticINFO->Set_nx2(in_nx2);
		ElasticINFO->Set_ny1(in_ny1);
		ElasticINFO->Set_ny2(in_ny2);
		ElasticINFO->Set_thetaCM(realTheta*TMath::RadToDeg());

		//Helium 6
		G4PrimaryParticle *PrimaryParticle_6He_EL = new G4PrimaryParticle(def6He);
		PrimaryParticle_6He_EL->SetKineticEnergy(lv6He_EL.e()-mass6He);		
		PrimaryParticle_6He_EL->SetMomentumDirection(lv6He_EL.vect().unit() );

		elasticVertex = new G4PrimaryVertex(VertexPosition,0);
		elasticVertex->SetPrimary(PrimaryParticle_2H_EL);
		elasticVertex->SetPrimary(PrimaryParticle_6He_EL);
		anEvent->AddPrimaryVertex(elasticVertex);
	
	
*/
//
//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	START OF INELastic PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/*
		G4LorentzVector lv6He_IN=lvBeam;

		G4LorentzVector lv2H_IN=lvTarget;
		G4LorentzVector lvCM_IN = lv6He_IN+lv2H_IN;

		G4ThreeVector boostVect_IN = lvCM_IN.boostVector();
		lv6He_IN.boost(-boostVect_IN);
		lv2H_IN.boost(-boostVect_IN);
		//Let's go inelastic now!
		G4double newEcm	= lv2H_IN.e() + lv6He_IN.e() - excitedStateEnergy_6He;
		G4double newE2H	= (1/(2*newEcm))*(newEcm*newEcm-mass6He*mass6He+lv2H_IN.m()*lv2H_IN.m());
		G4double newE6He	= (1/(2*newEcm))*(newEcm*newEcm+mass6He*mass6He-lv2H_IN.m()*lv2H_IN.m())+excitedStateEnergy_6He;
		G4double newMom2H	= sqrt(newE2H*newE2H-lv2H_IN.m()*lv2H_IN.m());
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
	//
	//
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//xxxxxxxxxxxxxx	INELASTIC PRIMARIES PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
		G4PrimaryParticle *PrimaryParticle_2H_IN = new G4PrimaryParticle(defTar); //0 PP
		PrimaryParticle_2H_IN->SetKineticEnergy(lv2H_IN.e()-lv2H_IN.m());
		PrimaryParticle_2H_IN->SetMomentumDirection(lv2H_IN.vect().unit());
		
		ParticleInfo *INelasticINFO = new ParticleInfo;
		PrimaryParticle_2H_IN->SetUserInformation(INelasticINFO);
		INelasticINFO->Set_LV_Beam(lvBeam);
		INelasticINFO->Set_LV_2H_CM(lv2H_CM_IN);
		INelasticINFO->Set_LV_6He_CM(lv6He_CM_IN);

		G4PrimaryParticle *PrimaryParticle_4He_IN = new G4PrimaryParticle(def4He);
		PrimaryParticle_4He_IN->SetKineticEnergy(lv4He_IN.e()-mass4He);
		PrimaryParticle_4He_IN->SetMomentumDirection(lv4He_IN.vect().unit());
		
		G4PrimaryParticle *PrimaryParticle_Neutron1_IN=new G4PrimaryParticle(defNeut);
		PrimaryParticle_Neutron1_IN->SetKineticEnergy(lvNeut1.e()-massNeut);
		PrimaryParticle_Neutron1_IN->SetMomentumDirection(lvNeut1.vect().unit());

		G4PrimaryParticle * PrimaryParticle_Neutron2_IN=new G4PrimaryParticle(defNeut);
		PrimaryParticle_Neutron2_IN->SetKineticEnergy(lvNeut2.e()-massNeut);
		PrimaryParticle_Neutron2_IN->SetMomentumDirection(lvNeut2.vect().unit());
		
		inelasticVertex = new G4PrimaryVertex(VertexPosition,0);
			
		inelasticVertex->SetPrimary(PrimaryParticle_2H_IN);				//0
		inelasticVertex->SetPrimary(PrimaryParticle_4He_IN);				//1
		//inelasticVertex->SetPrimary(PrimaryParticle_Neutron1_IN);		//2
		//inelasticVertex->SetPrimary(PrimaryParticle_Neutron2_IN);		//3
		anEvent->AddPrimaryVertex(inelasticVertex);
	
	*/

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	d,t kinematics	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

		G4LorentzVector lv6He_DT(lvBeam);
		G4LorentzVector lv2H_DT(lvTarget);
		G4LorentzVector lvCM_DT = lv6He_DT + lv2H_DT;

		G4ThreeVector boostVect_DT = lvCM_DT.boostVector();
		lv6He_DT.boost(-boostVect_DT);
		lv2H_DT.boost(-boostVect_DT);
		//Let's go d,t now!
		G4double Qdt = (mass2H + mass6He) - (mass3H + mass5He);
		G4double eneCM = lv2H_DT.e() + lv6He_DT.e();
		G4double kinECMBefore = lv2H_DT.e() - lv2H_DT.m() + lv6He_DT.e() - lv6He_DT.m();
		G4double kinECMAfter = kinECMBefore + Qdt;
		G4double kinE3H	= (kinECMAfter/(2.0*eneCM))*(kinECMAfter+2.0*mass5He);
		G4double ene3H	= kinE3H + mass3H;
		G4double kinE5He = (kinECMAfter/(2.0*eneCM))*(kinECMAfter+2.0*mass3H);
		G4double ene5He = kinE5He + mass5He;
		G4double newMom3H	= sqrt(ene3H*ene3H-mass3H*mass3H);
		G4double newMom5He = sqrt(ene5He*ene5He-mass5He*mass5He);

		G4ThreeVector v3H(lv2H_DT.vect());
		G4ThreeVector v5He(lv6He_DT.vect());
		v3H.setMag(newMom3H);
		v5He.setMag(newMom5He);
		G4LorentzVector lv3H;
		G4LorentzVector lv5He;
		lv3H.setVectMag(v3H, mass3H);
		lv3H.setTheta(theta_CM);
		lv3H.setPhi(phi_CM);
		realTheta = lv3H.vect().angle(lvBeam.vect());

		lv5He.setVectMag(v5He, mass5He);
		lv5He.setVect(-lv3H);
		G4LorentzVector lv3H_CM(lv3H);
		G4LorentzVector lv5He_CM(lv5He);

		lv3H.boost(boostVect_DT);
		lv5He.boost(boostVect_DT);


		//NEUTRON EMISSION - negligible lifetime of 5He allows for immidiate decay
		G4ThreeVector boostVect_5He = lv5He.boostVector();
		G4double dt_Q = mass5He - (mass4He + massNeut);
		G4double dt_eneCM = mass5He;
		G4double dt_kinECMBefore = 0.0;
		G4double dt_kinECMAfter = dt_kinECMBefore + dt_Q;
		G4double dt_kinENeut = (dt_kinECMAfter/(2.0*dt_eneCM))*(dt_kinECMAfter+2.0*mass4He);
		G4double dt_eneNeut = dt_kinENeut + massNeut;
		G4double dt_kinE4He = (dt_kinECMAfter/(2.0*dt_eneCM))*(dt_kinECMAfter+2.0*massNeut);
		G4double dt_ene4He = dt_kinE4He + mass4He;
		G4double dt_momNeut	= sqrt(dt_eneNeut*dt_eneNeut-massNeut*massNeut);
		G4double dt_mom4He = sqrt(dt_ene4He*dt_ene4He-mass4He*mass4He);

		//printf("dt_Q: %f\tkinENeut: %f\tkinE4He: %f\tmomNeut: %f\tmom4He: %f\n", dt_Q, dt_kinENeut, dt_kinE4He, dt_momNeut, dt_mom4He);

		//Setting neutron vectors
		//4He + n decay of 5He is classical - differences are negligible
		G4double thetaNCM = acos(CLHEP::RandFlat::shoot(-1.0,1.0));
		G4double phiNCM = CLHEP::RandFlat::shoot(0.0, 2.0*CLHEP::pi);

		G4ThreeVector dt_vNeut(0.0,0.0,1.0);
		G4ThreeVector dt_v4He(0.0,0.0,1.0);
		dt_vNeut.setMag(dt_momNeut);
		dt_v4He.setMag(dt_mom4He);
		G4LorentzVector dt_lvNeut;
		G4LorentzVector dt_lv4He;
		dt_lvNeut.setVectMag(dt_vNeut, massNeut);
		dt_lvNeut.setTheta(thetaNCM);

		dt_lv4He.setVectMag(dt_v4He, mass4He);
		dt_lv4He.setTheta(CLHEP::pi-thetaNCM);
		dt_lv4He.setPhi(phiNCM);

		dt_lvNeut.boost(boostVect_DT);
		dt_lv4He.boost(boostVect_5He);

		//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		//xxxxxxxxxxxxxx	DT PRIMARIES PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		//deuterium
		G4PrimaryParticle *PrimaryParticle_3H_DT = new G4PrimaryParticle(def3H);
		PrimaryParticle_3H_DT->SetKineticEnergy(lv3H.e()-lv3H.m());
		PrimaryParticle_3H_DT->SetMomentumDirection(lv3H.vect().unit());
		//I put all particles, that are not in PrimaryVertex into ParticleInfo of 2H
		//That is: lv2H_CM, lv6He_CM and lvBeam
		ParticleInfo *dtParticleInfo = new ParticleInfo;
		PrimaryParticle_3H_DT->SetUserInformation(dtParticleInfo);
		dtParticleInfo->Set_LV_Beam(lvBeam);
		dtParticleInfo->Set_LV_2H_CM(lv3H_CM);
		dtParticleInfo->Set_LV_6He_CM(lv5He_CM);
		dtParticleInfo->Set_LV_5He(lv5He);
		dtParticleInfo->Set_thetaCM(realTheta*TMath::RadToDeg());

		dtParticleInfo->Set_MWPC_1_X(in_MWPC_1_X);
		dtParticleInfo->Set_MWPC_2_X(in_MWPC_2_X);
		dtParticleInfo->Set_MWPC_1_Y(in_MWPC_1_Y);
		dtParticleInfo->Set_MWPC_2_Y(in_MWPC_2_Y);

		dtParticleInfo->Set_nx1(in_nx1);
		dtParticleInfo->Set_nx2(in_nx2);
		dtParticleInfo->Set_ny1(in_ny1);
		dtParticleInfo->Set_ny2(in_ny2);

		//Helium 5

		
		G4PrimaryParticle *PrimaryParticle_5He_DT = new G4PrimaryParticle(def5He);
		PrimaryParticle_5He_DT->SetKineticEnergy(lv5He.e()-mass5He);	
		PrimaryParticle_5He_DT->SetMomentumDirection(lv5He.vect().unit() );
		
		G4PrimaryParticle *PrimaryParticle_4He_DT = new G4PrimaryParticle(def4He);
		PrimaryParticle_4He_DT->SetKineticEnergy(dt_lv4He.e()-mass4He);	
		PrimaryParticle_4He_DT->SetMomentumDirection(dt_lv4He.vect().unit() );
		dtVertex = new G4PrimaryVertex(VertexPosition,0);
		dtVertex->SetPrimary(PrimaryParticle_3H_DT);
		dtVertex->SetPrimary(PrimaryParticle_4He_DT);
		anEvent->AddPrimaryVertex(dtVertex);
	

}
	
G4double PrimaryGeneratorAction::get_E(G4double E, G4double r, G4Material* mat)
{
	G4double R0=ELC->GetRange(E,def6He,mat);
	G4double R;
	R = R0 - r;
	G4double E1=ELC->GetKinEnergy(R,def6He,mat);
	return E1;
}

G4double PrimaryGeneratorAction::get_R(G4double E, G4Material* mat)
{
	G4double R0=ELC->GetRange(E,def6He,mat);
 
	return R0;
}
