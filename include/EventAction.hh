#ifndef EventAction_h
#define EventAction_h 1

<<<<<<< HEAD
#include "G4VUserPrimaryParticleInformation.hh"
#include "G4UserEventAction.hh"
#include "G4LorentzVector.hh"
#include "G4IonTable.hh"

#include "globals.hh"
#include "/home/guar/root/build/include/TTree.h"
#include "/home/guar/root/build/include/TFile.h"
#include "ParticleInfo.hh"
#include "/home/guar/root/build/include/TLorentzVector.h"


=======
#include "G4IonTable.hh"
#include "G4UserEventAction.hh"
#include "G4VUserPrimaryParticleInformation.hh"
#include <G4EmCalculator.hh>
#include "G4NistManager.hh"
#include "globals.hh"
#include "/home/guar/root/build/include/TTree.h"
#include "/home/guar/root/build/include/TFile.h"
#include "G4LorentzVector.hh"
//#include "/home/guar/root/build/include/TLorentzVector.h"
#include "/home/guar/root/build/include/TVector3.h"
#include "ParticleInfo.hh"
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217

/// Event action class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers:
/// - fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap
/// which are collected step by step via the functions
/// - AddAbs(), AddGap()

class EventAction : public G4UserEventAction
{
	public:
	EventAction();
<<<<<<< HEAD
	EventAction(TTree *T);
	virtual ~EventAction();

	virtual void	BeginOfEventAction(const G4Event *event);
	virtual void	EndOfEventAction(const G4Event *event);
	
=======
	EventAction(TTree* T);
	virtual ~EventAction();

	virtual void	BeginOfEventAction(const G4Event* event);
	virtual void	EndOfEventAction(const G4Event* event);
	

	G4Material* Deut_target;
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
	G4double CsIdeut[16];
	G4double SideutX[16];
	G4double SideutY[16];
	G4double CsIhe[16];
	G4double SiheX[16];
	G4double SiheY[16];
<<<<<<< HEAD

	G4double beamT;

	G4double sqlang, sqlde, sqletot, sqlesum, sqltheta, sqlphi;
	G4double sqrang, sqrde, sqretot, sqresum, sqrtheta, sqrphi, t_sqrang;

	G4double Xpos, Ypos, Zpos;
	G4double X6He, Y6He, Z6He;
	G4double X2H, Y2H, Z2H;


	G4double phiCM;
	G4double thetaCM;


	G4ParticleTable *particletable;
	G4IonTable *iontable;
	G4ParticleDefinition *def6He;
	G4ParticleDefinition *def4He;
	G4ParticleDefinition *def2H;

	TTree* tree;
	G4int	fsiliconHCID;
	G4int	fcesiumHCID;


	//From PrimaryVertex
	TVector3 *v2H;
	TVector3 *v6He;

	//From PrimaryVertex
	TLorentzVector *lv2H;
	TLorentzVector *lv6He;
	//From ParticleInfo
	TLorentzVector *lvBeam;
	TLorentzVector *lv2H_CM;
	TLorentzVector *lv6He_CM;
	//temporary holders for G4LorentzVector (aka CLHEP::HepLorentzVector)
	G4LorentzVector *tmp_lvBeam;
	G4LorentzVector *tmp_lv2H_CM;
	G4LorentzVector *tmp_lv6He_CM;

=======
	G4double thetaDeut;
	G4double phiDeut;
	G4double thetaHe;
	G4double phiHe;
	G4double x_helium;
	G4double y_helium;
	G4double x_deut;
	G4double y_deut;

	G4double Xpos;
	G4double Ypos;
	G4double Zpos;

	G4double X6He;
	G4double Y6He;
	G4double Z6He;

	G4double X2H;
	G4double Y2H;
	G4double Z2H;
	G4double mazz;
	G4double stuck;
	G4double Tdeut;
	G4double exp1;
	G4double The;
	G4double phiCM;
	G4double thetaCM;
	G4double phiBEAM;
	G4double theBEAM;
	G4double Tbeam;
	G4double T_CMdeut;
	G4double thetaCMdeut;
	G4double phiCMdeut;
	G4double BEAMenergy;
	G4double deutEDEP;
	G4double heEDEP;
	G4double labAng2H;
	G4double labAngHe;
	G4double exp2;
	G4double mass2H;
	G4double mass6He;
	G4double radThetaCM;
	G4double reTheta2H;
	G4double reTheta6He;

	G4ParticleTable* particletable;
	G4IonTable* iontable;
	G4ParticleDefinition * def6He;
	G4ParticleDefinition * def4He;
	G4ParticleDefinition * def2H;

	private:
	TTree* tree;
	G4EmCalculator* ELC;
	G4int	fsiliconHCID;
	G4int	fcesiumHCID;

	G4ThreeVector *v2H;
	G4ThreeVector *v6He;
	G4ThreeVector *vBEAM;
	G4ThreeVector *vDeutCM;	

	G4LorentzVector *lv2H;
	G4LorentzVector *lv6He;
	G4LorentzVector *lvBEAM;
	G4LorentzVector *lvDeutCM;	
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
};
#endif