//file:///home/zalewski/aku/geant4/include/EventAction.hh
#ifndef EventAction_h
#define EventAction_h 1

#include "G4VUserPrimaryParticleInformation.hh"
#include "G4UserEventAction.hh"
#include "G4LorentzVector.hh"
#include "G4IonTable.hh"

#include "globals.hh"
#include "/home/zalewski/root/install/include/TTree.h"
#include "/home/zalewski/root/install/include/TFile.h"
#include "ParticleInfo.hh"
#include "/home/zalewski/root/install/include/TLorentzVector.h"
#include "/home/zalewski/aku/analysis/constants.h"



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
	EventAction(TTree *T);
	virtual ~EventAction();

	virtual void	BeginOfEventAction(const G4Event *event);
	virtual void	EndOfEventAction(const G4Event *event);
	
	G4bool fEve2H;
	G4bool fEve6He;
	G4double cal_CsI_L[16];
	G4double cal_SQX_L[32];
	G4double cal_SQY_L[16];
	G4double cal_CsI_R[16];
	G4double cal_SQX_R[32];
	G4double cal_SQY_R[16];
	G4double cal_SQ300[16];
	G4double  tSQX_L[16];
	G4double  tSQX_R[16];

	G4double beamT;
	G4double tof;
	G4double corrthetaCM;
	G4double thetaCM;

	G4double sqlang, sqlde, sqletot, sqrang5He;
	G4double sqrang, sqrde, sqretot;

	G4int SipixelNo;
	G4int SistripNo;
	G4int SidetectorNo;
	G4int CsIpixelNo;
	G4int CsIstripNo;
	G4int CsIdetectorNo;

	G4double evx, evy, evz;
	G4double X6He, Y6He, Z6He;
	G4double X2H, Y2H, Z2H;


	G4double phiCM;


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
	TLorentzVector *lv5He;
	TLorentzVector *lv6He;
	//From ParticleInfo
	TLorentzVector *lvBeam;
	TLorentzVector *lv2H_CM;
	TLorentzVector *lv6He_CM;
	//temporary holders for G4LorentzVector (aka CLHEP::HepLorentzVector)
	G4LorentzVector *tmp_lvBeam;
	G4LorentzVector *tmp_lv2H_CM;
	G4LorentzVector *tmp_lv5He;
	G4LorentzVector *tmp_lv6He_CM;

	G4double MWPC_1_X, MWPC_1_Y;
	G4double MWPC_2_X, MWPC_2_Y;
	G4double nx1, nx2, ny1, ny2;

	G4int geo;
	G4int fSQX_R_strip, fSQX_L_strip;
	
	float progress = 0.0;
	int consoleWidth;
};
#endif
