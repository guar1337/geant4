//file:///home/zalewski/aku/geant4/include/EventAction.hh
#ifndef EventAction_h
#define EventAction_h 1

#include "G4VUserPrimaryParticleInformation.hh"
#include "G4UserEventAction.hh"
#include "G4LorentzVector.hh"
#include "G4IonTable.hh"

#include "globals.hh"
#include "/home/zalewski/root/build6_20_04/include/TTree.h"
#include "/home/zalewski/root/build6_20_04/include/TFile.h"
#include "ParticleInfo.hh"
#include "/home/zalewski/root/build6_20_04/include/TLorentzVector.h"



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
	G4double CsIdeut[16];
	G4double SideutX[32];
	G4double SideutY[16];
	G4double CsIhe[16];
	G4double SiheX[32];
	G4double SiheY[16];

	G4double beamT;

	G4double sqlang, sqlde, sqletot;
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

	G4float MWPC_1_X, MWPC_1_Y;
	G4float MWPC_2_X, MWPC_2_Y;
	G4float nx1, nx2, ny1, ny2;
	
	float progress = 0.0;
	int consoleWidth;
};
#endif
