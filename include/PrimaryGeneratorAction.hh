
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4PrimaryVertex.hh"
#include "G4ThreeVector.hh"
#include "/home/zalewski/root/install/include/TLorentzVector.h"
#include "G4ExceptionSeverity.hh"
#include "TTree.h"
#include "DetectorConstruction.hh"

#include <G4ParticleTable.hh>
#include <G4IonTable.hh>
#include <G4EmCalculator.hh>
#include <G4NistManager.hh>
#include <G4Material.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4Orb.hh>
#include <geomdefs.hh>
#include "ParticleInfo.hh"

#include "/home/zalewski/aku/analysis/constants.h"

class G4Event;
class G4Box;


	class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
	{
	public:
		PrimaryGeneratorAction();
		virtual ~PrimaryGeneratorAction();

		double MWPCrange(double position, int clusterMultiplicity);
		G4int geometryID;
		
		// method from the base class
		virtual void GeneratePrimaries(G4Event*);				

		inline bool
		is_target_losses(){return target_losses;}
 
		enum sPar {sMWPC_1_X, sMWPC_1_Y, sMWPC_2_X, sMWPC_2_Y, sTarPos, sLang1, sLang2, sLang3, sRang};
		std::vector<std::string> parNames = {"sMWPC_1_X", "sMWPC_1_Y", "sMWPC_2_X", "sMWPC_2_Y", "sTarPos", "sLang1", "sLang2", "sLang3", "sRang"};


	private:
		G4int entryShift = 0;
		G4int inputEventNo = 0;
		G4int inputTreeLoopCounter = 0;
		bool target_losses;
	 	TTree *inBeamTree;
		G4VSolid *gasCellSolid;
		G4Tubs *deutDiscTube;
		G4Sphere *deutSphere;
		G4Box *sphereCutoff;
		G4VSolid *deutCap;
		G4VSolid *tempDeuterUnion;
		G4double minZ, maxZ;

		double E_tar_loss;
		G4double tar_angle, tar_pos_Z;
		G4double MWPC_equivalent_of_Si;
		double Range;
		
		G4double MWPC_1_X, MWPC_1_Y, MWPC_1_Z;
		G4double MWPC_2_X, MWPC_2_Y, MWPC_2_Z;
		G4double in_nx1, in_nx2, in_ny1, in_ny2;

		G4double in_MWPC_1_X, in_MWPC_1_Y, in_MWPC_1_Z;
		G4double in_MWPC_2_X, in_MWPC_2_Y, in_MWPC_2_Z;

		G4double evX, evY, evZ;
		G4double dX, dY, dZ;
		G4double Tcoef;

		double mass1H;
		double mass2H;
		double mass3H;
		double mass4He;
		double mass5He;
		double mass6He;		
		double massTar;
		double massNeut;

		double excitedStateEnergy_6He;
		double massSum;
		double beam_spot_radius;
		double tar_thick;
		double Ek6He;
		double Ek2H;
		bool flag=true;

		G4double beam_T;
		
		G4ThreeVector MD2H;
		G4ThreeVector MD6He;
		G4Material *Deut_target;
		G4ThreeVector zAxis;
		G4double maxDist;
		G4double minDist;
		
	
		G4Material *silicon_material;
		TLorentzVector *in_lvBeam;
		G4PrimaryVertex *elasticVertex;
		G4PrimaryVertex *inelasticVertex;
		G4PrimaryVertex *dtVertex;
		G4ParticleTable *particletable;
		G4IonTable *iontable;
		G4ParticleDefinition *def1H;
		G4ParticleDefinition *def2H;
		//G4PrimaryParticle *PrimaryParticle_2H_EL;
		G4ParticleDefinition *def3H;
		//G4PrimaryParticle *PrimaryParticle_3H_DT;
		//G4PrimaryParticle *PrimaryParticle_3H_DT;
		G4ParticleDefinition *def4He;
		G4ParticleDefinition *def5He;
		//G4PrimaryParticle *PrimaryParticle_5He_DT;
		//G4PrimaryParticle *PrimaryParticle_5He_DT;
		G4ParticleDefinition *def6He;
		//G4PrimaryParticle *PrimaryParticle_6He_EL;
		G4ParticleDefinition *defTar;
		G4ParticleDefinition *defNeut;
		G4VUserPrimaryParticleInformation *partINFO;
		//ParticleInfo *ElasticINFO;
		//ParticleInfo *INelasticINFO;
		//ParticleInfo *dtParticleInfo;


		
		//G4ThreeVector VertexPosition;

		G4EmCalculator *ELC;
		double get_E(double E, double r, G4Material *mat);
		double get_R(double E, G4Material *mat);
	};

#endif
