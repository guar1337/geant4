
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4PrimaryVertex.hh"
#include "G4ThreeVector.hh"
#include "/home/zalewski/root/build6_20_04/include/TLorentzVector.h"
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

#include "/home/zalewski/aku/wrk/constants.h"

class G4Event;
class G4Box;


	class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
	{
	public:
		PrimaryGeneratorAction();
		virtual ~PrimaryGeneratorAction();

		double MWPCrange(double position, int clusterMultiplicity);
		
		// method from the base class
		virtual void GeneratePrimaries(G4Event*);				

		inline bool
		is_target_losses(){return target_losses;}
 


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
		
		G4float MWPC_1_X, MWPC_1_Y, MWPC_1_Z;
		G4float MWPC_2_X, MWPC_2_Y, MWPC_2_Z;
		G4float in_nx1, in_nx2, in_ny1, in_ny2;

		G4float in_MWPC_1_X, in_MWPC_1_Y, in_MWPC_1_Z;
		G4float in_MWPC_2_X, in_MWPC_2_Y, in_MWPC_2_Z;

		G4double evX, evY, evZ;
		G4double dX, dY, dZ;
		G4double Tcoef;

		double mass6He;
		double mass4He;
		double massTar;
		double massNeut;

		double excitedStateEnergy_6He;
		double massSum;
		double beam_spot_radius;
		double tar_thick;
		double Ek6He;
		double Ek2H;

		G4double beam_T;
		
		G4ThreeVector MD2H;
		G4ThreeVector MD6He;
		G4Material *Deut_target;
		G4ThreeVector zAxis;
		G4double maxDist;
		G4double minDist;
		
	
		G4Material *silicon_material;
		TLorentzVector *in_lvBeam;
		G4ParticleTable *particletable;
		G4IonTable *iontable;
		G4ParticleDefinition *def6He;
		G4ParticleDefinition *def4He;
		G4ParticleDefinition *defTar;
		G4ParticleDefinition *defNeut;
		G4VUserPrimaryParticleInformation *partINFO;

		
		//G4ThreeVector VertexPosition;
		G4PrimaryVertex *elasticVertex;
		G4PrimaryVertex *inelasticVertex;
		G4EmCalculator *ELC;
		double get_E(double E, double r, G4Material *mat);
		double get_R(double E, G4Material *mat);
	};

#endif
