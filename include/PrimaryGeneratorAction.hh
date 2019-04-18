
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4PrimaryVertex.hh"
#include "G4ThreeVector.hh"
<<<<<<< HEAD
#include "/home/guar/root/build/include/TLorentzVector.h"
#include "G4ExceptionSeverity.hh"
#include "TTree.h"
#include "DetectorConstruction.hh"

=======
#include "DetectorConstruction.hh"

//#include "G4ParticleGun.hh"
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4EmCalculator.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"

<<<<<<< HEAD
class G4Event;
class G4Box;

=======
//class G4ParticleGun;
class G4Event;
class G4Box;

#include <fstream>
#include <iostream>
#include <vector>

>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217

	class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
	{
	public:
<<<<<<< HEAD
		PrimaryGeneratorAction();
=======
		PrimaryGeneratorAction();		
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
		virtual ~PrimaryGeneratorAction();
		
		// method from the base class
		virtual void GeneratePrimaries(G4Event*);				

		inline bool
		is_target_losses(){return target_losses;}
 

		double mass6He;
		double mass4He;
		double mass1H;
		double mass2H;
		double massNeut;
		double mass10Li;

<<<<<<< HEAD
		double excitedStateEnergy_6He;
=======
		double Ex6He;
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
		double massSum;
		double beam_spot_radius;
		double tar_thick;
		double Ek6He;
		double Ek2H;
<<<<<<< HEAD

		G4double beam_T;
		
		G4ThreeVector MD2H;
		G4ThreeVector MD6He;
		G4Material *Deut_target;
	private:
		bool target_losses;
	 	TTree *inBeamTree;
=======
		
		G4ThreeVector MD2H;
		G4ThreeVector MD6He;
		G4Material* Deut_target;
	private:
		bool target_losses;
	 
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217

		double E_tar_loss;
		double Range;
		
		double Vertex_X, Vertex_Y, Vertex_Z;		
<<<<<<< HEAD
	
		TLorentzVector *in_lvBeam;
		G4ParticleTable *particletable;
		G4IonTable *iontable;
		G4ParticleDefinition *def6He;
		G4ParticleDefinition *def4He;
		G4ParticleDefinition *def2H;
		G4ParticleDefinition *defProt;
		G4ParticleDefinition *defNeut;
		G4ParticleDefinition *defAngel;
		G4VUserPrimaryParticleInformation *partINFO;
		
		//G4ThreeVector VertexPosition;
		G4PrimaryVertex *elasticVertex;
		G4PrimaryVertex *inelasticVertex;
		G4EmCalculator *ELC;
		double 
		get_E(double E, double r, G4Material *mat);
		double
		get_R(double E, G4Material *mat);
=======
		
		//G4ParticleGun*	fParticleGun; // pointer a to G4 gun class
		
		G4ParticleTable* particletable;
		G4IonTable* iontable;
		G4ParticleDefinition * def6He;
		G4ParticleDefinition * def4He;
		G4ParticleDefinition * def2H;
		G4ParticleDefinition * defProt;
		G4ParticleDefinition * defNeut;
		G4ParticleDefinition * defAngel;
		G4VUserPrimaryParticleInformation* partINFO;
		
		//G4ThreeVector VertexPosition;
		G4PrimaryVertex* elasticVertex;
		G4PrimaryVertex* inelasticVertex;
		G4EmCalculator* ELC;
		double 
		get_E(double E, double r, G4Material* mat);
		double
		get_R(double E, G4Material* mat);
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
	};

#endif