
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4PrimaryVertex.hh"
#include "G4ThreeVector.hh"
#include "DetectorConstruction.hh"

//#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4EmCalculator.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"

//class G4ParticleGun;
class G4Event;
class G4Box;

#include <fstream>
#include <iostream>
#include <vector>


	class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
	{
	public:
		PrimaryGeneratorAction();		
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

		double Ex6He;
		double massSum;
		double beam_spot_radius;
		double tar_thick;
		double Ek6He;
		double Ek2H;
		
		G4ThreeVector MD2H;
		G4ThreeVector MD6He;
		G4Material* Deut_target;
	private:
		bool target_losses;
	 

		double E_tar_loss;
		double Range;
		
		double Vertex_X, Vertex_Y, Vertex_Z;		
		
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
	};

#endif