#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include <G4String.hh>
#include <G4Tubs.hh>
#include <G4SubtractionSolid.hh>
#include <G4Cons.hh>
#include <G4Sphere.hh>
#include <G4UnionSolid.hh>

#include "siliconSD.hh"
#include "siliconHit.hh"
#include "cesiumSD.hh"
#include "cesiumHit.hh"
#include "parameters.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "/home/zalewski/aku/wrk/constants.h"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4UserLimits;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
	public:
		DetectorConstruction();
		virtual ~DetectorConstruction();
		virtual G4VPhysicalVolume* Construct();
		public:
		virtual void ConstructSDandField();
		void SetMaxStep (G4double);

		G4double deut_angle, deut_angle_5;
		G4double helium_angle, helium_angle_5;
		G4double target_angle;
		G4double tarPosition;

		G4double sqlang, sqrang;
		const G4double tar_x = 50.0*mm, tar_y = 40.0*mm;	//Target half-dimensions
		G4double tar_z;

		const G4double sql_dist=(170.0)*mm;
		const G4double sqr_dist=(250.0)*mm;
		
		const G4double world_sizeX = 120*cm;
		const G4double world_sizeY = 8*cm;
		const G4double world_sizeZ	= 120*cm;

		static constexpr bool gasTarget = false;
		G4LogicalVolume *log_pixelSi;
	 	G4LogicalVolume *log_crystalCsI;
		G4ThreeVector vPosition2HTelescope;
		G4ThreeVector vPosition6HeTelescope;
};

#endif