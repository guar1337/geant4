#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4String.hh"

#include "siliconSD.hh"
#include "siliconHit.hh"
#include "cesiumSD.hh"
#include "cesiumHit.hh"
#include "parameters.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"


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

<<<<<<< HEAD
		const G4float deut_angle	= 65.0*deg;
		const G4float helium_angle	= 15.0*deg;
		const G4float target_angle	= 45.0*deg;
		const G4float tar_x = 50.0*mm, tar_y = 40.0*mm, tar_z = 100.0*um;	//Target half-dimensions

		const G4float sql_dist=170.0*mm;
		const G4float sqr_dist=250.0*mm;



		const G4double world_sizeX = 120*cm;
		const G4double world_sizeY = 7*cm;
		const G4double world_sizeZ	= 120*cm;

		G4LogicalVolume *log_pixelSi;
	 	G4LogicalVolume *log_crystalCsI;
=======
	 G4LogicalVolume* Si_pixel_log;
	 G4LogicalVolume* CsI_crystal_log;
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
};

#endif