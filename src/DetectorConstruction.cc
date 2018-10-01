#include "DetectorConstruction.hh"

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{	
	G4bool checkOverlaps = true;
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();
	G4Material* Si_mat = nist->FindOrBuildMaterial("G4_Si");
	G4Material* CsI_mat = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
	
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//xxxxxxxx	DETECTORS PARAMETERS	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	G4float deut_angle=45;
	G4float helium_angle=15;

	G4float sql_dist=170;
	G4float sqr_dist=250;

	
	// Option to switch on/off checking of volumes overlaps
	//
	
			
	//	
	// World
	//
	G4double world_sizeX = 120*cm;
	G4double world_sizeY = 7*cm;
	G4double world_sizeZ	= 120*cm;
	G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
	
	
	G4Box* solidWorld =	
	new G4Box("World",	//its name
	0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);	//its size
	
	G4LogicalVolume* log_wrld =	
	new G4LogicalVolume(solidWorld,	//its solid
	world_mat,	//its material
	"log_World");	//its name
	
	G4VPhysicalVolume* phys_wrld = 
	new G4PVPlacement(0,	//no rotation
	G4ThreeVector(),	//at (0,0,0)
	log_wrld,	//its logical volume
	"phys_World",	//its name
	0,	//its mother	volume
	false,	//no boolean operation
	0,	//copy number
	checkOverlaps);	//overlaps checking
	log_wrld->SetVisAttributes(new G4VisAttributes(false));
	//	
	// Target
	// 

	G4Material* CD2 = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
	G4ThreeVector tar_pos_vect = G4ThreeVector(0, 0, 0);
	G4RotationMatrix* tar_rot_mtrx = new G4RotationMatrix();
	tar_rot_mtrx->rotateX(0.*deg);
	tar_rot_mtrx->rotateY(-45.*deg);
	tar_rot_mtrx->rotateZ(0.*rad); 
 
	// CD2 target	
	G4double tar_x =	2.*cm, tar_y = 2.5*cm, tar_z = 10.*um;
	G4Box* tar_solid =	
	new G4Box("foil", 
	tar_x,
	tar_y,
	tar_z);
	
	G4LogicalVolume* log_tar =	
	new G4LogicalVolume(tar_solid,	//its solid
	CD2,	//its material - I can take poly
	"log_tar");	//its name
	
	new G4PVPlacement(tar_rot_mtrx,	//no rotation
	tar_pos_vect,	//at position
	log_tar,	//its logical volume
	"phys_Target",	//its name
	log_wrld,	//its mother	volume
	false,	//no boolean operation
	0,	//copy number
	checkOverlaps);	//overlaps checking
	
	//log_tar->SetVisAttributes(new G4VisAttributes(false));






//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxx	BOXES - SOLIDS	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	// Det box
	G4double telescope_box_x =	35*mm, telescope_box_y = 35*mm, telescope_box_z = 35*mm;
	G4Box* telescope_box =	
	new G4Box("telescope_box", 
	telescope_box_x,
	telescope_box_y,
	telescope_box_z);

	// silicon strip box 
	G4double Si_strip_x =	32*mm, Si_strip_y = 2*mm, Si_strip_z = 0.5*mm;
	G4Box* Si_strip_box = 
	new G4Box("Si strip",		//name
				Si_strip_x,			//X
				Si_strip_y,			//Y
				Si_strip_z);		//Z


	// CsI strip box 
	G4double CsI_strip_x =	34*mm, CsI_strip_y = 8.25*mm, CsI_strip_z = 34*mm;
	G4Box* CsI_strip_box = 
	 new G4Box("CsI strip",			//name
			CsI_strip_x,		//X
			CsI_strip_y,		//Y
			CsI_strip_z);		//Z


	//pixel of Si box
	G4double Si_pixel_x =	2*mm, Si_pixel_y = 2*mm, Si_pixel_z = 0.5*mm;
	G4Box* Si_pixel_box = 
		new G4Box("Si pixel",		//name
		Si_pixel_x,			//X
		Si_pixel_y,			//Y
		Si_pixel_z);	 		//Z


	// Crystal box		//x and y are actually 8.25
	G4double CsI_crystal_x =	8.25*mm, CsI_crystal_y = 8.25*mm, CsI_crystal_z = 34*mm;
	G4Box* CsI_crystal_box =	
	new G4Box("CsI_crystal", 
	CsI_crystal_x,
	CsI_crystal_y,
	CsI_crystal_z);


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxx	LOGICAL VOLUMES	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	// Telescope LogVolume
	G4LogicalVolume* telescope_log=	
	new G4LogicalVolume(telescope_box,	//its solid
			world_mat,			//its material - I can take poly
			"telescope_log");		//its name

	// CsI strip LogVolume
	G4LogicalVolume* CsI_strip_log =	
	new G4LogicalVolume(CsI_strip_box,		//its solid
	world_mat,			//its material - I can take poly
	"CsI_strip_log");		//its name

	// Silicon strip LogVolume
	G4LogicalVolume* Si_strip_log =	
	new G4LogicalVolume(Si_strip_box,		//its solid
	world_mat,			//its material - I can take poly
	"Si_strip_log");		//its name

	// Crystal LogVolume
			CsI_crystal_log =	
	new G4LogicalVolume(CsI_crystal_box,		//its solid
	CsI_mat,			//its material - I can take poly
	"CsI_crystal_log");		//its name

	// Pixel LogVolume
			Si_pixel_log =	
	new G4LogicalVolume(Si_pixel_box,		//its solid
			Si_mat,				//its material
			"Si_pixel_log");		//its name


	// Visibility

telescope_log-> SetVisAttributes (G4VisAttributes::GetInvisible());
Si_strip_log-> SetVisAttributes (G4VisAttributes::GetInvisible());
CsI_strip_log-> SetVisAttributes (G4VisAttributes::GetInvisible());
//CsI_crystal_log-> SetVisAttributes (G4VisAttributes::GetInvisible());
//Si_pixel_log-> SetVisAttributes (G4VisAttributes::GetInvisible());


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxx	DETECTORS POSITIONING	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	//Deuterium telescope container positioning
G4ThreeVector d_tele_pos_vect = 
	G4ThreeVector(sin(deut_angle*deg)*(35+sql_dist)*mm,
	0,
	cos(deut_angle*deg)*(35+sql_dist)*mm);

	//Helium telescope container positioning
G4ThreeVector he_tele_pos_vect = 
	G4ThreeVector(-sin(helium_angle*deg)*(35+sqr_dist)*mm,
	0,
	cos(helium_angle*deg)*(35+sqr_dist)*mm);


	//Deuterium telescope container rotation
G4RotationMatrix* deut_tele_rot_mtrx = new G4RotationMatrix();
deut_tele_rot_mtrx->rotateX(0.*deg);
deut_tele_rot_mtrx->rotateY(-deut_angle*deg);
deut_tele_rot_mtrx->rotateZ(0.*rad); 


	//Helium telescope container rotation
G4RotationMatrix* he_tele_rot_mtrx = new G4RotationMatrix();
he_tele_rot_mtrx->rotateX(0.*deg);
he_tele_rot_mtrx->rotateY(helium_angle*deg);
he_tele_rot_mtrx->rotateZ(0.*rad);


	

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxx	DETECTORS PLACEMENT	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	// Telescopes placement

	// Deuterium telescope placement	
	new G4PVPlacement(deut_tele_rot_mtrx,	//rotation mtrx
	d_tele_pos_vect,		//at (0,0,0)
	telescope_log,		//its logical volume
	"helium telescope",		//its name
	log_wrld,		//its mother	volume
	false,				//no boolean operation
	0,					//copy number - for deuterium
	checkOverlaps);			//overlaps checking


	// Helium telescope placement
	new G4PVPlacement(he_tele_rot_mtrx,		//rotation mtrx
	he_tele_pos_vect,	 	//at (0,0,0)
	telescope_log,		//its logical volume
	"deuterium telescope",		//its name
	log_wrld,		//its mother	volume
	false,		//no boolean operation
	1,		//copy number - for helium
	checkOverlaps);		//overlaps checking

	// Strips placement

	// Silicon strips placement
for (int iii=0; iii<16; iii++)
	{
new G4PVPlacement(0,						//rotation mtrx
		G4ThreeVector(0,(30-iii*4)*mm, -34.5*mm),		//at (0,0,0)
	Si_strip_log,					//its logical volume
	"Si strip",					//its name
	telescope_log,		//its mother	volume
	false,			//no boolean operation
	iii,			//copy number
	checkOverlaps);			//overlaps checking
	}

	// CsI strips placement
for (int iii=0; iii<4; iii++)
	{
new G4PVPlacement(0,						//rotation mtrx
		G4ThreeVector(0,(24.75-iii*16.5)*mm, 1), 		//at (0,0,0)
	CsI_strip_log,				//its logical volume
	"CsI strip",				//its name
	telescope_log,		//its mother	volume
	false,			//no boolean operation
	iii,			//copy number
	checkOverlaps);			//overlaps checking
	}



	// Pixels/crystals placement

	// Pixels placement
for (int iii=0; iii<16; iii++)
	{
new G4PVPlacement(	0,						//rotation mtrx
			G4ThreeVector((30-iii*4)*mm,0,0), 		//at (0,0,0)
		Si_pixel_log,				//its logical volume
		"Si pixel",					//its name
		Si_strip_log,					 //its mother	volume
		false,			//no boolean operation
		iii,			//copy number
		checkOverlaps);			//overlaps checking

	}

	// Crystals placement
for (int iii=0; iii<4; iii++)
	{
new G4PVPlacement(	0,					//rotation mtrx
			G4ThreeVector((24.75-iii*16.5)*mm,0,0), 	//at (0,0,0)
		CsI_crystal_log,			//its logical volume
		"CsI crystal",				//its name
		CsI_strip_log,		//its mother	volume
		false,			//no boolean operation
		iii,			//copy number
		checkOverlaps);			//overlaps checking
	}


	//
	//always return the physical World
	//
	return phys_wrld;
}

	void DetectorConstruction::ConstructSDandField()
{

	G4SDManager* SDMan = G4SDManager::GetSDMpointer();
	G4String SDname;

	G4VSensitiveDetector* sensSilicon = new siliconSD(SDname="sensSilicon");
	SDMan->AddNewDetector(sensSilicon);
	Si_pixel_log->SetSensitiveDetector(sensSilicon);

	G4VSensitiveDetector* sensCesium = new cesiumSD(SDname="sensCesium");
	SDMan->AddNewDetector(sensCesium);
	CsI_crystal_log->SetSensitiveDetector(sensCesium);



}
/*
	void DetectorConstruction::SetMaxStep(G4double maxStep)
	{
	fStepLimit->SetMaxAllowedStep(maxStep);
	}	
*/