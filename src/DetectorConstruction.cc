#include "DetectorConstruction.hh"

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
	geometryID = cs::runNo/10;
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{

		if (geometryID==1)
		{
			deut_angle		= (65.0 + 0.3) *deg;
			helium_angle	= (15.0)*deg;
			target_angle	= -45.0*deg;
			tarPosition		= 10.0;
			tarThickness	= (80.0 + cs::tarThicknessShift)*um;
		}

		if (geometryID==2)
		{
			deut_angle		= (50.0) *deg;
			helium_angle	= (15.0)*deg;
			target_angle	= -6.0*deg;
			tarPosition		= 10.0;
			tarThickness	= 2.0 * (80.0 + cs::tarThicknessShift)*um;
		}

		if (geometryID==3)
		{
			deut_angle		= (35.0) *deg;
			helium_angle	= (15.0 - 1.0)*deg;
			target_angle	= 0.0*deg;
			tarPosition		= 10.0;
			tarThickness	= 2.0 * (80.0 + cs::tarThicknessShift)*um;
		}
		
	G4bool checkOverlaps = true;
	// Get nist material manager
	G4NistManager *nist = G4NistManager::Instance();
	G4Material *matSi = nist->FindOrBuildMaterial("G4_Si");
	G4Material *matCsI = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
	G4Material *matWrld = nist->FindOrBuildMaterial("G4_Galactic");
	G4Material *matCD2 = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
	G4Material *matSS = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
	G4Material *matMylar = nist->FindOrBuildMaterial("G4_MYLAR");
	//G4Material *gasDeut = new G4Material ("gasDeut", 1, 2, 0.800507*(gram/mm3), kStateGas, 30*kelvin, 0.5*bar);
	/*(const G4String &name, G4double z, G4double a, G4double density, G4State state=kStateUndefined, G4double temp, G4double pressure)*/
	


	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//xxxxxxxx	DETECTORS PARAMETERS	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	
	
	G4Box* worldBox =	new G4Box(	"World",	//its name
											0.5*world_sizeX,
											0.5*world_sizeY,
											0.5*world_sizeZ);	//its size
	
	G4LogicalVolume* log_wrld = new G4LogicalVolume(	worldBox,		//its solid
																		matWrld,		//its material
																		"log_World");	//its name
	
	G4VPhysicalVolume* phys_wrld = new G4PVPlacement(	0,							//no rotation
																		G4ThreeVector(),		//at (0,0,0)
																		log_wrld,				//its logical volume
																		"phys_World",			//its name
																		0,							//its mother	volume
																		false,					//no boolean operation
																		0,							//copy number
																		checkOverlaps);		//overlaps checking
																		log_wrld->SetVisAttributes(new G4VisAttributes(false));
	//
	// Target
	//
	const G4ThreeVector zeroVector(0.0,0.0,0.0);
	G4RotationMatrix *zeroRotMatrix = new G4RotationMatrix();
	const G4ThreeVector vTarPosition = G4ThreeVector{0.0, 0.0, tarPosition};

	G4RotationMatrix* rotMtrx_Target = new G4RotationMatrix();
	rotMtrx_Target->rotateY(target_angle);


	// CD2 target	1
	G4Box *solidTar = new G4Box(	"foil",	//target sizess
											tar_x/2.0,
											tar_y/2.0,
											tarThickness/2.0);
	
	G4LogicalVolume *logTar = new G4LogicalVolume(	solidTar,	//its solid
																	matCD2,		//its material - I can take poly
																	"logTar");	//its name
	
	new G4PVPlacement(rotMtrx_Target,	//no rotation
						vTarPosition,		//at position
						logTar,				//its logical volume
						"phys_Target",		//its name
						log_wrld,			//its mother	volume
						false,				//no boolean operation
						0,						//copy number
						checkOverlaps);	//overlaps checking



//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxx	BOXES - SOLIDS	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	// Det box
	const G4double boxTele_X =	35*mm, boxTele_Y = 35*mm, boxTele_Z = 35*mm;
	G4Box *boxTelescope = new G4Box(	"boxTelescope",
												boxTele_X,
												boxTele_Y,
												boxTele_Z);

	// silicon strip box 
	const G4double Si_strip_x =	0.5 * 32 * cs::widthStripX*mm, Si_strip_y = 0.5 * cs::widthStripY*mm, Si_strip_z = 0.5*mm;
	G4Box *box_stripSi = new G4Box(	"Si strip",		//name
												Si_strip_x,		//X
												Si_strip_y,		//Y
												Si_strip_z);	//Z


	// CsI strip box 
	const G4double rowCsI_x =	34*mm, rowCsI_y = 8.25*mm, rowCsI_z = 34*mm;
	G4Box* box_rowCsI = new G4Box("CsI strip",	//name
											rowCsI_x,		//X
											rowCsI_y,		//Y
											rowCsI_z);		//Z


	//pixel of Si box
	const G4double pixelSi_x = 0.5*cs::widthStripX*mm, pixelSi_y = 0.5*cs::widthStripY*mm, pixelSi_z = 0.5*mm;
	G4Box* box_pixelSi = new G4Box(	"pixeSi",		//name
												pixelSi_x,			//X
												pixelSi_y,			//Y
												pixelSi_z);	 		//Z


	// Crystal box		//x and y are actually 8.25
	const G4double crystalCsI_x =	8.25*mm, crystalCsI_y = 8.25*mm, crystalCsI_z = 34*mm;
	G4Box* box_crystalCsI = new G4Box(	"crystalCsI", 
													crystalCsI_x,
													crystalCsI_y,
													crystalCsI_z);


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxx	LOGICAL VOLUMES	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	// Telescope LogVolume
	G4LogicalVolume *logTelescope = new G4LogicalVolume(	boxTelescope,		//its solid
																			matWrld,				//its material - I can take poly
																			"telescope_log");	//its name

	// CsI strip LogVolume
	G4LogicalVolume *log_rowCsI = new G4LogicalVolume(		box_rowCsI,			//its solid
																			matWrld,				//its material - I can take poly
																			"log_rowCsI");		//its name

	// Silicon strip LogVolume
	G4LogicalVolume *log_stripSi = new G4LogicalVolume(	box_stripSi,		//its solid
																			matWrld,				//its material - I can take poly
																			"log_stripSi");	//its name

	//my sensitive detectors
	// Crystal LogVolume
	log_crystalCsI = new G4LogicalVolume	(box_crystalCsI,		//its solid
														matCsI,					//its material - I can take poly
														"crystalCsI_log");	//its name

	// Pixel LogVolume
	log_pixelSi = new G4LogicalVolume(	box_pixelSi,		//its solid
													matSi,				//its material
													"pixelSi_log");	//its name


	// Visibility

logTelescope->SetVisAttributes(G4VisAttributes::GetInvisible());
log_stripSi->SetVisAttributes(G4VisAttributes::GetInvisible());
log_rowCsI->SetVisAttributes(G4VisAttributes::GetInvisible());
//crystalCsI_log-> SetVisAttributes (G4VisAttributes::GetInvisible());
//pixelSi_log-> SetVisAttributes (G4VisAttributes::GetInvisible());


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxx	DETECTORS POSITIONING	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//Deuterium telescope container positioning
//For 5th geometry (with gas target there was detector shift by 6.7*mm at d=232 mm from rotation axis.
//Rotation axis is at distance of 132*mm from target (0.0, 0.0, 0.0) at 9.0*deg
sqlang = deut_angle;
sqrang = helium_angle;

vPosition2HTelescope = G4ThreeVector(	(sin(sqlang) * (boxTele_Z+sql_dist))* mm,
										cs::widthStripX,
										(cos(sqlang) * (boxTele_Z+sql_dist))* mm);

//Helium telescope container positioning
vPosition6HeTelescope = G4ThreeVector(	(-sin(sqrang) * (sqr_dist+boxTele_Z))* mm,
										0.0,
										(cos(sqrang) * (sqr_dist+boxTele_Z))* mm);

//Deuterium telescope virtContainer rotation
G4RotationMatrix *deutTeleRotMtrx = new G4RotationMatrix();
deutTeleRotMtrx->rotateX(0.0);
deutTeleRotMtrx->rotateY(-sqlang);
deutTeleRotMtrx->rotateZ(0.0); 


//Helium telescope virtContainer rotation
G4RotationMatrix *heTeleRotMtrx = new G4RotationMatrix();
heTeleRotMtrx->rotateX(0.0);
heTeleRotMtrx->rotateY(sqrang);
heTeleRotMtrx->rotateZ(0.0);	


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxx	DETECTORS PLACEMENT	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	// Telescopes placement

	// Deuterium telescope placement	
	new G4PVPlacement(	deutTeleRotMtrx,				//rotation mtrx
								vPosition2HTelescope,		//vector of center of the detector
								logTelescope,					//its logical volume
								"Deuterium telescope",		//its name
								log_wrld,						//its mother	volume
								false,							//no boolean operation
								0,									//copy number - for deuterium
								checkOverlaps);				//overlaps checking


	// Helium telescope placement
	new G4PVPlacement(	heTeleRotMtrx,					//rotation mtrx
								vPosition6HeTelescope,		//vector of center of the detector
								logTelescope,					//its logical volume
								"helium telescope",			//its name
								log_wrld,						//its mother	volume
								false,							//no boolean operation
								1,									//copy number - for helium
								checkOverlaps);				//overlaps checking

	// Strips placement

	// Silicon strips placement
for (int iii=0; iii<16; iii++)
{
	new G4PVPlacement(0,														//rotation mtrx
							G4ThreeVector(0,(cs::sqlYzero+iii*cs::widthStripY)*mm, -34.5*mm),	//at (0,0,0)
							log_stripSi,										//its logical volume
							"Si strip",											//its name
							logTelescope,										//its mother	volume
							false,												//no boolean operation
							iii,													//copy number
							checkOverlaps);									//overlaps checking
}

	// CsI strips placement
for (int iii=0; iii<4; iii++)
{
	new G4PVPlacement(0,														//rotation mtrx
							G4ThreeVector(0,(24.75-iii*16.5)*mm, 1), 	//at (0,0,0)
							log_rowCsI,											//its logical volume
							"CsI strip",										//its name
							logTelescope,										//its mother	volume
							false,												//no boolean operation
							iii,													//copy number
							checkOverlaps);									//overlaps checking
}



	// Pixels/crystals placement

	// Pixels placement
for (int iii=0; iii<32; iii++)
{
	new G4PVPlacement(0,															//rotation mtrx
							G4ThreeVector((cs::sqlXzero-iii*cs::widthStripX)*mm,0,0),	//at (0,0,0)
							log_pixelSi,											//its logical volume
							"Si pixel",												//its name
							log_stripSi,					 						//its mother	volume
							false,													//no boolean operation
							iii,														//copy number
							checkOverlaps);										//overlaps checking
}

	// Crystals placement
for (int iii=0; iii<4; iii++)
{
	new G4PVPlacement(0,														//rotation mtrx
							G4ThreeVector((24.75-iii*16.5)*mm,0,0), 	//at (0,0,0)
							log_crystalCsI,									//its logical volume
							"CsI crystal",										//its name
							log_rowCsI,											//its mother	volume
							false,												//no boolean operation
							iii,													//copy number
							checkOverlaps);									//overlaps checking
}


	//
	//always return the physical World
	//apparently

	return phys_wrld;
}

void DetectorConstruction::ConstructSDandField()
{
	G4SDManager *SensitiveDetectorManager = G4SDManager::GetSDMpointer();

	G4VSensitiveDetector *sensitiveDetectorSilicon = new siliconSD("sensitiveDetectorSi");
	SensitiveDetectorManager->AddNewDetector(sensitiveDetectorSilicon);
	log_pixelSi->SetSensitiveDetector(sensitiveDetectorSilicon);

	G4VSensitiveDetector *sensitiveDetectorCesium = new cesiumSD("sensisiveDetectorCsI");
	SensitiveDetectorManager->AddNewDetector(sensitiveDetectorCesium);
	log_crystalCsI->SetSensitiveDetector(sensitiveDetectorCesium);
}
/*
	void DetectorConstruction::SetMaxStep(G4double maxStep)
	{
	fStepLimit->SetMaxAllowedStep(maxStep);
	}	
*/
