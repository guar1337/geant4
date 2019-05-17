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
	G4NistManager *nist = G4NistManager::Instance();
	G4Material *matSi = nist->FindOrBuildMaterial("G4_Si");
	G4Material *matCsI = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
	G4Material *matWrld = nist->FindOrBuildMaterial("G4_Galactic");
	G4Material *matCD2 = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
	G4Material *matSS = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
	G4Material *matMylar = nist->FindOrBuildMaterial("G4_MYLAR");

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
	
	G4ThreeVector vTarPosition = G4ThreeVector{0.0, 0.0, 0.0};
	G4RotationMatrix* rotMtrx_Target = new G4RotationMatrix();
	rotMtrx_Target->rotateX(0.*deg);
	rotMtrx_Target->rotateY(-target_angle);
	rotMtrx_Target->rotateZ(0.*deg); 
if (!gasTarget)
{
	// CD2 target	
	G4Box *solidTar = new G4Box(	"foil",	//target sizess
											tar_x/2.0,
											tar_y/2.0,
											tar_z/2.0);
	
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

}
	// Gas target	
/*
1)		Deuterium
1.1)		Deuterium disc
1.2)		+Z part of sphere
1.3)		-Z part of sphere
2)		Exit flange
3)		Stainless Steel 6mkm foil (sphere)


*/
else
{
	G4Tubs *tarContainer = new G4Tubs("whole target container", 0.0*mm, 40.0*mm,			//inner & outer radius
																					25.0*mm,						//length
																					0.0*rad, CLHEP::twopi);	//starting & ending angle)
	G4ThreeVector zeroVector(0.0,0.0,0.0);
	G4RotationMatrix *zeroRotMatrix = new G4RotationMatrix();
	G4Tubs *deutDiscTube = new G4Tubs("deutDiscTube", 	0.0*mm, 12.5*mm,			//inner & outer radius
																		2.0*mm,						//length
																		0.0*rad, CLHEP::twopi);	//starting & ending angle

	G4Sphere *deutSphere = new G4Sphere("deutSphere", 	0.0*mm, 78.63*mm,			//inner & outer radius
																		0.0*rad, 4*CLHEP::pi,	//starting & ending phi
																		0.0*rad, 12.95*2*deg);	//starting & ending theta

	G4Box *sphereCutoff = new G4Box("sphereCutoff", 100.0*mm, 100.0*mm, 77.63*mm);
	G4VSolid *deutCap = new G4SubtractionSolid("deutSphere-sphereCutoff", deutSphere, sphereCutoff);
	G4ThreeVector tempUnionDeutVector(0.0,0.0,(-76.63+1.0)*mm);
	G4VSolid *tempDeuterUnion = new G4UnionSolid("deutDiscTube+deutCap", deutDiscTube, deutCap, zeroRotMatrix, tempUnionDeutVector);
	G4RotationMatrix *deutCapRot = new G4RotationMatrix();
	deutCapRot->rotateY(CLHEP::pi);
	G4ThreeVector secondCapShiftVect(0.0,0.0,(76.63-1.0)*mm);

	G4VSolid *gasCellSolid = new G4UnionSolid("deutCapS+deutDiscTube", tempDeuterUnion, deutCap, deutCapRot, secondCapShiftVect);
	//There is deuterium gas solid (disc + 2 "caps")

	G4Tubs *tarFlangeTube = new G4Tubs("tarFlangeTube", 	12.5*mm, 26*mm,			//inner & outer radius
																			3.5*mm,						//halflength
																			0.0*rad, CLHEP::twopi);	//starting & ending angle

	G4Cons *tarFlangeCutoff = new G4Cons("tarFlangeCutoff",	0.0*mm, 12.5*mm,				//inner & outer radius1 
																				0.0*mm, 18.5*mm,				//inner & outer radius2 
																				3.0*mm,							//halflength
																				0.0*rad, CLHEP::twopi);		//starting & ending angle
	
	
	G4ThreeVector cutoffConeShift(0.0,0.0,0.5*mm);
	G4VSolid *tarFlangeSolid = new G4SubtractionSolid("tarFlangeTube-tarFlangeCutoff", tarFlangeTube, tarFlangeCutoff, zeroRotMatrix, cutoffConeShift);
	//Target flange

	G4Sphere *stainlessSphere = new G4Sphere("stainlessSphere", 	78.63*mm, 78.636*mm,			//inner & outer radius
																						0.0*rad, 4*CLHEP::pi,	//starting & ending phi
																						0.0*rad, 12.95*2*deg);	//starting & ending theta

	G4Box *sphereSSCutoff = new G4Box("sphereCutoff", 100.0*mm, 100.0*mm, 77.63*mm);
	G4VSolid *SSCapSolid = new G4SubtractionSolid("stainlessSphere-sphereSSCutoff", stainlessSphere, sphereSSCutoff);
	//##########################################################################################################
	//#################################  DEUTERIUM TARGET COVER  ###############################################
	//##########################################################################################################
	G4Tubs *tarCoverBody = new G4Tubs("tarCoverBody", 		24.0*mm, 39.5*mm,			//inner & outer radius
																			20.0*mm,						//halflength
																			0.0*rad, CLHEP::twopi);	//starting & ending angle

	G4Tubs *tarCoverBodyCutoff = new G4Tubs("tarCoverBody", 	0.0*mm, 35.5*mm,			//inner & outer radius
																				12.0*mm,						//halflength
																				0.0*rad, CLHEP::twopi);	//starting & ending angle
	G4Tubs *tarCoverBodyCutoffShort = new G4Tubs("tarCoverBody", 	0.0*mm, 27.5*mm,			//inner & outer radius
																						1.0*mm,						//halflength
																						0.0*rad, CLHEP::twopi);	//starting & ending angle

	G4VSolid *tarCoverShell = new G4SubtractionSolid("tarCoverBody-tarCoverBodyCutoff", tarCoverBody, tarCoverBodyCutoff, zeroRotMatrix, zeroVector);
	G4Cons *tarCoverSmallerCone = new G4Cons("tarCoverSmallerCone",	0.0*mm, 24.0*mm,				//inner & outer radius1 
																							0.0*mm, 27.0*mm,				//inner & outer radius2 
																							1.5*mm,							//halflength
																							0.0*rad, CLHEP::twopi);		//starting & ending angle

	G4Cons *tarCoverBiggerCone = new G4Cons("tarCoverBiggerCone",	0.0*mm, 27.5*mm,				//inner & outer radius1
																						0.0*mm, 30.5*mm,				//inner & outer radius2 
																						1.5*mm,							//halflength
																						0.0*rad, CLHEP::twopi);		//starting & ending angle
	
	G4ThreeVector coverConeShift1(0.0, 0.0, 14.5);
	G4ThreeVector coverTubeShift(0.0, 0.0, 17.0);
	G4ThreeVector coverConeShift2(0.0, 0.0, 18.5);
	
	G4VSolid *tarCoverMinusCone1 = new G4SubtractionSolid("tarCoverShell-tarCoverSmallerCone", tarCoverShell, tarCoverSmallerCone, zeroRotMatrix, coverConeShift1);
	G4VSolid *tarCoverMinusTube = new G4SubtractionSolid("tarCoverMinusCone1-tarCoverBodyCutoffShort", tarCoverMinusCone1, tarCoverBodyCutoffShort, zeroRotMatrix, coverTubeShift);
	G4VSolid *tarCoverSolid = new G4SubtractionSolid("tarCoverMinusTube-tarCoverBiggerCone", tarCoverMinusTube, tarCoverBiggerCone, zeroRotMatrix, coverConeShift2);
	//crosssection of target cover
	/*
	G4Box *makeCrosssection = new G4Box("makeCrosssection", 100.0*mm, 100.0*mm, 100.0*mm);
	G4ThreeVector makeCrosssectionVect1(0.0,100.5,0.0);
	G4ThreeVector makeCrosssectionVect2(0.0,-100.5,0.0);
	G4VSolid *tarCoverMinusConeCrossSect1 = new G4SubtractionSolid("tarCoverMinusCone2-makeCrosssection1", tarCoverSolid, makeCrosssection, zeroRotMatrix, makeCrosssectionVect1);
	G4VSolid *tarCoverMinusConeCrossSect2 = new G4SubtractionSolid("tarCoverMinusCone2-makeCrosssection2", tarCoverMinusConeCrossSect1, makeCrosssection, zeroRotMatrix, makeCrosssectionVect2);
	*/
	G4Tubs *mylarFoilSolid = new G4Tubs("mylarFoil",
													0.0*mm, 27.5*mm,			//inner & outer radius
													3.5*um,						//halflength
													0.0*rad, CLHEP::twopi);
	
	G4LogicalVolume *targetLog = new G4LogicalVolume(	tarContainer,	//its solid
																		matWrld,		//its material - I can take poly
																		"mylarFoilLog");	//its name

	G4LogicalVolume *gasCellLog = new G4LogicalVolume(	gasCellSolid,	//its solid
																	matCD2,		//its material - I can take poly
																	"gasCellLog");	//its name

	G4LogicalVolume *tarFlangeLog = new G4LogicalVolume(tarFlangeSolid,	//its solid
																		matSS,		//its material - I can take poly
																		"tarFlangeLog");	//its name

	G4LogicalVolume *SScapLog = new G4LogicalVolume(	SSCapSolid,	//its solid
																		matSS,		//its material - I can take poly
																		"SScapLog");	//its name

	G4LogicalVolume *tarCoverLog = new G4LogicalVolume(tarCoverSolid,	//its solid
																		matSS,		//its material - I can take poly
																		"tarCoverLog");	//its name

	G4LogicalVolume *mylarFoilLog = new G4LogicalVolume(mylarFoilSolid,	//its solid
																		matMylar,		//its material - I can take poly
																		"mylarFoilLog");	//its name

	//logTar->SetVisAttributes(new G4VisAttributes(false));
	targetLog->SetVisAttributes(new G4VisAttributes(false));

	new G4PVPlacement(rotMtrx_Target,	//no rotation
							vTarPosition,		//at position
							targetLog,			//its logical volume
							"Target container",		//its name
							log_wrld,			//its mother	volume
							false,				//no boolean operation
							0,						//copy number
							checkOverlaps);	//overlaps checking

	new G4PVPlacement(zeroRotMatrix,	//no rotation
							vTarPosition,		//at position
							gasCellLog,				//its logical volume
							"Gas volume",		//its name
							targetLog,			//its mother	volume
							false,				//no boolean operation
							0,						//copy number
							checkOverlaps);	//overlaps checking

	new G4PVPlacement(zeroRotMatrix,	//no rotation
							G4ThreeVector(0.0,0.0,-75.63),		//at position
							SScapLog,				//its logical volume
							"SS foil",		//its name
							targetLog,			//its mother	volume
							false,				//no boolean operation
							0,						//copy number
							checkOverlaps);	//overlaps checking

	new G4PVPlacement(zeroRotMatrix,	//no rotation
							G4ThreeVector(0.0,0.0,5.5),		//at position
							tarFlangeLog,				//its logical volume
							"Target flange",		//its name
							targetLog,			//its mother	volume
							false,				//no boolean operation
							0,						//copy number
							checkOverlaps);	//overlaps checking

	new G4PVPlacement(zeroRotMatrix,	//no rotation
							G4ThreeVector(0.0,0.0,0.0),		//at position
							tarCoverLog,				//its logical volume
							"Target cover",		//its name
							targetLog,			//its mother	volume
							false,				//no boolean operation
							0,						//copy number
							checkOverlaps);	//overlaps checking

	new G4PVPlacement(zeroRotMatrix,	//no rotation
							G4ThreeVector(0.0,0.0,16.0),		//at position
							mylarFoilLog,				//its logical volume
							"mylar foil 3.5um",		//its name
							targetLog,			//its mother	volume
							false,				//no boolean operation
							0,						//copy number
							checkOverlaps);	//overlaps checking
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxx	BOXES - SOLIDS	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	// Det box
	const G4double boxTelescope_x =	35*mm, boxTelescope_y = 35*mm, boxTelescope_z = 35*mm;
	G4Box *boxTelescope = new G4Box(	"boxTelescope",
												boxTelescope_x,
												boxTelescope_y,
												boxTelescope_z);

	// silicon strip box 
	const G4double Si_strip_x =	29*mm, Si_strip_y = 1.8125*mm, Si_strip_z = 0.5*mm;
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
	const G4double pixelSi_x = 1.8125*mm, pixelSi_y = 3.625*mm, pixelSi_z = 0.5*mm;
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
//xxxxxxxxxxxx	DETECTORS POSITIONING	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//Deuterium telescope container positioning
G4ThreeVector vPosition_2H_telescope = G4ThreeVector(	sin(deut_angle)*(boxTelescope_x+sql_dist)*mm,
																		0,
																		cos(deut_angle)*(boxTelescope_x+sql_dist)*mm);

//Helium telescope container positioning
G4ThreeVector vPosition_6He_telescope = G4ThreeVector(	-sin(helium_angle)*(boxTelescope_x+sqr_dist)*mm,
																			0,
																			cos(helium_angle)*(boxTelescope_x+sqr_dist)*mm);


//Deuterium telescope container rotation
G4RotationMatrix *deut_tele_rot_mtrx = new G4RotationMatrix();
deut_tele_rot_mtrx->rotateX(0.0);
deut_tele_rot_mtrx->rotateY(-deut_angle);
deut_tele_rot_mtrx->rotateZ(0.0); 


//Helium telescope container rotation
G4RotationMatrix *he_tele_rot_mtrx = new G4RotationMatrix();
he_tele_rot_mtrx->rotateX(0.0);
he_tele_rot_mtrx->rotateY(helium_angle);
he_tele_rot_mtrx->rotateZ(0.0);	

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxx	DETECTORS PLACEMENT	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	// Telescopes placement
	// Deuterium telescope placement	
	new G4PVPlacement(	deut_tele_rot_mtrx,			//rotation mtrx
								vPosition_2H_telescope,		//vector of center of the detector
								logTelescope,					//its logical volume
								"helium telescope",			//its name
								log_wrld,						//its mother	volume
								false,							//no boolean operation
								0,									//copy number - for deuterium
								checkOverlaps);				//overlaps checking


	// Helium telescope placement
	new G4PVPlacement(	he_tele_rot_mtrx,				//rotation mtrx
								vPosition_6He_telescope,	//vector of center of the detector
								logTelescope,					//its logical volume
								"deuterium telescope",		//its name
								log_wrld,						//its mother	volume
								false,							//no boolean operation
								1,									//copy number - for helium
								checkOverlaps);				//overlaps checking

	// Strips placement

	// Silicon strips placement
for (int iii=0; iii<16; iii++)
{
	new G4PVPlacement(	0,														//rotation mtrx
								G4ThreeVector(0,(-27.1875+iii*3.625)*mm, -34.5*mm),	//at (0,0,0)
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
	new G4PVPlacement(	0,														//rotation mtrx
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
	new G4PVPlacement(	0,												//rotation mtrx
								G4ThreeVector((28.09375-iii*1.8125)*mm,0,0),	//at (0,0,0)
								log_pixelSi,								//its logical volume
								"Si pixel",									//its name
								log_stripSi,					 			//its mother	volume
								false,										//no boolean operation
								iii,											//copy number
								checkOverlaps);							//overlaps checking
}

	// Crystals placement
for (int iii=0; iii<4; iii++)
{
	new G4PVPlacement(	0,														//rotation mtrx
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