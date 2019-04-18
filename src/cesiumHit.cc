#include "cesiumHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"

#include <iomanip>
//class G4AttDef;
//class G4AttValue;
G4Allocator<cesiumHit> cesiumHitAllocator;

cesiumHit::cesiumHit()
 : G4VHit(),
	pixelNo(-1),
	stripNo(-1),
	detectorNo(-1),
	checkNo(-1),
	position(-1),
	momentum(-1),
	energy(-1),
	particle()
{}

cesiumHit::~cesiumHit()
{}

cesiumHit::cesiumHit(const cesiumHit &right)
: G4VHit()
{
stripNo	= right.pixelNo;
stripNo	= right.stripNo;
stripNo	= right.detectorNo;	
stripNo	= right.checkNo;
position	= right.position;
momentum	= right.momentum;
energy	= right.energy;
particle	= right.particle;
}

const cesiumHit& cesiumHit::operator=(const cesiumHit &right)
{
pixelNo	= right.pixelNo;
stripNo	= right.stripNo;
detectorNo	= right.detectorNo;
checkNo	= right.checkNo;
position	= right.position;
momentum	= right.momentum;
energy		= right.energy;
particle	= right.particle;
return *this;
}

G4int cesiumHit::operator==(const cesiumHit &right) const
{
	return ( this == &right ) ? 1 : 0;
}

void cesiumHit::Draw()
{
	//G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	//if(pVVisManager)
	{
	/*
	G4Circle circle(fPos);
	circle.SetScreenSize(4.);
	circle.SetFillStyle(G4Circle::filled);
	G4Colour colour(1.,0.,0.);
	G4VisAttributes attribs(colour);
	circle.SetVisAttributes(attribs);
	pVVisManager->Draw(circle);
	*/
	}
}

void cesiumHit::Print()
{
G4cout<<pixelNo<<G4endl;
	/*G4cout
	<< "	trackID: " << fTrackID << " chamberNb: " << fChamberNb
	<< "Edep: "
	<< std::setw(7) << G4BestUnit(fEdep,"Energy")
	<< " Position: "
	<< std::setw(7) << G4BestUnit( fPos,"Length")
	<< G4endl;*/
}