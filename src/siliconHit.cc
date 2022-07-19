#include "siliconHit.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"

#include <iomanip>
G4Allocator<siliconHit> siliconHitAllocator;

siliconHit::siliconHit()
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

siliconHit::~siliconHit()
{}

siliconHit::siliconHit(const siliconHit & right)
: G4VHit()
{
pixelNo        = right.pixelNo;
stripNo        = right.stripNo;
detectorNo     = right.detectorNo;
stripNo        = right.checkNo;
position       = right.position;
momentum       = right.momentum;
energy         = right.energy;
particle       = right.particle;
position       = right.position;
}

const siliconHit& siliconHit::operator=(const siliconHit& right)
{
pixelNo        = right.pixelNo;
stripNo        = right.stripNo;
detectorNo     = right.detectorNo;
checkNo        = right.checkNo;
position       = right.position;
momentum       = right.momentum;
energy         = right.energy;
particle       = right.particle;
position       = right.position;
  return *this;
}

G4int siliconHit::operator==(const siliconHit &right) const
{
  return ( this == &right ) ? 1 : 0;
}

void siliconHit::Draw()
{
}

void siliconHit::Print()
{
}