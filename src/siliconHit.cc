#include "siliconHit.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"

#include <iomanip>
G4Allocator<siliconHit> siliconHitAllocator;

<<<<<<< HEAD
=======
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
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

<<<<<<< HEAD
siliconHit::~siliconHit()
{}

=======
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

siliconHit::~siliconHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
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

<<<<<<< HEAD
=======

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
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

<<<<<<< HEAD
=======
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
G4int siliconHit::operator==(const siliconHit &right) const
{
  return ( this == &right ) ? 1 : 0;
}

<<<<<<< HEAD
=======
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
void siliconHit::Draw()
{
}

<<<<<<< HEAD
void siliconHit::Print()
{
}
=======
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void siliconHit::Print()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
