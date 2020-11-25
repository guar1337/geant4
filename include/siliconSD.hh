#ifndef siliconSD_h
#define siliconSD_h 1

#include "G4VSensitiveDetector.hh"
#include "siliconHit.hh"

class G4Step;           //for tracking event
class G4HCofThisEvent;  //Hit collection of This Event

///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero 
/// energy deposit.

class siliconSD : public G4VSensitiveDetector       
{
public:
	siliconSD(const G4String &name);
	virtual ~siliconSD();
  
	// methods from base class
	virtual void   Initialize(G4HCofThisEvent *HCofEven);
	virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory *history);
	virtual void	EndOfEvent(G4HCofThisEvent *HCofEven);

private:
	siliconHitsCollection *fHitsCollection;
	G4int fHCID;

};
#endif
