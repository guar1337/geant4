#ifndef cesiumSD_h
#define cesiumSD_h 1

#include "G4VSensitiveDetector.hh"
#include "cesiumHit.hh"
#include <vector>

class G4Step;           //for tracking event
class G4HCofThisEvent;  //Hit collection of This Event
class cesiumSD : public G4VSensitiveDetector       
{
  public:
    cesiumSD(const G4String& name);
    virtual ~cesiumSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* HCofEven);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* HCofEven);

  private:
   cesiumHitsCollection *fHitsCollection;
   G4int fHCID;

};


#endif
