#include "cesiumSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"



cesiumSD::cesiumSD(const G4String& name)	 //det name
 : G4VSensitiveDetector(name),
	fHitsCollection(nullptr),
	fHCID(-1)
{
	//G4String HCname;
	collectionName.insert("cesiumColl");	 //creating hit collection
}

cesiumSD::~cesiumSD() 
{}

void cesiumSD::Initialize(G4HCofThisEvent* HCofEven)
{
  // Create hits collection
fHitsCollection = new cesiumHitsCollection(SensitiveDetectorName, collectionName[0]); 
if (fHCID<0) 
  { 
	fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); 
  }
  // Add this collection in HCofEven
HCofEven->AddHitsCollection(fHCID, fHitsCollection); 
}

G4bool cesiumSD::ProcessHits(G4Step* step, G4TouchableHistory*)

{  
  G4StepPoint* preStep = step->GetPreStepPoint();
auto edep=step->GetTotalEnergyDeposit();
  if (edep==0.) return true;
  cesiumHit* newHit = new cesiumHit();


  G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());

  newHit->SetPixelNo(touchable->GetReplicaNumber(0) );
  newHit->SetStripNo(touchable->GetReplicaNumber(1) );
  newHit->SetDetectorNo(touchable->GetReplicaNumber(2) );
  newHit->SetParticle(step->GetTrack()->GetDefinition() );
  newHit->SetEnergy(step->GetTotalEnergyDeposit());
  //G4cout<<(step->GetTotalEnergyDeposit())<<" =Edep"<<G4endl;
  fHitsCollection->insert( newHit );
  return true;

  //newHit->SetPosition(step->GetPreStepPoint()->GetPosition() );
  //newHit->SetMomentum( step->GetPreStepPoint()->GetMomentum() );


}



void cesiumSD::EndOfEvent(G4HCofThisEvent*)
{
/*
  //if ( verboseLevel>1 ) { 
	G4int nofHits = fHitsCollection->entries();
	G4cout << G4endl
	<< "-------->Hits Collection: in this event they are " << nofHits 
	 << " hits in the tracker chambers: " << G4endl;
	//for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
*/
}


