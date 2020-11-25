#include "siliconSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
	

siliconSD::siliconSD(const G4String &name)	//det name
 : G4VSensitiveDetector(name),
	fHitsCollection(nullptr),
	fHCID(-1)
{
	collectionName.insert("siliconColl");	//creating hit collection
}
siliconSD::~siliconSD() 
{}
void siliconSD::Initialize(G4HCofThisEvent *HCofEven)
{
	// Create hits collection
	fHitsCollection = new siliconHitsCollection(SensitiveDetectorName, collectionName[0]); 
	if (fHCID<0) 
	{ 
		fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); 
	}
	// Add this collection in HCofEven
HCofEven->AddHitsCollection(fHCID, fHitsCollection); 
}

G4bool siliconSD::ProcessHits(G4Step *step, G4TouchableHistory *)
{	
	siliconHit *newHit = new siliconHit();
	auto edep=step->GetTotalEnergyDeposit();
	if (edep==0.0) return true;
	G4StepPoint *preStep = step->GetPreStepPoint();
	G4TouchableHistory *touchable = (G4TouchableHistory*)(preStep->GetTouchable()); 

	newHit->SetPixelNo(touchable->GetReplicaNumber(0));
	newHit->SetStripNo(touchable->GetReplicaNumber(1));
	newHit->SetDetectorNo(touchable->GetReplicaNumber(2));
	newHit->SetCheckNo(touchable->GetReplicaNumber(3));
	newHit->SetPos(step->GetPreStepPoint()->GetPosition());
	//newHit->SetMomentum( step->GetPreStepPoint()->GetMomentum() );
	newHit->SetEnergy(step->GetTotalEnergyDeposit());
	newHit->SetParticle(step->GetTrack()->GetDefinition() );
	fHitsCollection->insert(newHit);

	return true;
}


void siliconSD::EndOfEvent(G4HCofThisEvent *)
{
 /*
	if ( verboseLevel>1 ) { 
	G4int nofHits = fHitsCollection->entries();
	G4cout << G4endl
	<< "-------->Hits Collection: in this event they are " << nofHits 
	<< " hits in the tracker chambers: " << G4endl;
	//for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
	}
*/
}