#ifndef siliconHit_h
#define siliconHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ParticleDefinition.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and Pos of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class siliconHit : public G4VHit
{
	public:
	siliconHit();
	siliconHit(const siliconHit&);
	virtual ~siliconHit();		

	// operators
	const siliconHit &operator=(const siliconHit &right);
	G4int operator==(const siliconHit &right) const;

	inline void *operator new(size_t);
	inline void operator delete(void* hit);

	// methods from base class
	virtual void Draw();
	virtual void Print();

	// Set methods
	inline void SetPixelNo(G4int pixel)								{ pixelNo=pixel; }
	inline G4int GetPixelNo()											{ return pixelNo; }

	inline void SetStripNo(G4int strip)								{ stripNo=strip; }
	inline G4int GetStripNo()											{ return stripNo; }

	inline void SetDetectorNo(G4int detector)						{ detectorNo=detector; }
	inline G4int GetDetectorNo()										{ return detectorNo; }

	inline void SetCheckNo(G4int check)								{ checkNo=check; }
	inline G4int GetCheckNo()											{ return checkNo; }

	inline void SetPos(G4ThreeVector pos)							{ position=pos; }
	inline G4ThreeVector GetPos()										{ return position; }

	inline void SetMomentum(G4ThreeVector mom)					{ momentum = mom; }
	inline G4ThreeVector GetMomentum()								{ return momentum; }

	inline void SetEnergy(G4double ene)								{ energy = ene; }
	inline G4double GetEnergy()										{ return energy; }
	
	inline void SetParticle(G4ParticleDefinition* pdef)		{ particle = pdef; }
	inline G4ParticleDefinition* GetParticle()					{ return particle; }

	void add(const siliconHit &right);
//stripNo;

//Pos;
//momentum;
//energy;
//particle;	
private:
	G4int	pixelNo;
	G4int	stripNo;
	G4int	detectorNo;
	G4int	checkNo;
	G4ThreeVector position;
	G4ThreeVector momentum;
	G4double	energy;
	G4ParticleDefinition* particle;
};

typedef G4THitsCollection<siliconHit> siliconHitsCollection;
extern G4Allocator<siliconHit> siliconHitAllocator;

inline void *siliconHit::operator new(size_t)
{
	void *hit;
	hit = (void *) siliconHitAllocator.MallocSingle();
	return hit;
}

inline void siliconHit::operator delete(void *hit)
{
	siliconHitAllocator.FreeSingle((siliconHit*) hit);
}

#endif