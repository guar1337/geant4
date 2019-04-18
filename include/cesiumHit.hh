#ifndef cesiumHit_h
#define cesiumHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ParticleDefinition.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

class cesiumHit : public G4VHit
{
<<<<<<< HEAD
public:
	cesiumHit();
	cesiumHit(const cesiumHit &right);
	~cesiumHit();

	// operators
	const cesiumHit& operator=(const cesiumHit &right);
	G4int operator==(const cesiumHit &right) const;

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
	
	inline void SetPosition(G4ThreeVector pos)					{ position=pos; }
	inline G4ThreeVector GetPosition()								{ return position; }
	
	inline void SetMomentum(G4ThreeVector mom)					{ momentum = mom; }
	inline G4ThreeVector GetMomentum()								{ return momentum; }
	
	inline void SetEnergy(G4double ene)								{ energy = ene; }
	inline G4double GetEnergy()										{ return energy; }
	
	inline void SetParticle(G4ParticleDefinition* pdef)		{ particle = pdef; }
	inline G4ParticleDefinition* GetParticle()					{ return particle; }

	
private:
	G4int	pixelNo;
	G4int	stripNo;
	G4int	detectorNo;
	G4int	checkNo;
	G4ThreeVector position;
	G4ThreeVector momentum;
	G4double	energy;
	G4double	Ekin;
	G4ParticleDefinition *particle;
=======
  public:
    cesiumHit();
    cesiumHit(const cesiumHit &right);
    ~cesiumHit();

    // operators
    const cesiumHit& operator=(const cesiumHit &right);
    G4int operator==(const cesiumHit &right) const;

    inline void *operator new(size_t);
    inline void  operator delete(void* hit);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
      inline void SetPixelNo(G4int pixel)                   { pixelNo=pixel; }
      inline G4int GetPixelNo()                             { return pixelNo; }

      inline void SetStripNo(G4int strip)                   { stripNo=strip; }
      inline G4int GetStripNo()                             { return stripNo; }

      inline void SetDetectorNo(G4int detector)             { detectorNo=detector; }
      inline G4int GetDetectorNo()                          { return detectorNo; }

      inline void SetCheckNo(G4int check)                   { checkNo=check; }
      inline G4int GetCheckNo()                             { return checkNo; }
      
      inline void SetPosition(G4ThreeVector pos)            { position=pos; }
      inline G4ThreeVector GetPosition()                    { return position; }
      
      inline void SetMomentum(G4ThreeVector mom)            { momentum = mom; }
      inline G4ThreeVector GetMomentum()                    { return momentum; }
      
      inline void SetEnergy(G4double ene)                   { energy = ene; }
      inline G4double GetEnergy()                           { return energy; }
      
      inline void SetParticle(G4ParticleDefinition* pdef)   { particle = pdef; }
      inline G4ParticleDefinition* GetParticle()            { return particle; }

    
  private:
      G4int         pixelNo;
      G4int         stripNo;
      G4int         detectorNo;
      G4int         checkNo;
      G4ThreeVector position;
      G4ThreeVector momentum;
      G4double      energy;
      G4double      Ekin;
      G4ParticleDefinition* particle;

>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
};


typedef G4THitsCollection<cesiumHit> cesiumHitsCollection;
<<<<<<< HEAD
extern G4Allocator<cesiumHit> cesiumHitAllocator;

inline void* cesiumHit::operator new(size_t)
{
	void *hit;
	hit = (void *) cesiumHitAllocator.MallocSingle();
	return hit;
}

inline void cesiumHit::operator delete(void *hit)
{
	cesiumHitAllocator.FreeSingle((cesiumHit*) hit);
=======

extern G4Allocator<cesiumHit> cesiumHitAllocator;


inline void* cesiumHit::operator new(size_t)
{
    void *hit;
    hit = (void *) cesiumHitAllocator.MallocSingle();
    return hit;
}


inline void cesiumHit::operator delete(void *hit)
{
    cesiumHitAllocator.FreeSingle((cesiumHit*) hit);
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
}

#endif