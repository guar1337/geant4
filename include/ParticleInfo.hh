#ifndef ParticleInfo_h
#define ParticleInfo_h 1

#include "G4VUserPrimaryParticleInformation.hh"
#include "globals.hh"
#include "G4LorentzVector.hh"
class ParticleInfo : public G4VUserPrimaryParticleInformation {

public:
ParticleInfo();
virtual ~ParticleInfo();
virtual void Print() const;
//Set methods	

inline void Set_LV_Beam(G4LorentzVector lv_Beam);
inline void Set_LV_2H_CM(G4LorentzVector lv2H_CM);
inline void Set_LV_6He_CM(G4LorentzVector lv_6He_CM);

inline G4LorentzVector Get_LV_Beam();
inline G4LorentzVector Get_LV_2H_CM();
inline G4LorentzVector Get_LV_6He_CM();
	
private:
G4LorentzVector ParticleInfo_LV_Beam;
G4LorentzVector ParticleInfo_LV_2H_CM;
G4LorentzVector ParticleInfo_LV_6He_CM;
};

#endif