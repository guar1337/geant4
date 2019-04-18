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
<<<<<<< HEAD
inline void Set_LV_Beam(G4LorentzVector lvBeam)
{
	ParticleInfo_LV_Beam=lvBeam;
}
inline G4LorentzVector Get_LV_Beam()
{
	return ParticleInfo_LV_Beam;
}

inline void Set_LV_2H_CM(G4LorentzVector lv2H_CM)
{
	ParticleInfo_LV_2H_CM = lv2H_CM;
}
inline G4LorentzVector Get_LV_2H_CM()
{
	return ParticleInfo_LV_2H_CM;
}

inline void Set_LV_6He_CM(G4LorentzVector lv6He_CM)
{
	ParticleInfo_LV_6He_CM=lv6He_CM;
}
inline G4LorentzVector Get_LV_6He_CM()
{
	return ParticleInfo_LV_6He_CM;
}
/*
inline void Set_LV_Beam(G4LorentzVector lvBeam);
inline void Set_LV_2H_CM(G4LorentzVector lv2H_CM);
inline void Set_LV_6He_CM(G4LorentzVector lv6He_CM);
=======

inline void Set_LV_Beam(G4LorentzVector lv_Beam);
inline void Set_LV_2H_CM(G4LorentzVector lv2H_CM);
inline void Set_LV_6He_CM(G4LorentzVector lv_6He_CM);
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217

inline G4LorentzVector Get_LV_Beam();
inline G4LorentzVector Get_LV_2H_CM();
inline G4LorentzVector Get_LV_6He_CM();
<<<<<<< HEAD
*/
=======
	
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
private:
G4LorentzVector ParticleInfo_LV_Beam;
G4LorentzVector ParticleInfo_LV_2H_CM;
G4LorentzVector ParticleInfo_LV_6He_CM;
};

#endif