#include "ParticleInfo.hh"

<<<<<<< HEAD
ParticleInfo::ParticleInfo():
   ParticleInfo_LV_Beam(0),
   ParticleInfo_LV_2H_CM(0),
	ParticleInfo_LV_6He_CM(0)
=======
ParticleInfo::ParticleInfo()
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
{}

ParticleInfo::~ParticleInfo()
{}

<<<<<<< HEAD
void ParticleInfo::Print() const 
{
}
=======
inline void ParticleInfo::Set_LV_Beam(G4LorentzVector lv_Beam)
{
	ParticleInfo_LV_Beam = lv_Beam;
}

inline void ParticleInfo::Set_LV_2H_CM(G4LorentzVector lv2H_CM)
{
	ParticleInfo_LV_2H_CM = lv2H_CM;
}

inline void ParticleInfo::Set_LV_6He_CM(G4LorentzVector lv_6He_CM)
{
	ParticleInfo_LV_6He_CM = lv_6He_CM;
}

inline G4LorentzVector ParticleInfo::Get_LV_Beam()
{
	return ParticleInfo_LV_Beam;
}

inline G4LorentzVector ParticleInfo::Get_LV_2H_CM()
{
	return ParticleInfo_LV_2H_CM;
}

inline G4LorentzVector ParticleInfo::Get_LV_6He_CM()
{
	return ParticleInfo_LV_6He_CM;
}

void ParticleInfo::Print() const 
{}
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
