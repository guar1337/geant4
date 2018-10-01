#include "ParticleInfo.hh"

ParticleInfo::ParticleInfo()
{}

ParticleInfo::~ParticleInfo()
{}

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
