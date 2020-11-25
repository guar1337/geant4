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

inline void Set_MWPC_1_X(G4float MWPC_1_X)
{
	ParticleInfo_MWPC_1_X=MWPC_1_X;
}
inline G4float Get_MWPC_1_X()
{
	return ParticleInfo_MWPC_1_X;
}

inline void Set_MWPC_2_X(G4float MWPC_2_X)
{
	ParticleInfo_MWPC_2_X=MWPC_2_X;
}
inline G4float Get_MWPC_2_X()
{
	return ParticleInfo_MWPC_2_X;
}

inline void Set_MWPC_1_Y(G4float MWPC_1_Y)
{
	ParticleInfo_MWPC_1_Y=MWPC_1_Y;
}
inline G4float Get_MWPC_1_Y()
{
	return ParticleInfo_MWPC_1_Y;
}

inline void Set_MWPC_2_Y(G4float MWPC_2_Y)
{
	ParticleInfo_MWPC_2_Y=MWPC_2_Y;
}
inline G4float Get_MWPC_2_Y()
{
	return ParticleInfo_MWPC_2_Y;
}

inline void Set_nx1(G4float nx1)
{
	ParticleInfo_nx1=nx1;
}
inline G4float Get_nx1()
{
	return ParticleInfo_nx1;
}

inline void Set_nx2(G4float nx2)
{
	ParticleInfo_nx2=nx2;
}
inline G4float Get_nx2()
{
	return ParticleInfo_nx2;
}

inline void Set_ny1(G4float ny1)
{
	ParticleInfo_ny1=ny1;
}
inline G4float Get_ny1()
{
	return ParticleInfo_ny1;
}

inline void Set_ny2(G4float ny2)
{
	ParticleInfo_ny2=ny2;
}
inline G4float Get_ny2()
{
	return ParticleInfo_ny2;
}



/*
inline void Set_LV_Beam(G4LorentzVector lvBeam);
inline void Set_LV_2H_CM(G4LorentzVector lv2H_CM);
inline void Set_LV_6He_CM(G4LorentzVector lv6He_CM);

inline G4LorentzVector Get_LV_Beam();
inline G4LorentzVector Get_LV_2H_CM();
inline G4LorentzVector Get_LV_6He_CM();
*/
private:
G4LorentzVector ParticleInfo_LV_Beam;
G4LorentzVector ParticleInfo_LV_2H_CM;
G4LorentzVector ParticleInfo_LV_6He_CM;
G4float ParticleInfo_MWPC_1_X;
G4float ParticleInfo_MWPC_2_X;
G4float ParticleInfo_MWPC_1_Y;
G4float ParticleInfo_MWPC_2_Y;
G4float ParticleInfo_nx1;
G4float ParticleInfo_nx2;
G4float ParticleInfo_ny1;
G4float ParticleInfo_ny2;

};

#endif