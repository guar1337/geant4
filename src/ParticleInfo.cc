#include "ParticleInfo.hh"

ParticleInfo::ParticleInfo():
   ParticleInfo_LV_Beam(0),
   ParticleInfo_LV_2H_CM(0),
	ParticleInfo_LV_6He_CM(0),
	ParticleInfo_MWPC_1_X(0),
	ParticleInfo_MWPC_2_X(0),
	ParticleInfo_MWPC_1_Y(0),
	ParticleInfo_MWPC_2_Y(0),
	ParticleInfo_nx1(0),
	ParticleInfo_nx2(0),
	ParticleInfo_ny1(0),
	ParticleInfo_ny2(0)
{}

ParticleInfo::~ParticleInfo()
{}

void ParticleInfo::Print() const 
{
}
