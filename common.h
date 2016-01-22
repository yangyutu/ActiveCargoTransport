#pragma once
#include <string>

struct Parameter{
	int N, dim, trajOutputInterval;
	double radius, dt, diffu_t, diffu_r, Bpp, Os_pressure, L_dep, cutoff,kappa;
	int controlStep, equilibrateStep;
	std::string iniConfig, filetag, targetConfig;
        int seed;
        int motionFlag, shapeFlag, cargoTransFlag;
        
        int landmarkFlag;
        double landmarkLength;
        double landmarkDist;
        double blockCost;
};


