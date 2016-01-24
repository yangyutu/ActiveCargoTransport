#pragma once
#include <string>

struct Parameter{
	int N, dim, trajOutputInterval;
	double radius, dt, diffu_t, diffu_r, Bpp, Os_pressure, L_dep, cutoff,kappa;
	int controlStep, equilibrateStep;
	std::string iniConfig, filetag, targetConfig;
        int seed;
        int motionFlag, shapeFlag, cargoTransFlag;
        
        int assignViaEud;
        int landmarkFlag;
        double landmarkLength;
        double landmarkDist;
        double blockCost;
        int blockThresh;
        int obstacleFlag;
        std::string obstacleFilename;
        int assignmentMethod; // 1. single optimal cost based assignment 2. Eud distance based assignment
        // 3. shortest path landmark method based on assignment
};


