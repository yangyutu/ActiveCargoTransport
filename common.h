#pragma once
#include <string>

struct Parameter{
	int N, dim, trajOutputInterval;
	double radius, dt, diffu_t, diffu_r, Bpp, Os_pressure, L_dep, cutoff,kappa;
	int controlStep, equilibrateStep;
	std::string iniConfig, filetag, targetConfig;
        int seed;
        int motionFlag, shapeFlag, cargoTransFlag;
        int particleCellListFlag, obstacleCellListFlag, noControlFlag;
        
        double cellListCutoff, cellListBox_x, cellListBox_y, cellListBox_z;
        int cellListDim, cellListMaxCount;
        int assignViaEud;
        int landmarkFlag;
        double landmarkLength, landmarkMin;
        double landmarkDist;
        double blockCost;
        int blockThresh;
        int obstacleFlag;
        std::string obstacleFilename;
        int assignmentMethod; // 1. single optimal cost based assignment 2. Eud distance based assignment
        // 3. shortest path landmark method based on assignment
};
class CoorPair{
public:
	int x;
	int y;



	CoorPair(){};

	CoorPair(int x0,int y0){x=x0;y=y0;}

};

typedef struct
{
	std::size_t operator() (const CoorPair & CP) const {
		std::size_t h1=std::hash<int>()(CP.x);
		std::size_t h2 = std::hash<int>()(CP.y);
		return h1^(h2<<1);
	}
}CoorPairHash;

typedef struct
{
	bool operator() (const CoorPair & CP1,const CoorPair & CP2) const {
		return (CP1.x==CP2.x)&&(CP1.y==CP2.y);
	}
}CoorPairEqual;

