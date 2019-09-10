#pragma once
#include <string>

struct Parameter{
	int N, dim, trajOutputInterval;
	double radius, dt, diffu_t, diffu_r, Bpp, Os_pressure, L_dep, cutoff,kappa;
	int controlStep, equilibrateStep,nCycles;
	std::string iniConfig, filetag, targetConfig;
        int seed;
        int motionFlag,  cargoTransFlag; //1 move horizontally 2 vertically 3 following a path
        int shapeFlag; // 1 normal shape formation 2 sequential shape formation
        int particleCellListFlag, obstacleCellListFlag, noControlFlag;
        
        int collective_MoveStep, collective_RestoreStep, CollectiveMoveCycle;
        int cargoCaptureStep;
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
        int dynamicTargetFlag;
        double controlTimeInterval, targetDiffuseRatio,  targetVelocityRatio, assignmentTimeInterval;
        double maxVelocity, velocityChangePoint;
        int selfAvoidanceFlag;
        int binaryVelocityFlag;
        int cargoInteractingFlag;
        double targetCenter[3];
        double setConstantV;
        int setConstantVFlag;
        
        int targetHistoryFlag, targetHistorySaveInterval, targetHistoryLength;
        double targetMoveThresh;
        double totalCost;
        
        // transporter selection criterion
        int transporter_nb_thresh;
        double transporter_dist_thresh, transport_angle_thresh; // angle are in unit of degree
        
        // velocity map
        std::string velocityMapName;
        
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


