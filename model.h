#pragma once
#include<vector>
#include<memory>
#include<random>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
//#include "CellList.h"

class CellList; 
typedef std::shared_ptr<CellList> cellList_ptr;
class Model {
public:

    struct particle {
        double r[3],F[3],ori_vec[3][3];
        double phi;
        double theta;
        int u;
        int targetIdx;
        double cost;
        int availControl;
        int nbcount;
        int inlier;
        int marked;
	int landmarkIdx;
        double targetPos[3];
        std::vector<int> nbLandmark;
        std::vector<double> nbLandmarkDist;
        double EudDistToTarget;
	double ShortestPathDistToTarget;
        
        particle(double x = 0, double y = 0, double z = 0){
            r[0]=x;r[1]=y;r[2]=z;
        }
    };
    typedef std::shared_ptr<particle> particle_ptr;
    typedef std::vector<particle_ptr> state;
   
    Model(){}
    Model(cellList_ptr cell);
    ~Model() {trajOs.close();
    opOs.close(); osTarget.close();
    }
    void run();
    void run(int steps);
    void createInitialState();
    state getCurrState(){return particles;}
    int getDimP(){return dimP;}
    void setControl(int c);
    state getTargets(){return targets;}
    double dt(){return dt_;}
    int np(){return numP;}
    double calHausdorff();
    double calPsi6();
    double calRg();
private:
    void calForces();
    void calForcesHelper(int i, int j, double F[3]);
    bool cellListFlag;
    std::shared_ptr<CellList> cellList;
    int dimP;
    static const double kb, T, vis;
    int numP;
    double radius, radius_nm;
    double LJ,rm;
    double Bpp; //2.29 is Bpp/a/kT
    double Kappa; // here is kappa*radius
    double Os_pressure;
    double L_dep; // 0.2 of radius size, i.e. 200 nm
    double combinedSize;
    double eps[3][3][3];
    std::vector<double> velocity={0.0,2.0e-6,5.0e-6};
    state particles, targets;
    std::vector<int> control;
    std::string iniFile;
    double dt_, cutoff, mobility, diffusivity_r, diffusivity_t;
    std::default_random_engine rand_generator;
    std::shared_ptr<std::normal_distribution<double>> rand_normal;
    int trajOutputInterval;
    int timeCounter,fileCounter;
    std::ofstream trajOs, opOs, osTarget;
    std::string filetag;
 //   std::vector<Model::particle> targets;
    void outputTrajectory(std::ostream& os);
    void outputOrderParameter(std::ostream& os);
    void readxyz(const std::string filename);
    void updateBodyFrameVec();
    void readTarget(std::string filename);
    void getPermutator();
};



