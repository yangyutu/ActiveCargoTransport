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

    struct pos{
        double r[3];
        pos(double x = 0, double y = 0, double z = 0){
            r[0]=x;r[1]=y;r[2]=z;
        }
    };
    
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
	std::vector<int> landmarkIdx;
        double targetPos[3];
        std::vector<int> nbLandmark;
        std::vector<double> nbLandmarkDist;
        bool targetIsLandmark;
        bool targetIsTarget;
        double EudDistToTarget;
	double ShortestPathDistToTarget;
        
        particle(double x = 0, double y = 0, double z = 0){
            r[0]=x;r[1]=y;r[2]=z;
        }
    };
    typedef std::shared_ptr<particle> particle_ptr;
    typedef std::vector<particle_ptr> state;
    typedef std::vector<std::shared_ptr<pos>> posArray;
   
    Model();
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
    posArray getObstacles(){return obstacles;}
    double dt(){return dt_;}
    int np(){return numP;}
    double calHausdorff();
    double calPsi6();
    double calRg();
    double calEudDeviation();
private:
    void calForces();
    void calForcesHelper(int i, int j, double F[3]);
    void calObsForcesHelper(int i, int j, double F[3]);
    
    bool cellListFlag;
    std::shared_ptr<CellList> cellList, obsCellList;
    int dimP;
    static const double kb, T, vis;
    int numP, numObstacles;
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
    posArray obstacles; 
    std::vector<int> control;
    std::string iniFile;
    double dt_, cutoff, mobility, diffusivity_r, diffusivity_t;
    std::default_random_engine rand_generator;
    std::shared_ptr<std::normal_distribution<double>> rand_normal;
    int trajOutputInterval;
    int timeCounter,fileCounter;
    std::ofstream trajOs, opOs, osTarget;
    std::string filetag;
    void outputTrajectory(std::ostream& os);
    void outputOrderParameter(std::ostream& os);
    void readxyz(const std::string filename);
    void updateBodyFrameVec();
    void readTarget(std::string filename);
    void readObstacle();
    void getPermutator();
};


