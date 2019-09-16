#pragma once
#include<vector>
#include<memory>
#include<random>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
//#include "CellList.h"
#include <armadillo>
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
        double u;
        int targetIdx;
        double cost;
        double availControl;
        int nbcount;
        int inlier;
        int marked;
	std::vector<int> landmarkIdx;
        double targetPos[3];
        std::vector<int> nbLandmark;
        std::vector<double> nbLandmarkDist;
        bool targetIsLandmark;
        bool targetIsTarget;
        bool notReached;
        int transporterFlag;
        
        double EudDistToTarget;
	double ShortestPathDistToTarget;
        double energy_accumlator1,energy_accumlator2,energy_accumlator3,energy_accumlator4, eneregy_maintain_accumlator;
        double instant_output_work, instant_conserve_work, instant_input_work;
        double friction_accumulator1, friction_accumulator2,friction_accumulator3, friction_accumulator4;
        double eneregy_useful_maintain_accumlator, v_projection;
        double Fx,Fy, Vx,Vy; // Fx, Fy are the drifting velocity
        double potential;
        double maintainratio;
        double energy_input_transport,energy_useful_transport;
        double energy_vsp_input_maintain, energy_vsp_useful_maintain, energy_vsp_useful_maintain_positive;
        double energy_vsp_input_transport, energy_vsp_useful_transport;
        double energy_instant_vsp_input, energy_instant_vsp_transport, energy_instant_vsp_maintain;
        
        
        
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
    void calForcesHelper_DLAO(double ri[3], double rj[3], double F[3],int i, int j, double& pot);
    void calForcesHelper_DL(double ri[3], double rj[3], double F[3],int i, int j);
    
    bool cellListFlag;
    std::shared_ptr<CellList> cellList, obsCellList;
    int dimP;
    static const double kb, T, vis;
    int numP, numObstacles;
    double radius, radius_nm;
    double LJ,rm,accumTargetMove;
    double Bpp; //2.29 is Bpp/a/kT
    double Kappa; // here is kappa*radius
    double Os_pressure, Os_pressure_origin;
    double L_dep; // 0.2 of radius size, i.e. 200 nm
    double combinedSize;
    double eps[3][3][3];
    std::vector<double> velocity={0.0,5.0e-6,5.0e-6}; // here is for simpication of binary actuation
//    std::vector<double> velocity={0.0, 5.0e-6};
    int numControl;
    particle targetCenter, targetCenter_avg, previousTargetCenter;
    state particles, targets,initialDistToCenter;
    arma::mat targetCenter_history;
    long long targetCenter_historyCounter;
    posArray obstacles; 
    std::vector<int> control;
    std::string iniFile;
    double dt_, cutoff, mobility, diffusivity_r, diffusivity_t;
    std::default_random_engine rand_generator;
    std::shared_ptr<std::normal_distribution<double>> rand_normal;
    int trajOutputInterval;
    long long timeCounter,fileCounter;
    std::ofstream trajOs, opOs, osTarget, osCargo;
    std::string filetag;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> Uniformdistribution{0.0,1.0};
    std::normal_distribution<double> Normaldistribution{0.0,1.0};
    void outputTrajectory(std::ostream& os);
    void outputOrderParameter(std::ostream& os);
    void readxyz(const std::string filename);
    void updateBodyFrameVec();
    void readTarget(std::string filename);
    void readObstacle();
    void getPermutator();
};



