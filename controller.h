#pragma once
#include<vector>
#include<map>
#include <random>
#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <set>
#include "common.h"
#include "boost/multi_array.hpp"
#include "model.h"
#include "lemon/list_graph.h"
#include <array>
#include <unordered_map>
class Controller {
public:

    Controller(Model::state s, Model::state targets, Model::posArray obstacles0);
    ~Controller() {
    }
    typedef std::vector<int> control;
    enum obstacleType{
      staticObs, dynamicObs
    };
    double calAssignment(Model::state s, Model::state targets,int dimP);
    void calControl(Model::state s, Model::state targets, int dimP);
    

    void getErrorDist();
    void translate_2d(double phi,Model::state s);
    void rotate_2d(Model::state s);
    void alignTarget_t(Model::state s, Model::state targets);
//    void alignTarget_r(Model::state s, Model::state targets);
    void alignTarget_rt(Model::state s, Model::state targets);
    void translateCargo_2d(double phi, Model::state s);
    void translateCargoFollowPath_2d(Model::state s);
    void alignCargo(Model::state s,Model::state t);
    void constructNotReachedSet(Model::state s);
    double calSeqAssignment(Model::state s, Model::state targets, int expand);
    double calAssignmentSeqViaShortestPath(Model::state s,Model::state targets, int expand);
 
    void buildTargetGraph();
    int getNumMarked(){return numMarked;}
    double getDeviation(){return deviation;}
    void readVelocityMap(std::string filename);
//    void register_2d(Model::state s, Model::state targets);
    struct ColliInfo{
        double vSet[3];
        double colliThresh;
    };
    
private:
    
    lemon::ListGraph targetG, landmarkG;
    std::vector<lemon::ListGraph::Node> nodes_t, nodes_l; 
    std::vector<lemon::ListGraph::Edge> edges_t, edges_l;
    std::shared_ptr<lemon::ListGraph::NodeMap<int>> marked_t,index_t;
    //std::shared_ptr<lemon::ListGraph::NodeMap<Model::particle_ptr>> landmark_pos;
    std::vector<Model::pos> landmarkPos;
    //std::shared_ptr<lemon::ListGraph::NodeMap<std::array<double,3>>> landmark_pos2;
	std::shared_ptr<lemon::ListGraph::EdgeMap<double>> internalLength, length;
    int numMarked,numSurface, landmarkLength, numLandmark, numTargets;
    double landmarkDist;
    double deviation;
    std::vector<int> markedIdx;
    ColliInfo colliInfo;
    std::random_device rd;
    typedef boost::multi_array<double, 2> Array2D_type;
    std::shared_ptr<Array2D_type> maps[3];
    
//    std::uniform_int_distribution<> dis(0,1);
    void calAvoidance2d(Model::state s);
    int dimP, numP;
    double radius;
    int x_binNum,y_binNum;
    double del_x,del_y,del_z,min_x,min_y,min_z;
    void readErrorMap();
    void readPolicyMap();
    void readTargets();
    void calAvoidance2d_simpleCollision(Model::state s);
    void calAvoidance2d_simpleCollision_continousV(Model::state s);
    double getCostFromMap(double x,double y, double z, int mapIndex);
    void calWeightCenter(Model::state s, double center[3],int flag);
    void calInlier(Model::state s);
    void readCostMap();
    void calEudDist(Model::state s);
    void expandTargets();

    void readObtacle();
    void constructLandmark();
    void assignLandmarkIdx(Model::state s, double scale);
    double calExtraCost(Model::state s, double r1[3],double scale1, double r2[3], double scale2);
    void calShortestPathDistBetweenLandmarks(Model::state s);
    void calShortestPathDistBetweenST(Model::state s, Model::state targets);
    void calShortestPathHelper(Model::state s, Model::state targets);
  
    std::vector<std::vector<double>> shortestPathDistLandmarkMat, shortestPathDistSTMat;
    std::vector<Model::particle> landmarks;
    std::vector<int> assignment;
    std::vector<double> availControl;
    std::set<int> notReachedSet;

    Model::state targets_, s_;
    double blockCost;
    
    void constructObstacles();
    void constructDynamicObstacles(Model::state s);
    bool isOverlapObstacle(int x, int y);
    bool isOverlapDynamicObstacle(int x, int y);
    bool isOverlapObstacle(double x, double y);
    bool isPathIntersectObstacle(double x, double y, double newx, double newy, Controller::obstacleType obsType);
    Model::posArray obstacles;
    std::unordered_set<CoorPair,CoorPairHash,CoorPairEqual> obstacleSet;
    std::unordered_set<CoorPair,CoorPairHash,CoorPairEqual> dynamicObstacleSet;
    
    
    double calAssignment2d(Model::state s, Model::state targets);
    double calAssignmentVisEudCost(Model:: state s, Model::state targets);
    double calAssignmentViaShortestPath(Model::state s,Model::state targets);
   void calControl3d(Model::state s, Model::state targets);
    double calAssignment3d(Model::state s, Model::state targets);
    void calControl2d(Model::state s, Model::state targets);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> Uniformdistribution{0.0,1.0};
    std::normal_distribution<double> Normaldistribution{0.0,1.0};
    // read velocity map
    std::unordered_map<CoorPair,double,CoorPairHash,CoorPairEqual> velocityMap;
    
};
