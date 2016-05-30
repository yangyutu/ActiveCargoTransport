#include<iostream>
#include "simulator.h"
#include "controller.h"
#include "model.h"
#include "CellList.h"
#include "common.h"

Parameter parameter;

void readParameter();

int main(){

    readParameter();
    int N = 85;
    int dim = 2;
    double radius = 1e-6;
//    cellList_ptr cell(new CellList(3.0*radius,2,10,150.0*radius,150.0*radius,50.0*radius));
    std::shared_ptr<Model> m(new Model());
    std::shared_ptr<Controller> c(new Controller(m->getCurrState(),m->getTargets(),m->getObstacles()));
    Simulator simulator(m,c);
    
    for (int i = 0; i < parameter.nCycles; i++) {
        if (parameter.shapeFlag == 1) {
            simulator.shapeForming();
        } else if (parameter.shapeFlag == 2) {
            simulator.shapeForming_seq();
        } else if (parameter.motionFlag == 1) {
            simulator.translate_2d();
        } else if (parameter.motionFlag == 2) {
            simulator.rotate_2d();
        } else if (parameter.cargoTransFlag == 1) {
            simulator.cargoTransport_2d();
        } else if (parameter.noControlFlag == 1) {
            m->createInitialState();
            m->getCurrState()[0]->u = 2;
            m->run(parameter.controlStep);
        }
    }
    return 0;
}

void readParameter(){
    std::string line;
    std::ifstream runfile;
    runfile.open("run.txt");
    getline(runfile, line);
    runfile >> parameter.N;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.dim;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.radius;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.nCycles;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.controlStep;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.equilibrateStep;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.dt;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.controlTimeInterval;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.diffu_t;    
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.diffu_r;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.maxVelocity;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.velocityChangePoint;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.Bpp;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.Os_pressure;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.L_dep;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.cutoff;   
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.kappa;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.seed;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.selfAvoidanceFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.assignmentMethod;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.noControlFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.shapeFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.motionFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.CollectiveMoveCycle >> parameter.collective_MoveStep >> parameter.collective_RestoreStep;
    getline(runfile, line);
    getline(runfile, line); 
    runfile >> parameter.cargoTransFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.cargoCaptureStep;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.trajOutputInterval;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.landmarkLength >> parameter.landmarkMin;;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.landmarkDist;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.blockCost;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.obstacleFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.dynamicTargetFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.targetDiffuseRatio >> parameter.targetVelocityRatio >> parameter.cargoInteractingFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.targetCenter[0] >> parameter.targetCenter[1] >> parameter.targetCenter[2];
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.particleCellListFlag >> parameter.obstacleCellListFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.binaryVelocityFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.setConstantVFlag >> parameter.setConstantV;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.targetHistoryFlag >> parameter.targetHistorySaveInterval >> parameter.targetHistoryLength;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.cellListCutoff >> parameter.cellListDim >> parameter.cellListMaxCount >>parameter.cellListBox_x
            >> parameter.cellListBox_y>> parameter.cellListBox_z;
    getline(runfile, line);
    getline(runfile, line);
    getline(runfile, parameter.iniConfig);
    getline(runfile, line);
    getline(runfile, parameter.filetag);
    getline(runfile, line);
    getline(runfile, parameter.targetConfig);
    getline(runfile, line);
    getline(runfile, parameter.obstacleFilename);

}

