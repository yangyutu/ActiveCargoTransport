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
    cellList_ptr cell(new CellList(3.0*radius,2,10,150.0*radius,150.0*radius,50.0*radius));
    std::shared_ptr<Model> m(new Model(nullptr));
    std::shared_ptr<Controller> c(new Controller(m->getTargets()));
    Simulator simulator(m,c);
    if (parameter.shapeFlag == 1){
        simulator.shapeForming();
    } else if(parameter.shapeFlag == 2){
        simulator.shapeForming_seq();
    }else if(parameter.motionFlag==1){
        simulator.translate_2d();
    } else if(parameter.motionFlag==2){
        simulator.rotate_2d();
    } else if(parameter.cargoTransFlag==1){
        simulator.cargoTransport_2d();
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
    runfile >> parameter.controlStep;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.equilibrateStep;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.dt;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.diffu_t;    
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.diffu_r;
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
    runfile >> parameter.assignViaEud;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.shapeFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.motionFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.cargoTransFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.trajOutputInterval;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.landmarkFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.landmarkLength;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.landmarkDist;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.blockCost;
    getline(runfile, line);
    getline(runfile, line);
    getline(runfile, parameter.iniConfig);
    getline(runfile, line);
    getline(runfile, parameter.filetag);
    getline(runfile, line);
    getline(runfile, parameter.targetConfig);


}
