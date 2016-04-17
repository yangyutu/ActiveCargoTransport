#include<memory>
#include "model.h"
#include "controller.h"

class Simulator{
public:
    Simulator(){}
    Simulator(std::shared_ptr<Model> model0,std::shared_ptr<Controller> controller0);
    ~Simulator(){}
    void shapeForming();
    void shapeForming_seq();
    void translate_2d();
    void rotate_2d();
    void cargoTransport_2d();
    
private:
    std::shared_ptr<Model> model;
    std::shared_ptr<Controller> controller;
    int nstep_equilibrate, nstep_control, motionCycle, move_motionStep, move_recoverStep;
    int assignmentFrequency;
    int controlFrequency;
};