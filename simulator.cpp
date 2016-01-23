#include"simulator.h"
#include"common.h"

extern Parameter parameter;

Simulator::Simulator(std::shared_ptr<Model> model0,std::shared_ptr<Controller> controller0):
                    model(model0),controller(controller0){
    controlFrequency = 1.0/model->dt();
//    controlFrequency = 1000;
    assignmentFrequency = 1;
    nstep_control = parameter.controlStep;
    nstep_equilibrate = parameter.equilibrateStep;
//    motionCycle = parameter.motionCycle;
//    motionFlag = parameter.motionFlag;
    motionCycle = 70;
    move_motionStep = 2;
    move_recoverStep = 20;
}

void Simulator::shapeForming(){
    double totalCost;   
    model->createInitialState();
    totalCost = controller->calAssignment(model->getCurrState(),model->getTargets(),model->getDimP());
    controller->calControl(model->getCurrState(),model->getTargets(),model->getDimP());
    int iter;
    std::cout << totalCost << std::endl;
    for(int s=0; s < nstep_control; s++){
        if ((s+1)%assignmentFrequency == 0){
            totalCost = controller->calAssignment(model->getCurrState(),model->getTargets(),model->getDimP());
            std::cout << s << "\t" <<totalCost << std::endl;
        }
        controller->calControl(model->getCurrState(),model->getTargets(),model->getDimP());        
        model->run(controlFrequency);
    }
    model->setControl(0);
    for(int s = 0; s < nstep_equilibrate; s++){
        model->run(controlFrequency);
    }
}


void Simulator::shapeForming_seq(){
    double totalCost;   
    model->createInitialState();
    controller->buildTargetGraph();
    totalCost = controller->calSeqAssignment(model->getCurrState(),model->getTargets(),0);
    controller->calControl(model->getCurrState(),model->getTargets(),model->getDimP());
    int iter;
    std::cout << totalCost << std::endl;
    for(int s=0; s < nstep_control; s++){
        if ((s+1)%assignmentFrequency == 0){
            totalCost = controller->calSeqAssignment(model->getCurrState(),model->getTargets(),1);
            std::cout << s << "\t" <<totalCost <<"\t" << controller->getNumMarked()<<"\t" <<controller->getDeviation()<< std::endl;
        }
        controller->calControl(model->getCurrState(),model->getTargets(),model->getDimP());        
        model->run(controlFrequency);
    }
    model->setControl(0);
    for(int s = 0; s < nstep_equilibrate; s++){
        model->run(controlFrequency);
    }
}


void Simulator::translate_2d(){
    double totalCost;   
    model->createInitialState();
    totalCost = controller->calAssignment(model->getCurrState(),model->getTargets(),model->getDimP());
    controller->calControl(model->getCurrState(),model->getTargets(),model->getDimP());
    int iter;
    std::cout << totalCost << std::endl;
    for(int c = 0; c < motionCycle; c++){
        for (int s = 0; s < move_motionStep; s++) {
            
            controller->translate_2d(0.0, model->getCurrState());
            model->run(controlFrequency);
            controller->alignTarget_rt(model->getCurrState(),model->getTargets());
//            controller->register_2d(model->getCurrState(),model->getTargets());
        }
        
        for (int s = 0; s < move_recoverStep; s++) {
            if ((s + 1) % assignmentFrequency == 0) {
                totalCost = controller->calAssignment(model->getCurrState(), model->getTargets(), model->getDimP());
                std::cout << totalCost << std::endl;
            }
            controller->calControl(model->getCurrState(), model->getTargets(), model->getDimP());
            model->run(controlFrequency);
            controller->alignTarget_rt(model->getCurrState(),model->getTargets());
        }
    
    }
}


void Simulator::cargoTransport_2d(){
    double totalCost;   
    model->createInitialState();
    totalCost = controller->calAssignment(model->getCurrState(),model->getTargets(),model->getDimP());
    controller->calControl(model->getCurrState(),model->getTargets(),model->getDimP());
    int iter;
    std::cout << totalCost << std::endl;
    
    
    // first do a capture 
    int captureStep = 100;
     for (int s = 0; s < captureStep; s++) {
            controller->alignCargo(model->getCurrState(),model->getTargets());
            if ((s + 1) % assignmentFrequency == 0) {
                totalCost = controller->calAssignment(model->getCurrState(), model->getTargets(), model->getDimP());
                std::cout << totalCost << std::endl;
            }
            controller->calControl(model->getCurrState(), model->getTargets(), model->getDimP());
            model->run(controlFrequency);

        }
    
    
    for(int c = 0; c < motionCycle; c++){



        for (int s = 0; s < move_motionStep; s++) {
            
            controller->translateCargo_2d(0.0, model->getCurrState());
            model->run(controlFrequency);
            controller->alignCargo(model->getCurrState(),model->getTargets());
//            controller->register_2d(model->getCurrState(),model->getTargets());
        }
        
        
        for (int s = 0; s < move_recoverStep; s++) {
            if ((s + 1) % assignmentFrequency == 0) {
                totalCost = controller->calAssignment(model->getCurrState(), model->getTargets(), model->getDimP());
                std::cout << totalCost << std::endl;
            }
            controller->calControl(model->getCurrState(), model->getTargets(), model->getDimP());
            model->run(controlFrequency);
            controller->alignCargo(model->getCurrState(),model->getTargets());
        }
    
    }
}

void Simulator::rotate_2d(){
    double totalCost;   
    model->createInitialState();
    totalCost = controller->calAssignment(model->getCurrState(),model->getTargets(),model->getDimP());
    controller->calControl(model->getCurrState(),model->getTargets(),model->getDimP());
    int iter;
    std::cout << totalCost << std::endl;
    for(int c = 0; c < motionCycle; c++){
        for (int s = 0; s < move_motionStep; s++) {
            
            controller->rotate_2d( model->getCurrState());
            model->run(controlFrequency);
            controller->alignTarget_rt(model->getCurrState(),model->getTargets());
//            controller->register_2d(model->getCurrState(),model->getTargets());
        }
        
        for (int s = 0; s < move_recoverStep; s++) {
            if ((s + 1) % assignmentFrequency == 0) {
                totalCost = controller->calAssignment(model->getCurrState(), model->getTargets(), model->getDimP());
                std::cout << totalCost << std::endl;
            }
            controller->calControl(model->getCurrState(), model->getTargets(), model->getDimP());
            model->run(controlFrequency);
            controller->alignTarget_rt(model->getCurrState(),model->getTargets());
//            controller->register_2d(model->getCurrState(),model->getTargets());
        }
    
    }
}
