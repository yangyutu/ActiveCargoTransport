#include"simulator.h"
#include"common.h"

extern Parameter parameter;

Simulator::Simulator(std::shared_ptr<Model> model0,std::shared_ptr<Controller> controller0):
                    model(model0),controller(controller0){
    controlFrequency = parameter.controlTimeInterval/model->dt();
    assignmentFrequency = 1;
    nstep_control = parameter.controlStep;
    nstep_equilibrate = parameter.equilibrateStep;
    motionCycle = parameter.CollectiveMoveCycle;
    move_motionStep = parameter.collective_MoveStep;
    move_recoverStep = parameter.collective_RestoreStep;
}

void Simulator::shapeForming(){
    double &totalCost = parameter.totalCost;   
    model->createInitialState();
    totalCost = controller->calAssignment(model->getCurrState(),model->getTargets(),model->getDimP());
    controller->calControl(model->getCurrState(),model->getTargets(),model->getDimP());
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
    double &totalCost = parameter.totalCost;   
    model->createInitialState();
    controller->buildTargetGraph();
    if (parameter.assignmentMethod == 1){
        totalCost = controller->calSeqAssignment(model->getCurrState(),model->getTargets(),0);
    } else  if (parameter.assignmentMethod == 3){
        totalCost = controller->calAssignmentSeqViaShortestPath(model->getCurrState(),model->getTargets(),0);
    }
        controller->calControl(model->getCurrState(),model->getTargets(),model->getDimP());
    int iter;
    std::cout << totalCost << std::endl;
    for(int s=0; s < nstep_control; s++){
        if ((s+1)%assignmentFrequency == 0){
            if (parameter.assignmentMethod == 1){
                totalCost = controller->calSeqAssignment(model->getCurrState(),model->getTargets(),1);
            } else if (parameter.assignmentMethod == 3){
                totalCost = controller->calAssignmentSeqViaShortestPath(model->getCurrState(),model->getTargets(),1);
    
            }
            std::cout << s << "\t totalCost  " <<totalCost <<"\t number of current targets: " << controller->getNumMarked()
                    <<"\t current average deviations:  " <<controller->getDeviation()<< std::endl;
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
    double &totalCost= parameter.totalCost;   
    model->createInitialState();
    totalCost = controller->calAssignment(model->getCurrState(),model->getTargets(),model->getDimP());
    controller->calControl(model->getCurrState(),model->getTargets(),model->getDimP());
    std::cout << totalCost << std::endl;
    for(int c = 0; c < motionCycle; c++){
        for (int s = 0; s < move_motionStep; s++) {
            
            controller->translate_2d(0.0, model->getCurrState());
            model->run(controlFrequency);
            controller->alignTarget_rt(model->getCurrState(),model->getTargets());
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
    double &totalCost= parameter.totalCost;   
    model->createInitialState();
//    totalCost = controller->calAssignment(model->getCurrState(),model->getTargets(),model->getDimP());
//    controller->calControl(model->getCurrState(),model->getTargets(),model->getDimP());
//    int iter;
//    std::cout << totalCost << std::endl;
    
    
    // first do a capture 
    int captureStep = parameter.cargoCaptureStep;
     for (int s = 0; s < captureStep; s++) {
            controller->alignCargo(model->getCurrState(),model->getTargets());
            if ((s + 1) % assignmentFrequency == 0) {
                totalCost = controller->calAssignment(model->getCurrState(), model->getTargets(), model->getDimP());
                std::cout << "capture step: "<< s  << "\t" <<totalCost << std::endl;
            }
            controller->calControl(model->getCurrState(), model->getTargets(), model->getDimP());
            model->getCurrState()[0]->u = 0;
            model->run(controlFrequency);
        }
    
    if (parameter.motionFlag){
    for(int c = 0; c < motionCycle; c++){
        for (int s = 0; s < move_motionStep; s++) {            
//            controller->translateCargo_2d(0.0, model->getCurrState());
            controller->translateCargoFollowPath_2d(model->getCurrState());
            model->run(controlFrequency);
            controller->alignCargo(model->getCurrState(),model->getTargets());
            std::cout << "cargo move step: " << "cycle: " << c << "step: " << s << std::endl;
        }
        
        
        for (int s = 0; s < move_recoverStep; s++) {
            if ((s + 1) % assignmentFrequency == 0) {
                totalCost = controller->calAssignment(model->getCurrState(), model->getTargets(), model->getDimP());
                std::cout << "cargo recover step: " << "cycle: " << c << "step: " << s << "\t" <<totalCost << std::endl;
            }
            controller->calControl(model->getCurrState(), model->getTargets(), model->getDimP());
            model->run(controlFrequency);
            controller->alignCargo(model->getCurrState(),model->getTargets());
        }
    
    }
    }
}

void Simulator::rotate_2d(){
    double &totalCost= parameter.totalCost;   
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
