#include "model.h"
#include "CellList.h"
#include "common.h"
#include<algorithm>
double const Model::T = 293.0;
double const Model::kb = 1.38e-23;
double const Model::vis = 1e-3;

extern Parameter parameter;

Model::Model(){
    rand_normal = std::make_shared<std::normal_distribution<double>>(0.0, 1.0);

   
    filetag = parameter.filetag;
    iniFile = parameter.iniConfig;
    numP = parameter.N;
    dimP = parameter.dim;
    radius = parameter.radius;
    dt_ = parameter.dt;
    diffusivity_t = parameter.diffu_t;// this corresponds the diffusivity of 1um particle
    diffusivity_r = parameter.diffu_r; // this correponds to rotation diffusity of 1um particle
    Bpp = parameter.Bpp * kb * T * 1e9; //2.29 is Bpp/a/kT
    Kappa = parameter.kappa; // here is kappa*radius
    Os_pressure = parameter.Os_pressure * kb * T * 1e9;
    L_dep = parameter.L_dep; // 0.2 of radius size, i.e. 200 nm
    radius_nm = radius*1e9;
    combinedSize = (1+L_dep)*radius_nm;
    mobility = diffusivity_t/kb/T;
    trajOutputInterval = parameter.trajOutputInterval;
    fileCounter = 0;
    cutoff = parameter.cutoff;
    numControl = this->velocity.size();
    this->rand_generator.seed(parameter.seed);

    for(int i = 0; i < numP; i++){
        particles.push_back(particle_ptr(new Model::particle));
        targets.push_back(particle_ptr(new Model::particle));
    }
    
    
    this->getPermutator();
    
    cellListFlag = false;
    if (parameter.particleCellListFlag){
        this->cellList = std::make_shared<CellList>(parameter.cellListCutoff*radius,
                parameter.cellListDim,parameter.cellListMaxCount,
                parameter.cellListBox_x*radius, parameter.cellListBox_y*radius,
                parameter.cellListBox_z*radius);
        cellListFlag = true;
    }
    if (parameter.obstacleCellListFlag){
        this->obsCellList = std::make_shared<CellList>(parameter.cellListCutoff*radius,
                parameter.cellListDim,parameter.cellListMaxCount,
                parameter.cellListBox_x*radius, parameter.cellListBox_y*radius,
                parameter.cellListBox_z*radius);
    }
    
    if (parameter.obstacleFlag) {
        this->readObstacle();
        this->obsCellList = std::make_shared<CellList>(parameter.cellListCutoff*radius,
                parameter.cellListDim,parameter.cellListMaxCount,
                parameter.cellListBox_x*radius, parameter.cellListBox_y*radius,
                parameter.cellListBox_z*radius);
        int builtCount = this->obsCellList->buildList(obstacles);
            
            if (builtCount!=numObstacles){
                std::cout << "build imcomplete" << std::endl;
            }

    }
    
}

void Model::run() {
    if (this->timeCounter == 0 || ((this->timeCounter + 1) % trajOutputInterval == 0)) {
        // here we need to calculate the initial energy
        for (int i = 0; i < numP; i++) {
            for (int k = 0; k < dimP; k++) {
               particles[i]->potential = 0.0;
            }
        }
    

        for (int i = 0; i < numP - 1; i++) {
            for (int j = i + 1; j < numP; j++) {
                double pot;
                double F[3];
                calForcesHelper_DLAO(particles[i]->r, particles[j]->r, F,i, j,pot);
                particles[i]->potential += 0.5*pot;
                particles[j]->potential += 0.5*pot;
            }
        }
        
        
        this->outputTrajectory(this->trajOs);
        this->outputOrderParameter(this->opOs);
    }
    

    
    if (cellListFlag){
            int builtCount = cellList->buildList(particles);
            if (builtCount!=numP){
                std::cout << "build imcomplete" << std::endl;
            }
    }
    calForces();

    
    if (dimP == 2){
        for (int i = 0; i < numP; i++) {
            
            double random_x = sqrt(2.0 * diffusivity_t / dt_) * (*rand_normal)(rand_generator);

        //    random_x = 0.0;
            
            double random_y = sqrt(2.0 * diffusivity_t / dt_) * (*rand_normal)(rand_generator);
            
        //    random_y = 0.0;
            double ux = mobility * particles[i]->F[0]  +
                        parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi) 
                    +   random_x;
            
            
            
            
            double uy = mobility * particles[i]->F[1]  +
                        parameter.maxVelocity*particles[i]->u * sin(particles[i]->phi) 
                    +   random_y;
            
            
            particles[i]->r[0] += ux * dt_;
            particles[i]->r[1] += uy * dt_;
        
            particles[i]->energy_accumlator1 += (parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi)*ux + 
                        parameter.maxVelocity*particles[i]->u * sin(particles[i]->phi)*uy)*dt_;
            

            
            particles[i]->energy_accumlator2 += (ux*ux + uy*uy)*dt_;
            
            
            particles[i]->energy_accumlator3 += (ux* mobility * particles[i]->F[0]+
                    uy* mobility * particles[i]->F[1] )*dt_;
                        
            
            particles[i]->energy_accumlator4 += (ux* random_x+
                    uy* random_y)*dt_;
            
            if (!particles[i]->transporterFlag){// counting maintaining input energy
                particles[i]->eneregy_maintain_accumlator += (parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi)*ux + 
                        parameter.maxVelocity*particles[i]->u * sin(particles[i]->phi)*uy)*dt_;
            }
            particles[i]->Fx = mobility * particles[i]->F[0];
            particles[i]->Fy = mobility * particles[i]->F[1];
            
            particles[i]->Vx = ux;
            particles[i]->Vy = uy;
            
            particles[i]->phi += sqrt(2.0 * diffusivity_r * dt_) * (*rand_normal)(rand_generator);
            
            
            if( particles[i]->Vx > 0){
                particles[i]->instant_output_work = ux*(mobility * particles[i]->F[0]+parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi))*dt_;
            }
            
            particles[i]->instant_conserve_work = (ux* mobility * particles[i]->F[0]+
                    uy* mobility * particles[i]->F[1] )*dt_;
            particles[i]->instant_input_work = (parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi)*ux + 
                        parameter.maxVelocity*particles[i]->u * sin(particles[i]->phi)*uy)*dt_;
            
            particles[i]->friction_accumulator1 += (particles[i]->Fx + parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi))*dt_;
            particles[i]->friction_accumulator2 += ux*dt_;
            particles[i]->friction_accumulator3 += (parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi))*dt_;
            particles[i]->friction_accumulator4 += (random_x)*dt_;
           
            
            double motor_target_vector[2];
            
            motor_target_vector[0] = particles[i]->targetPos[0]*radius - particles[i]->r[0];
            motor_target_vector[1] = particles[i]->targetPos[1]*radius - particles[i]->r[1];
            double vector_norm;
            vector_norm = pow(motor_target_vector[0],2) + pow(motor_target_vector[1],2);
            vector_norm = sqrt(vector_norm);
            motor_target_vector[0] = motor_target_vector[0]/vector_norm;
            motor_target_vector[1] = motor_target_vector[1]/vector_norm;
            
            
            double v_projection = ux*motor_target_vector[0]+uy*motor_target_vector[1];
            particles[i]->v_projection = 0.0;
            double v_sp_projection = parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi)*motor_target_vector[0] + 
            parameter.maxVelocity*particles[i]->u * sin(particles[i]->phi)*motor_target_vector[1];
            
            double vsp_x,vsp_y;
            vsp_x = parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi);
            vsp_y = parameter.maxVelocity*particles[i]->u * sin(particles[i]->phi);
            
            
            double ratio = cos(particles[i]->phi)*motor_target_vector[0] + sin(particles[i]->phi)*motor_target_vector[1];
//            particles[i]->energy_vsp_input_maintain += (vsp_x*vsp_x + vsp_y*vsp_y)*dt_;
            
            particles[i]->energy_instant_vsp_input = vsp_x*vsp_x + vsp_y*vsp_y;
            
            if (!particles[i]->transporterFlag){
//                particles[i]->eneregy_useful_maintain_accumlator += (parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi)*ux + 
//                        parameter.maxVelocity*particles[i]->u * sin(particles[i]->phi)*uy)*dt_*ratio*ratio;
                particles[i]->eneregy_useful_maintain_accumlator += (v_projection*v_sp_projection)*dt_;
                particles[i]->v_projection = ratio;
                particles[i]->maintainratio = ratio*v_projection/sqrt(ux*ux + uy*uy);
     //           if (ratio < 0){ // here is to prevent overshoot
     //               particles[i]->u = 0.0;
     //           }
                particles[i]->energy_vsp_input_maintain += (vsp_x*vsp_x + vsp_y*vsp_y)*dt_;
                particles[i]->energy_vsp_useful_maintain += (v_sp_projection*v_sp_projection)*dt_;
                
                particles[i]->energy_instant_vsp_maintain = (v_sp_projection*v_sp_projection);
                particles[i]->energy_instant_vsp_transport = 0.0;
                
                
                
            } else{
                particles[i]->energy_input_transport += (parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi)*ux + 
                        parameter.maxVelocity*particles[i]->u * sin(particles[i]->phi)*uy)*dt_;
                particles[i]->energy_useful_transport += (parameter.maxVelocity*particles[i]->u * cos(particles[i]->phi)*ux)*dt_;
            
                particles[i]->energy_vsp_input_transport += (vsp_x*vsp_x + vsp_y*vsp_y)*dt_;
                particles[i]->energy_vsp_useful_transport += (vsp_x*vsp_x)*dt_;
                
                particles[i]->energy_instant_vsp_maintain = 0;
                particles[i]->energy_instant_vsp_transport = vsp_x*vsp_x;
            
            }
        
        
        
        }
        
        
        if (parameter.dynamicTargetFlag){
            

            
            
            double randDist[2], disp[2];
            randDist[0] = (*rand_normal)(rand_generator);
            randDist[1] = (*rand_normal)(rand_generator);
            
    //        randDist[0] = 0.0;
    //        randDist[1] = 0.0;
            
            
            disp[0] = parameter.targetDiffuseRatio*mobility * targetCenter.F[0] * dt_ 
                    + (parameter.targetVelocityRatio*parameter.maxVelocity*dt_ + 
                    sqrt(2.0 * diffusivity_t*parameter.targetDiffuseRatio * dt_) * randDist[0]);
            disp[1] = parameter.targetDiffuseRatio*mobility * targetCenter.F[1] * dt_ 
                    +sqrt(2.0 * diffusivity_t*parameter.targetDiffuseRatio * dt_) * randDist[1];
            
            
            // disp is in the unit of m
            targetCenter.r[0] += disp[0];
            targetCenter.r[1] += disp[1];
            
            
                  
            

            
            if (parameter.targetHistoryFlag == 1){// whether to use target history based tracking, use time basis criterion
                
                // first we need to update our target history buffer
                if ( (this->timeCounter + 1) % parameter.targetHistorySaveInterval == 0){
                    int idx = this->targetCenter_historyCounter++%parameter.targetHistoryLength;
                    this->targetCenter_history(idx,0) = targetCenter.r[0];
                    this->targetCenter_history(idx,1) = targetCenter.r[1];
                }
                
                arma::mat avg = arma::mean(targetCenter_history);// here calculate the average target position
                this->targetCenter_avg.r[0] = avg(0,0);
                this->targetCenter_avg.r[1] = avg(0,1);
                this->targetCenter_avg.r[2] = avg(0,2);
                for (int i = 0; i < numP; i++) {
                    targets[i]->r[0] = (parameter.targetCenter[0] + initialDistToCenter[i]->r[0]) + avg(0,1)/radius;
                    targets[i]->r[1] = (parameter.targetCenter[1] + initialDistToCenter[i]->r[1]) + avg(0,2)/radius;
                }
            }else if (parameter.targetHistoryFlag == 2){ // use distance based criterion to update target position
                double dispX, dispY;
                dispX = this->previousTargetCenter.r[0] - targetCenter.r[0];
                dispY = this->previousTargetCenter.r[1] - targetCenter.r[1];        
                accumTargetMove = sqrt(pow(dispX/radius,2) + pow(dispY/radius,2));
                if (accumTargetMove > parameter.targetMoveThresh){
                    
                    for (int i = 0; i < numP; i++) {
                        targets[i]->r[0] = (parameter.targetCenter[0] + initialDistToCenter[i]->r[0]) + targetCenter.r[0]/radius;
                        targets[i]->r[1] = (parameter.targetCenter[1] + initialDistToCenter[i]->r[1]) + targetCenter.r[1]/radius;
                    }
                    this->previousTargetCenter.r[0] = targetCenter.r[0];
                    this->previousTargetCenter.r[1] = targetCenter.r[1];
                    
                }
            
            
            }else{
            
                for (int i = 0; i < numP; i++) {
                    targets[i]->r[0] += disp[0] / radius;
                    targets[i]->r[1] += disp[1] / radius;
                }
            }
    
        }   
    } else if(dimP == 3){
        std::cerr << "3d equation of motion is depreciated in the moment!" << std::endl;
        exit(3);
        for (int i = 0; i < numP; i++) {
            for(int j = 0; j < dimP; j++){
                particles[i]->r[j] += mobility * particles[i]->F[j] * dt_ +
                        velocity[particles[i]->u] * particles[i]->ori_vec[0][j]* dt_
                    +   sqrt(2.0 * diffusivity_t * dt_) * (*rand_normal)(rand_generator);
            }
        // update the orientation vector via spherical surface diffusion
            double phitemp = sqrt(2.0 * diffusivity_r * dt_) * (*rand_normal)(rand_generator);
            double thetatemp = sqrt(2.0 * diffusivity_r * dt_) * (*rand_normal)(rand_generator);
            for(int j = 0; j < dimP; j++){
                particles[i]->ori_vec[0][j] += particles[i]->ori_vec[1][j]*phitemp + 
                    particles[i]->ori_vec[2][j]*thetatemp;
            }
        }
        this->updateBodyFrameVec();  
        
    }
    this->timeCounter++;
    
}

void Model::run(int steps){
    for (int i = 0; i < steps; i++){
	run();
    }
}

void Model::setControl(int c){
    for (int i = 0; i < numP; i++) {
        particles[i]->u = c;
    }

}
void Model::readTarget(std::string filename){
    std::ifstream is;
    is.open(filename);
    std::string line;
    std::stringstream linestream;
    double dum;
    for (int i = 0; i < numP; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> targets[i]->r[0];
        linestream >> targets[i]->r[1];
        linestream >> targets[i]->r[2];
        linestream >> dum;
        linestream >> dum;
        linestream >> targets[i]->marked;
        // now do the target shift respect to the target center
        targets[i]->r[0] += parameter.targetCenter[0];
        targets[i]->r[1] += parameter.targetCenter[1];
        targets[i]->r[2] += parameter.targetCenter[2];
        
    }
    
 }

// this force calculation includes double layer repulsion and depletion attraction 
void Model::calForcesHelper_DLAO(double ri[3], double rj[3], double F[3],int i,int j,double& pot) {
    double r[dimP], dist;

    dist = 0.0;
    for (int k = 0; k < dimP; k++) {
        F[k] = 0.0;
        r[k] = (rj[k] - ri[k]) / radius;
        dist += pow(r[k], 2.0);
    }
    dist = sqrt(dist);
    if (dist < 2.0) {
        std::cerr << "overlap " << i << "\t" << j << "\t"<< this->timeCounter << "dist: " << dist <<std::endl;
        dist = 2.06;
    }
    if (dist < cutoff) {
        double Fpp = -4.0/3.0*
        Os_pressure*M_PI*(-3.0/4.0*pow(combinedSize,2.0)+3.0*dist*dist/16.0*radius_nm*radius_nm);
        Fpp += -Bpp * Kappa * exp(-Kappa*(dist-2.0));
        for (int k = 0; k < dimP; k++) {
            F[k] = Fpp * r[k] / dist;

        }
        //    Bpp = parameter.Bpp * kb * T * 1e9; //2.29 is Bpp/a/kT, a in unit of nm
        pot = parameter.Bpp * radius_nm * exp(-Kappa*(dist-2.0));
        double norm_comb = 1 +  L_dep;
        double ao = -4.0/3.0*parameter.Os_pressure*pow(radius_nm,3) *M_PI*(pow(norm_comb,3) -3.0/4.0*pow(norm_comb,2) *dist+1.0/16.0* pow(dist,3));
        pot += ao;
    }
}

// this force calculation only includes double layer repulsion 
void Model::calForcesHelper_DL(double ri[3], double rj[3], double F[3],int i, int j) {
    double r[dimP], dist;

    dist = 0.0;
    for (int k = 0; k < dimP; k++) {
        F[k] = 0.0;
        r[k] = (rj[k] - ri[k]) / radius;
        dist += pow(r[k], 2.0);
    }
    dist = sqrt(dist);
    if (dist < 2.0) {
        std::cerr << "overlap " << i << "\t with " << j << "\t"<< this->timeCounter << "dist: " << dist <<std::endl;
        dist = 2.06;
    }
    if (dist < cutoff) {
        double Fpp = -Bpp * Kappa * exp(-Kappa*(dist-1.9));
        
        for (int k = 0; k < dimP; k++) {
            F[k] = Fpp * r[k] / dist;
        }
    }
}

void Model::calForces() {
    double r[dimP], dist, F[3];
    for (int i = 0; i < numP; i++) {
        for (int k = 0; k < dimP; k++) {
            particles[i]->F[k] = 0.0;
            particles[i]->potential = 0.0;
        }
    }
    
    if(!cellListFlag){

        for (int i = 0; i < numP - 1; i++) {
            for (int j = i + 1; j < numP; j++) {
                double pot;
                calForcesHelper_DLAO(particles[i]->r, particles[j]->r, F,i, j,pot);
                particles[i]->potential += 0.5*pot;
                particles[j]->potential += 0.5*pot;
                
                for (int k = 0; k < dimP; k++) {
                    particles[i]->F[k] += F[k];
                    particles[j]->F[k] += -F[k];
                }
            }
        }
        
        if (parameter.cargoInteractingFlag){
            // Note here we treat the target center as the cargo
            for (int k = 0; k < dimP; k++) {
                targetCenter.F[k] = 0.0;
            }
            
            for (int j = 0; j < numP; j++) {
                calForcesHelper_DL(targetCenter.r, particles[j]->r, F,-1, j);
                for (int k = 0; k < dimP; k++) {
                    targetCenter.F[k] += F[k];
                    particles[j]->F[k] += -F[k];
                }
            }
        
        }
        
    } else{
        
        std::cerr << "this force calculation using cellList is depreciated in this repo! " << std::endl;
        exit(6);
        for (int i = 0; i < numP; i++) {
            std::vector<int> nblist = 
            cellList->getNeighbors(particles[i]->r[0],particles[i]->r[1],particles[i]->r[2]);
            int nblistSize = nblist.size();
            for (int j = 0; j < nblistSize; j++){
                if (i!=nblist[j]){
                    double pot;
                    calForcesHelper_DLAO(particles[i]->r, particles[nblist[j]]->r, F,i, nblist[j],pot);
                    for (int k = 0; k < dimP; k++) {
                        particles[i]->F[k] += F[k];
                    }
                }
            } 
        }    
    }
    
    if(parameter.obstacleFlag) {
        for (int i = 0; i < numP; i++) {
            std::vector<int> nblist =
                    obsCellList->getNeighbors(particles[i]->r[0], particles[i]->r[1], particles[i]->r[2]);
            int nblistSize = nblist.size();
            for (int j = 0; j < nblistSize; j++) {
                calForcesHelper_DL(particles[i]->r, obstacles[nblist[j]]->r, F, i, nblist[j]);
                for (int k = 0; k < dimP; k++) {
                    particles[i]->F[k] += F[k];
                }
            }

        }
        
        // Note here we treat the target center as the cargo
        std::vector<int> nblist =
                    obsCellList->getNeighbors(targetCenter.r[0], targetCenter.r[1], targetCenter.r[2]);
        int nblistSize = nblist.size();
        for (int j = 0; j < nblistSize; j++) {
            calForcesHelper_DL(targetCenter.r, obstacles[nblist[j]]->r, F, -1, nblist[j]);
            for (int k = 0; k < dimP; k++) {
                targetCenter.F[k] += F[k];
            }
        }
    
    }
    
}
    


void Model::createInitialState(){

    this->readxyz(iniFile);
    this->readTarget(parameter.targetConfig);
    this->targetCenter.r[0] = parameter.targetCenter[0]*radius;
    this->targetCenter.r[1] = parameter.targetCenter[1]*radius;
    this->targetCenter.r[2] = parameter.targetCenter[2]*radius;
    this->previousTargetCenter.r[0] = parameter.targetCenter[0]*radius;
    this->previousTargetCenter.r[1] = parameter.targetCenter[1]*radius;
    this->previousTargetCenter.r[2] = parameter.targetCenter[2]*radius;    
    
    
    for (int i = 0; i < numP; i++){
        initialDistToCenter.push_back(particle_ptr(new Model::particle));
        initialDistToCenter[i]->r[0] = this->targets[i]->r[0] - parameter.targetCenter[0];
        initialDistToCenter[i]->r[1] = this->targets[i]->r[1] - parameter.targetCenter[1];
        initialDistToCenter[i]->r[2] = this->targets[i]->r[2] - parameter.targetCenter[2];
    }
    
//  initialize target history 
    this->targetCenter_historyCounter = 0;
    this->targetCenter_history.set_size(parameter.targetHistoryLength,3);
    for (int i = 0; i < parameter.targetHistoryLength; i++){
        for (int j = 0; j < 3; j++){
            this->targetCenter_history(i,j) = this->targetCenter.r[j];
        
        }
    }
    // initialize target move tracker
    accumTargetMove = 0.0;
    
    
    
    
    std::stringstream ss;
    std::cout << "model initialize at round " << fileCounter << std::endl;
    ss << this->fileCounter++;
    if (trajOs.is_open()) trajOs.close();
    if (opOs.is_open()) opOs.close();
    if (osCargo.is_open()) osCargo.close();
    this->trajOs.open(filetag + "xyz_" + ss.str() + ".txt");
    this->opOs.open(filetag + "op" + ss.str() + ".txt");
    this->osTarget.open(filetag +"target"+ss.str() + ".txt");
    this->osCargo.open(filetag +"cargo"+ss.str() + ".txt");
    this->timeCounter = 0;

}

void Model::outputTrajectory(std::ostream& os) {

	/*
	columnIndex    dataName
	0				motor ID
	1               x
	2               y
	3				z
	4				phi
	5				theta
	6				cost
	7				u (from 0 to 1)
	8				targetIdx
	9				ShortestPathDistToTarget
	10				availControl
	11				t_x
	12				t_y
	13				t_z
	14				time
	15				Transporter Flag
	16				energy_accum1
	17				energy_accum2
	18				energy_accum3
	19				energy_accum4
	20				Fx
	21				Fy
	22				Vx
	23				Vy
	24				potential
	25				energy_maintain_accum
	26				instant_input_work
	27				instant_conserve_work
	28				instant_output_work 
	29				friction_accumulator1
	30				friction_accumulator2
	31				friction_accumulator3
	32				friction_accumulator4
	33				eneregy_useful_maintain_accumlator
	34				v_projection
	35				maintainratio
	36				energy_input_transport
	37				energy_useful_transport
	38				energy_vsp_input_maintain accumulated
	39				energy_vsp_useful_maintain accumulated 
	40				energy_vsp_input_transport accumulated
	41				energy_vsp_useful_transport accumulated
	42				energy_instant_vsp_input   this is instant input energy
	43				energy_instant_vsp_transport
	44				energy_instant_vsp_maintain

	
	
	
	*/


    for (int i = 0; i < numP; i++) {
        os << i << "\t";
        this->osTarget << i << "\t";
        for (int j = 0; j < 3; j++){
            os << particles[i]->r[j]/radius << "\t";
        }
        // first four columns id, x, y z
        for (int j = 0; j < 3; j++){
            osTarget << targets[i]->r[j] << "\t";
        }
        osTarget << targets[i]->marked << "\t";
        osTarget << std::endl;
        
        os << particles[i]->phi<< "\t";
        os << particles[i]->theta<< "\t";
        os << particles[i]->cost<< "\t";
        os << particles[i]->u<< "\t";
        os << particles[i]->targetIdx<< "\t";
        os << particles[i]->ShortestPathDistToTarget<< "\t";
        os << particles[i]->availControl<<"\t";
        os << targets[i]->r[0]<< "\t";
        os << targets[i]->r[1]<< "\t";
        os << targets[i]->r[2]<<"\t";
        os << this->timeCounter*this->dt_ << "\t";
        os << particles[i]->transporterFlag<<"\t";
        os << particles[i]->energy_accumlator1<<"\t";
        os << particles[i]->energy_accumlator2<<"\t";
        os << particles[i]->energy_accumlator3<<"\t";
        os << particles[i]->energy_accumlator4<<"\t";
        os << particles[i]->Fx<<"\t";
        os << particles[i]->Fy<<"\t";
        os << particles[i]->Vx<<"\t";
        os << particles[i]->Vy<<"\t";
        os << particles[i]->potential<<"\t";
        os << particles[i]->eneregy_maintain_accumlator<<"\t";
        os << particles[i]->instant_input_work << "\t";
        os << particles[i]->instant_conserve_work << "\t";
        os << particles[i]->instant_output_work << "\t";
        os << particles[i]->friction_accumulator1 << "\t";
        os << particles[i]->friction_accumulator2 << "\t";
        os << particles[i]->friction_accumulator3 << "\t";
        os << particles[i]->friction_accumulator4 << "\t";       
        os << particles[i]->eneregy_useful_maintain_accumlator << "\t";
        os << particles[i]->v_projection << "\t";
        os << particles[i]->maintainratio << "\t";
        os << particles[i]->energy_input_transport << "\t";
        os << particles[i]->energy_useful_transport << "\t";
        os << particles[i]->energy_vsp_input_maintain << "\t";
        os << particles[i]->energy_vsp_useful_maintain << "\t";
 //       os << particles[i]->energy_vsp_useful_maintain_positive << "\t";
        os << particles[i]->energy_vsp_input_transport << "\t";
        os << particles[i]->energy_vsp_useful_transport << "\t";
        os << particles[i]->energy_instant_vsp_input << "\t";
        os << particles[i]->energy_instant_vsp_transport << "\t";
        os << particles[i]->energy_instant_vsp_maintain << "\t";
        os << std::endl;
    }
    for (int j = 0; j < 3; j++){
        this->osCargo << targetCenter.r[j]/radius << "\t";
    }
    for (int j = 0; j < 3; j++){
        this->osCargo << targetCenter_avg.r[j]/radius << "\t";
    }
    this->osCargo << this->timeCounter*this->dt_ << "\t";
    this->osCargo << std::endl;
    
}

void Model::outputOrderParameter(std::ostream& os) {

    os << this->timeCounter << "\t";
//    os << this->calHausdorff() << "\t";
//    os << this->calPsi6() << "\t";
    os << this->calRg() << "\t";
    os << this->calEudDeviation() << "\t";
    os << parameter.totalCost << "\t";
    os << parameter.totalCost/parameter.N << std::endl;
}


void Model::readxyz(const std::string filename) {
    std::ifstream is;
    is.open(filename.c_str());
    std::string line;
    double dum;
    for (int i = 0; i < numP; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> particles[i]->r[0];
        linestream >> particles[i]->r[1];
        linestream >> particles[i]->r[2];
        linestream >> particles[i]->phi;
        linestream >> particles[i]->theta;
    }
    for (int i = 0; i < numP; i++) {
        particles[i]->r[0] *=radius;
        particles[i]->r[1] *=radius;
        particles[i]->r[2] *=radius;
        particles[i]->ori_vec[0][0] = cos(particles[i]->phi)*sin(particles[i]->theta);
        particles[i]->ori_vec[0][1] = sin(particles[i]->phi)*sin(particles[i]->theta);
        particles[i]->ori_vec[0][2] = cos(particles[i]->theta);
        
        // initialize energy input accumulator
        particles[i]->energy_accumlator1 = 0.0;
        particles[i]->energy_accumlator2 = 0.0;
        particles[i]->energy_accumlator3 = 0.0;
        particles[i]->energy_accumlator4 = 0.0;
        particles[i]->eneregy_maintain_accumlator = 0.0;
        particles[i]->instant_input_work = 0.0;
        particles[i]->instant_conserve_work = 0.0;
        particles[i]->instant_output_work =0.0;
        particles[i]->friction_accumulator1 = 0.0;
        particles[i]->friction_accumulator2 =0.0;
        particles[i]->friction_accumulator3 = 0.0;
        particles[i]->friction_accumulator4 =0.0;
        particles[i]->eneregy_useful_maintain_accumlator = 0.0;
        particles[i]->v_projection = 0.0;
        particles[i]->Fx = 0.0;
        particles[i]->Fy = 0.0;
        particles[i]->maintainratio = 0.0;
        particles[i]->energy_input_transport = 0.0;
        particles[i]->energy_useful_transport = 0.0;
        
        particles[i]->energy_vsp_input_maintain =0.0;
        particles[i]->energy_vsp_useful_maintain =0.0;
        particles[i]->energy_vsp_useful_maintain_positive =0.0;
        particles[i]->energy_vsp_input_transport = 0.0;
        particles[i]->energy_vsp_useful_transport = 0.0;
    }
    this->updateBodyFrameVec();
    
    if (cellListFlag){
        cellList->buildList(particles);
    }
    
    is.close();
}

void Model::updateBodyFrameVec(){
    
    double thresh = 0.99999;
    double norm;
    
      for (int ii = 0; ii < numP; ii++){
          // set n2, n3 to zero
          for (int i = 0; i < 3; i++){
              particles[ii]->ori_vec[1][i] = 0.0;
              particles[ii]->ori_vec[2][i] = 0.0;
              
          }
          
          norm = 0.0;
          for(int i = 0; i < 3; i++){
              norm = norm + pow(particles[ii]->ori_vec[0][i],2.0);
          }
                    norm = sqrt(norm);
          for(int i = 0; i < 3; i++){
              particles[ii]->ori_vec[0][i] /= norm;
          }          
          
          
          // first do some safe-guard to prevent numerical instability
          if(particles[ii]->ori_vec[0][2] >= thresh){
              particles[ii]->ori_vec[0][0] = 0.0;
              particles[ii]->ori_vec[0][1] = 0.0;
              particles[ii]->ori_vec[0][2] = 1.0;
          }
      
          if(particles[ii]->ori_vec[0][2] <= -thresh){
              particles[ii]->ori_vec[0][0] = 0.0;
              particles[ii]->ori_vec[0][1] = 0.0;
              particles[ii]->ori_vec[0][2] = -1.0;
          }
      

           
          // first consider degenrate case that n1 = e_z
          if(abs(particles[ii]->ori_vec[0][2]) >= thresh){
             particles[ii]->ori_vec[1][0]=0;
             particles[ii]->ori_vec[1][1]=1.0;
             particles[ii]->ori_vec[1][2]=0;
          }
          else {
              for (int i = 0; i < 3; i++){
//                  for (int j = 0; j < 3; j++){
                      for (int k = 0; k < 3; k++){
                          particles[ii]->ori_vec[1][i] += eps[i][2][k]*particles[ii]->ori_vec[0][k];
                      
                      
                      }
//                  }
            }
            norm = 0.0;
            for (int i = 0; i < 3; i++) {
                norm = norm + pow(particles[ii]->ori_vec[1][i],2.0);
            }
            norm = sqrt(norm);
            for (int i = 0; i < 3; i++) {
                particles[ii]->ori_vec[1][i] /= norm;
            }
          }
          for (int i = 0; i < 3; i++){
                  for (int j = 0; j < 3; j++){
                      for (int k = 0; k < 3; k++){
                          particles[ii]->ori_vec[2][i] += 
                                  eps[i][j][k]*particles[ii]->ori_vec[1][j]
                                  *particles[ii]->ori_vec[0][k];
                      
                      
                      }
                  }
             }
                
    //                for(int i = 0; i < 3; i++){
    //                    for(int j = 0; j <3;j++){
    //                        std::cout << particles[ii]->ori_vec[i][j] << std::endl;
                        
    //                    }
    //                }
     }
}

double Model::calHausdorff(){
    double mindist, maxdist;
    maxdist = 0;
    for(int i = 0; i < numP; i++){
        mindist = 100000;
        for (int j = 0; j < numP; j++){
            double dist = 0;
            for (int k = 0 ; k < 3; k++){
                dist += pow(targets[j]->r[k] - particles[i]->r[k],2.0);                       
            }
            dist = sqrt(dist);
            if (dist < mindist){
                mindist = dist;
            }
        }
        if (mindist > maxdist){
            maxdist = mindist;        
        }
    }
    return maxdist;
}

double Model::calRg(){
    //      calculate Rg
    double xmean = 0;
    double ymean = 0;
    double zmean = 0;

    for (int i = 0; i < numP; i++) {
        xmean = xmean + particles[i]->r[0];
        ymean = ymean + particles[i]->r[1];
        zmean = zmean + particles[i]->r[2];
    }
    xmean /= numP;
    ymean /= numP;
    zmean /= numP;
    double rgmean = 0.0;
    for (int i = 0; i < numP; i++) {
        rgmean = rgmean + (particles[i]->r[0] - xmean)*(particles[i]->r[0] - xmean);
        rgmean = rgmean + (particles[i]->r[1] - ymean)*(particles[i]->r[1] - ymean);
        rgmean = rgmean + (particles[i]->r[2] - zmean)*(particles[i]->r[2] - zmean);
    }
    rgmean /= numP;

    rgmean = sqrt(rgmean)/radius;
    
    return rgmean;
}

double Model::calEudDeviation(){
    //      calculate Rg
    double dev = 0;
    for (int i = 0; i < numP; i++){
        dev += this->particles[i]->ShortestPathDistToTarget;
    }
    dev /= numP;
    return dev;
}


double Model::calPsi6(){
    
    double rmin=2.7;
    std::vector<double> psir(numP,0.0),psii(numP,0.0);
    for (int i = 0; i < numP; i++) {
        int nb = 0;
        for (int j = 0; j < numP; j++) {
            if (i != j) {
                double rxij = particles[j]->r[0] - particles[i]->r[0];
                double ryij = particles[j]->r[1] - particles[i]->r[1];
                double RP = sqrt(rxij * rxij + ryij * ryij)/radius;
                if (RP < rmin) {
                    nb += 1;
                    double theta = std::atan2(ryij, rxij);
                    psir[i] += cos(6 * theta);
                    psii[i] += sin(6 * theta);
                }
            }
                      
        } 
            if (nb > 0) {
                psir[i] /=  nb;
                psii[i] /=  nb;
            }
    }
    double psi6 = 0;
    double accumpsi6r = 0;
    double accumpsi6i = 0;
    for (int i = 0; i < numP; i++) {

        accumpsi6r = accumpsi6r + psir[i];
        accumpsi6i = accumpsi6i + psii[i];
    }
    accumpsi6r = accumpsi6r / numP;
    accumpsi6i = accumpsi6i / numP;
    psi6 = sqrt(accumpsi6r * accumpsi6r + accumpsi6i * accumpsi6i);
    return psi6;
}

void Model::getPermutator() {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                eps[i][j][k] == 0;
                if (i == 0 && j == 1 && k == 2) {
                    eps[i][j][k] = 1;
                }
                if (i == 1 && j == 2 && k == 0) {
                    eps[i][j][k] = 1;
                }
                if (i == 2 && j == 0 && k == 1) {
                    eps[i][j][k] = 1;
                }
                if (i == 2 && j == 1 && k == 0) {
                    eps[i][j][k] = -1;
                }
                if (i == 0 && j == 2 && k == 1) {
                    eps[i][j][k] = -1;
                }
                if (i == 1 && j == 0 && k == 2) {
                    eps[i][j][k] = -1;
                } 
            }
        }
    }
}


void Model::readObstacle(){
    std::ifstream is;
    is.open(parameter.obstacleFilename);

    std::string line;
    double dum;
    double r[3];
    while (getline(is, line)){
        std::stringstream linestream(line);
        
        linestream >> dum;
        linestream >> r[0];
        linestream >> r[1];
        linestream >> r[2];
        obstacles.push_back(std::make_shared<pos>(r[0]*radius,r[1]*radius,r[2]*radius));
    } 
    numObstacles = obstacles.size();
    is.close();
}

