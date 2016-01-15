#include "CellList.h"
#include "model.h"

void test(double x, double y, double z);
int main(){
    int N = 2;
    Model::state s;
    for (int i = 0; i < N; i++){
        s.push_back(Model::particle_ptr(new Model::particle));
        s[i]->r[0] = i - N/2.0;
        s[i]->r[1] = i - N/2.0;
        s[i]->r[2] = i - N/2.0;
    }
    for (int i = 0; i < s.size(); i++) {
        std::cout <<s[i]->r[0] << "\t"<< s[i]->r[1] <<"\t" <<s[i]->r[2] << std::endl;
    }
    
    
    CellList cell(1,2,10,4,4,4);
    cell.buildList(s);
//    cell.printParticleList();
    cell.printCellList();
    
    for (int i = 0; i < N; i++){
        std::vector<int> idx = cell.getNeighbors(i-N/2.0,i-N/2.0,i-N/2.0);
        std::cout << i << "\t";
        for (int j = 0; j < idx.size(); j++){
            std::cout << idx[j] << "\t"; 
        }
        std::cout << std::endl;
    }
//    cell.buildList(s);
    return 0;
}

