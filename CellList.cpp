#include "CellList.h"

CellList::CellList(double cutDist, int dim0, int maxCount0, double box_x0,
        double box_y0,double box_z0) :
cutDistance(cutDist), dim(dim0), maxCount(maxCount0), boxSize_x(box_x0),
boxSize_y(box_y0), boxSize_z(box_z0){

    del_x = cutDistance;
    del_y = cutDistance;
    del_z = cutDistance;
    nbin_x = (int) (boxSize_x / del_x) + 1;
    nbin_y = (int) (boxSize_y / del_y) + 1;
    nbin_z = (int) (boxSize_z / del_z) + 1;
    if (dim == 2) {
        nbin_z = 1;
    }
    nbin = nbin_x * nbin_y * nbin_z;
    min_x = -0.5 * boxSize_x;
    min_y = -0.5 * boxSize_y;
    min_z = -0.5 * boxSize_z;
    this->setup();
}

CellList::Idx_3d CellList::coordToIdx(double x, double y, double z) const{

    int z_bin, x_bin, y_bin;

    if (dim == 2) {
        z_bin = 0;
        x_bin = (int) ((x - min_x) / del_x);
        y_bin = (int) ((y - min_y) / del_y);
    } else {
        x_bin = (int) ((x - min_x) / del_x);
        y_bin = (int) ((y - min_y) / del_y);
        z_bin = (int) ((z - min_z) / del_z);
    }
    Idx_3d idx(z_bin, y_bin, x_bin);
    return idx;
}

std::vector<int> CellList::getNeighbors(double x, double y, double z) {

    Idx_3d idx = coordToIdx(x,y,z);
    std::vector<int> res;
    int idx_x0 = std::get<2>(idx);
    int idx_y0 = std::get<1>(idx);
    int idx_z0 = std::get<0>(idx);
    int nb_idx = (*oneDIdx)[idx_z0][idx_y0][idx_x0];
    for (auto& cellIdx : cellNeighborIdxList[nb_idx]) {
        int idx_x = std::get<2>(cellIdx);
        int idx_y = std::get<1>(cellIdx);
        int idx_z = std::get<0>(cellIdx);
   
        int count = (*cellListCount)[idx_z][idx_y][idx_x];
        for (int i = 0; i < count; i++) {
            res.push_back((*cellList)[idx_z][idx_y][idx_x][i]);
        }
    }
    return res;
}
// build cellList

int CellList::buildList(const Model::state &s) {
    int totalCount = 0;
//    threeDIdx.clear();
    for(int i = 0; i < nbin_z; i++){
        for (int j = 0; j < nbin_y; j++){
            for(int k = 0; k < nbin_x; k++){
                (*cellListCount)[i][j][k] = 0;
        }
    }
    }
    for (int i = 0; i < s.size(); i++) {
        Idx_3d idx = coordToIdx(s[i]->r[0], s[i]->r[1], s[i]->r[2]);

        int idx_x = std::get<2>(idx);
        int idx_y = std::get<1>(idx);
        int idx_z = std::get<0>(idx);
        
        if(idx_x >= nbin_x) idx_x=nbin_x-1;
        if(idx_y >= nbin_y) idx_y=nbin_y-1;
        if(idx_z >= nbin_z) idx_z=nbin_z-1;
        if(idx_x < 0 ) idx_x = 0;
        if(idx_y < 0 ) idx_y = 0;
        if(idx_z < 0 ) idx_z = 0;
        
        
        (*cellList)[idx_z][idx_y][idx_x][(*cellListCount)[idx_z][idx_y][idx_x]] = i;
        (*cellListCount)[idx_z][idx_y][idx_x]++;        
        totalCount++;
    }
    return totalCount;
//    std::cout <<totalCount<< "particle binned" << std::endl;
}

int CellList::buildList(const Model::posArray &s) {
    int totalCount = 0;
//    threeDIdx.clear();
    for(int i = 0; i < nbin_z; i++){
        for (int j = 0; j < nbin_y; j++){
            for(int k = 0; k < nbin_x; k++){
                (*cellListCount)[i][j][k] = 0;
        }
    }
    }
    for (int i = 0; i < s.size(); i++) {
        Idx_3d idx = coordToIdx(s[i]->r[0], s[i]->r[1], s[i]->r[2]);

        int idx_x = std::get<2>(idx);
        int idx_y = std::get<1>(idx);
        int idx_z = std::get<0>(idx);
        
        if(idx_x >= nbin_x) idx_x=nbin_x-1;
        if(idx_y >= nbin_y) idx_y=nbin_y-1;
        if(idx_z >= nbin_z) idx_z=nbin_z-1;
        if(idx_x < 0 ) idx_x = 0;
        if(idx_y < 0 ) idx_y = 0;
        if(idx_z < 0 ) idx_z = 0;
        
        
        (*cellList)[idx_z][idx_y][idx_x][(*cellListCount)[idx_z][idx_y][idx_x]] = i;
        (*cellListCount)[idx_z][idx_y][idx_x]++;        
        totalCount++;
    }
    return totalCount;
//    std::cout <<totalCount<< "particle binned" << std::endl;
}

void CellList::setup() {

    std::array<Array4D_type::index, 4> dim1 = {nbin_z, nbin_y, nbin_x, maxCount};
    std::array<Array4D_type::index, 3> dim2 = {nbin_z, nbin_y, nbin_x};
    cellList = std::make_shared<Array4D_type>(dim1);
    oneDIdx = std::make_shared<Array3D_type>(dim2);
    cellListCount = std::make_shared<Array3D_type>(dim2);

    int count = 0;
    // now build the cellNeighborIdxList
    for (int i = 1; i < nbin_x - 1; i++) {
        for (int j = 1; j < nbin_y - 1; j++) {
            if (dim == 2) {
                cellNeighborIdxList.push_back(std::list<Idx_3d>());
                (*oneDIdx)[0][j][i] = count;   
                threeDIdx.push_back(std::make_tuple(0,j,i));
                for (int ii = -1; ii < 2; ii++) {
                    for (int jj = -1; jj < 2; jj++) {
                        int idx_x = i + ii;
                        int idx_y = j + jj;
                        int idx_z = 0;
                        Idx_3d idx(idx_z, idx_y, idx_x);

                        cellNeighborIdxList[count].push_back(idx);
                        
                    }
                }
                count++;
            } else {

                for (int k = 1; k < nbin_z - 1; k++) {
                    cellNeighborIdxList.push_back(std::list<Idx_3d>());
                    (*oneDIdx)[k][j][i] = count;
                    threeDIdx.push_back(std::make_tuple(k,j,i));
                    for (int ii = -1; ii < 2; ii++) {
                        for (int jj = -1; jj < 2; jj++) {
                            for (int kk = -1; kk < 2; kk++) {
                                int idx_x = i + ii;
                                int idx_y = j + jj;
                                int idx_z = k + kk;
                                Idx_3d idx(idx_z, idx_y, idx_x);

                                cellNeighborIdxList[count].push_back(idx);

                            }
                        }
                    }
//                    std::cout << "nb size: " << count << "\t" << cellNeighborIdxList[count].size() << std::endl;
                    count++;
                }
                
                
            }
        }
    }

}

void CellList::printCellList() const{
    for(int nb_idx = 0; nb_idx < nbin; nb_idx++){
        std::cout << "cell idx:" << "\t" << std::get<2>(threeDIdx[nb_idx])
                << "\t" << std::get<1>(threeDIdx[nb_idx]) << "\t"
                << std::get<0>(threeDIdx[nb_idx]) << std::endl;
        for (auto& cellIdx : cellNeighborIdxList[nb_idx]) {
            int idx_x = std::get<2>(cellIdx);
            int idx_y = std::get<1>(cellIdx);
            int idx_z = std::get<0>(cellIdx);
            int count = (*cellListCount)[idx_z][idx_y][idx_x];
            std::cout << "neighborlist: " <<idx_x << " \t" << idx_y << "\t" << idx_z << std::endl;
        }
    }
}


void CellList::printParticleList() const{
    for(int i = 0; i < threeDIdx.size(); i++){
        Idx_3d cellIdx = threeDIdx[i];
        int idx_x = std::get<2>(cellIdx);
        int idx_y = std::get<1>(cellIdx);
        int idx_z = std::get<0>(cellIdx);
        std::cout << i<< "\t" <<idx_x << "\t" << idx_y << "\t" << idx_z << std::endl;
       
    }
}


void CellList::getParticleIdx(double x, double y, double z, int idx[3]) const{

    Idx_3d idx2 = coordToIdx(x,y,z);
    idx[0] = std::get<2>(idx2);
    idx[1] = std::get<1>(idx2);
    idx[2] = std::get<0>(idx2);
}


void CellList::printCellContent(int idx[3]) const{
  
    int count = (*cellListCount)[idx[2]][idx[1]][idx[0]];
        for (int i = 0; i < count; i++) {
            std::cout << (*cellList)[idx[2]][idx[1]][idx[0]][i] << std::endl;
        }
    



}