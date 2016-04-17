CC = gcc
CXX = g++

HOME=/home/yuguangyang/
ARMA_INCLUDE=-I$(HOME)Downloads/armadillo-5.100.2/include
RL_INCLUDE=-I../../include
GTEST_INCLUDE=-I$(HOME)workspace/libs/gtest-1.7.0/include
BOOST_INCLUDE=-I/opt/boost/boost_1_57_0
PROTO_INCLUDE=-I/usr/local/include

GTEST_PATH=-L$(HOME)workspace/libs/gtest-1.7.0/mybuilds
RL_PATH=-L../../lib
PROTO_PATH=-L/usr/local/lib

DEBUGFLAG=-DDEBUG -g3
RELEASEFLAG= -O3 -march=native -DARMA_NO_DEBUG
CXXFLAGS=  -std=c++0x -I$(HOME)Copy/workspace/munkres-cpp/src $(BOOST_INCLUDE) -D__LINUX  -I./libicp $(ARMA_INCLUDE) -DARMA_DONT_USE_WRAAPER
#CXXFLAGS += $(DEBUGFLAG)
#CXXFLAGS += $(RELEASEFLAG)
LINKOPTFLAGS= -O3 -flto=4 -fwhole-program
LDFLAG= -L./libicp -licp -llapack -lblas -lgfortran -lquadmath -pthread

OBJ=HungarianAlg.o controller.o model.o simulator.o Driver.o CellList.o
OBJ2=cellList_test.o CellList.o
all:test.exe test_cell.exe
test.exe: $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAG) 
	
test_static.exe: $(OBJ)
	$(CXX) -o $@ $^ -static $(LDFLAG) 
test_cell.exe:$(OBJ2)
	$(CXX) -o $@ $^ $(LDFLAG) 

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $(RELEASEFLAG) $^

clean:
	rm *.o *.exe
