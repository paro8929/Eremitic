g++ -I. -fopenmp -std=c++11  MC-PbPb.cpp -o MC-PbPb.exe
g++ -I. -fopenmp -std=c++11  MC-pPb.cpp -o MC-pPb.exe
g++ -I. -fopenmp -I Eigen Hermit-3d.cpp  -lgsl -lgslcblas -o Hermit-3d.exe
g++ -I. -fopenmp -I Eigen Hermit-BoostInv.cpp  -lgsl -lgslcblas -o Hermit-BoostInv.exe
