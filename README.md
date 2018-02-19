# Eremitic
Non-hydrodynamic system dynamics

Needs: Eigen library
Needs: gsl library
Nice but not necessary: openmp library

Contains: 
compile.sh: script for compiling executables
MC-xx: Monte Carlo sampling for generating initial hot-spot locations from nuclear geometries
Hermit-xx: Numerical integration of first order eremitic result from inital hot spot locations

Executables:
./MC-pPb.exe: generates hot-spot locations "nucpos_xxx" for central pPb collisions in the folder "out"
./MC-PbPb.exe: generates hot-spot locations "nucpos_xxx" for 30-40% PbPb collisions in the folder "out"

./Hermit-3d.exe params-3d.txt: generates results obtained from zeroth/first order eremitic expansion f1 in folder "data" using the run parameters specified in "params-3d.txt"
./Hermit-BoostInv.exe params-BoostInv.txt nucpos_xxx.dat: generates results obtained from zeroth/first order eremitic expansion f1 assuming boost invariance in folder "data" using the run parameters specified in "params-BoostInv.txt"
