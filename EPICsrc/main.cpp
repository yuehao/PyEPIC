#include <iostream>
#include"beam.h"
#include<random>
#include<mpi.h>

using namespace std;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double betx=0.26,bety=0.13;
  double alfx=0.0,alfy=0.0;
  double ex=3.41e-9,ey=3.41e-9;
  double sigx=sqrt(betx*ex), sigy=sqrt(bety*ey), sigz=0.16, sigd=1.0e-4;
  double E0=15e9, E1=275e9, m0=0.511e6, m1=938.0e6;
  COneTurnMap mapx(betx, alfx, 0.675,0);
  COneTurnMap mapy(bety, alfy, 0.685,0);
  Crf rf(0,-1,0,0,0);

  std::random_device seed;
  std::mt19937 mt(seed());

  beam b0,b1;
  b0.SetN(3.0e11, 100000);b1.SetN(3.0e11, 100000);
  b0.SetSize(sigx,sigy,sigz,sigd);b1.SetSize(sigx,sigy,sigz,sigd);
  b0.SetEmit(ex,ey);b1.SetEmit(ex,ey);
  b0.SetEnergy(E0,m0,-1);b1.SetEnergy(E1,m1,1);
  b0.SetSlice(4);b1.SetSlice(32);

  b0.Ini6D(mapx, mapy, rf, mt);
  b1.Ini6D(mapx, mapy, rf, mt);

  double t0=MPI_Wtime();
  b0.set_longitudinal_slices2();
  b1.set_longitudinal_slices2();
  double lum=beambeam_pass(b1,b0);
  double t1=MPI_Wtime();
  if(rank==0){
	std::cout<<t1-t0<<std::endl;
	std::cout<<lum<<std::endl;
  }

  b0.print_coord("hello");



  MPI_Finalize();
  return 0;
}
