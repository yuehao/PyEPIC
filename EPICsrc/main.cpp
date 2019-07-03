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
  double alfx=-1.0,alfy=1.0;
  double ex=3.41e-9,ey=3.41e-9;
  double E0=15e9, E1=275e9, m0=0.511e6, m1=938.0e6;
  double nux=0.675, nuy=0.685, nuz=0.0043;
  double chromx=0.0, chromy=0.0;
  double N0=3.0e11, N1=3.0e11;

  COneTurnMap mapx(betx, alfx, nux, chromx);
  COneTurnMap mapy(bety, alfy, nuy, chromy);
  Crf rf(nuz);
  lattice_radiation_property rad;
  rad.SetEmit(ex, ey);
  rad.SetDamping(500, 500, 500);
  rad.SetLongitudinal(0.16, 1.0e-4);
  rad.SetTransverse(mapx, mapy);

  std::random_device seed;
  std::mt19937 mt(seed());

  /*
  beam b0,b1;
  b0.SetN(N0, 100000);b1.SetN(N1, 100000);
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
  */


  beam b0;
  b0.SetN(N0,100000);
  b0.SetSize(10*rad.xsize,5*rad.ysize,3*rad.zsize,6*rad.deltaE);
  b0.Ini6D(mapx,mapy,rf,mt);
  int turns=11;
  if(argc>1)
	turns=std::stoi(argv[1]); 
  while(--turns>0){
	b0.OneTurn2(mapx,mapy,rf,rad,mt);
	if(rank==0)
	  std::cout<<turns<<std::endl;
  }
  b0.print_coord("hello");


  MPI_Finalize();
  return 0;
}
