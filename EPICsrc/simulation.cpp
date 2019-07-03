#include <iostream>
#include"beam.h"
#include<random>
#include<mpi.h>
#include<cmath>
#include<fstream>

using namespace std;
using std::ofstream;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::random_device seed;
  std::mt19937 mt(seed());

  beam electron, proton;
  double half_crossing_angle=12.5e-3;
  unsigned long nMacros_ele=51200;
  unsigned long nMacros_pro=204800;

  /*************************************************************************************************************************************/
  /* set up parameters for electron*/
  double beta_ele_IP_x=0.63, beta_ele_IP_y=0.104, beta_ele_cc_x=800.0;
  //double L_ele_from_IP_to_cc=sqrt((beta_ele_cc_x-beta_ele_IP_x)*beta_ele_IP_x);
  //double phase_ele_cc_x=atan(L_ele_from_IP_to_cc/beta_ele_IP_x); 
  //double alpha_ele_cc_x=L_ele_from_IP_to_cc/beta_ele_IP_x;
  double phase_ele_cc_x=M_PI/2.0, alpha_ele_cc_x=0.0;

  electron.SetSize(112.0e-6,25.0e-6,1.9e-2,5.5e-4); //sigma_x, sigma_y, sigma_z, sigma_delta
  electron.SetN(3.0e11,nMacros_ele);
  electron.SetEnergy(10.0e9,0.510997e6,-1.0);
  electron.SetSlice(4);

  COneTurnMap map_ele_x(beta_ele_IP_x, 0.0, 0.08, 1.0); //bet, alf, nu, chrom
  COneTurnMap map_ele_y(beta_ele_IP_y, 0.0, 0.06, 1.0); //bet, alf, nu, chrom
  Crf rf_ele(6.9e-2); //nu_z

  lattice_radiation_property rad_ele;
  rad_ele.SetEmit(20.0e-9, 4.9e-9);
  rad_ele.SetDamping(4000., 4000., 2000.);
  rad_ele.SetLongitudinal(1.9e-2, 1.0e-4);
  rad_ele.SetTransverse(map_ele_x, map_ele_y);

  CCrabCavity cc_ele_b(1,563.0e6,-half_crossing_angle,1.0,-1); //is before collision, frequency, half crossing angle, deviation, harmonic
  cc_ele_b.set_optics(beta_ele_cc_x, -alpha_ele_cc_x, 0.0, 0.0, -phase_ele_cc_x, //at cc: beta, alpha, eta, etap, phase advance to IP
	  beta_ele_IP_x, 0.0, 0.0, 0.0); //at IP: beta, alpha, eta, etap
  CCrabCavity cc_ele_a(-1,563.0e6,-half_crossing_angle,1.0,-1); //is before collision, frequency, half crossing angle, deviation, harmonic
  cc_ele_a.set_optics(beta_ele_cc_x, alpha_ele_cc_x, 0.0, 0.0, phase_ele_cc_x, //at cc: beta, alpha, eta, etap, phase advance to IP
	  beta_ele_IP_x, 0.0, 0.0, 0.0); //at IP: beta, alpha, eta, etap

  /*************************************************************************************************************************************/
  /* set up parameters for proton */
  double beta_pro_IP_x=0.90, beta_pro_IP_y=0.059, beta_pro_cc_x=800.0;
  //double L_pro_from_IP_to_cc=sqrt((beta_pro_cc_x-beta_pro_IP_x)*beta_pro_IP_x);
  //double phase_pro_cc_x=atan(L_pro_from_IP_to_cc/beta_pro_IP_x); 
  //double alpha_pro_cc_x=L_pro_from_IP_to_cc/beta_pro_IP_x;
  double phase_pro_cc_x=M_PI/2.0, alpha_pro_cc_x=0.0;

  proton.SetSize(112.0e-6,22.5e-6,0.07,6.5e-4);
  proton.SetN(1.05e11, nMacros_pro);
  proton.SetEnergy(275.0e9,938.27231e6,1.0);
  proton.SetSlice(32);

  COneTurnMap map_pro_x(beta_pro_IP_x, 0.0, 0.18 , 2.0); //bet, alf, nu, chrom
  COneTurnMap map_pro_y(beta_pro_IP_y, 0.0, 0.175, 2.0); //bet, alf, nu, chrom
  Crf rf_pro(0.005);
  
  lattice_radiation_property rad_pro; //default: no damping and excitation

  CCrabCavity cc_pro_b(1,200.0e6,-half_crossing_angle,1.0,-1); //is before collision, frequency, half crossing angle, deviation, harmonic
  cc_pro_b.set_optics(beta_pro_cc_x, -alpha_pro_cc_x, 0.0, 0.0, -phase_pro_cc_x, //at cc: beta, alpha, eta, etap, phase advance to IP
	  beta_pro_IP_x, 0.0, 0.0, 0.0); //at IP: beta, alpha, eta, etap
  CCrabCavity cc_pro_a(-1,200.0e6,-half_crossing_angle,1.0,-1); //is before collision, frequency, half crossing angle, deviation, harmonic
  cc_pro_a.set_optics(beta_pro_cc_x, alpha_pro_cc_x, 0.0, 0.0, phase_ele_cc_x, //at cc: beta, alpha, eta, etap, phase advance to IP
	  beta_pro_IP_x, 0.0, 0.0, 0.0); //at IP: beta, alpha, eta, etap

  /*************************************************************************************************************************************/
  /* generate initial distribution */
  electron.Ini6D(map_ele_x, map_ele_y, rf_ele, mt);
  proton.Ini6D(map_pro_x, map_pro_y, rf_pro, mt);

  /* track */
  unsigned turns=50000,index=0,output_lines=100,current=0;
  vector<double> luminosity(turns);
  ofstream out;
  if(rank==0)
	out.open("luminosity");

  double t0=MPI_Wtime();
  while(index!=turns){
	electron.crab_kick(cc_ele_b);
	electron.LorentzTransform(half_crossing_angle, 1);
	proton.crab_kick(cc_pro_b);
	proton.LorentzTransform(half_crossing_angle,1);

	electron.set_longitudinal_slices2();
	proton.set_longitudinal_slices2();
	luminosity[index]=beambeam_pass(proton,electron);

	electron.LorentzTransform(half_crossing_angle, -1);
	electron.crab_kick(cc_ele_a);
	proton.LorentzTransform(half_crossing_angle, -1);
	proton.crab_kick(cc_pro_a);

	electron.OneTurn2(map_ele_x, map_ele_y, rf_ele, rad_ele, mt);
	proton.OneTurn2(map_pro_x, map_pro_y, rf_pro, rad_pro, mt);

	electron.findzBoundary();
	proton.findzBoundary();

	if(rank==0 && (index+1)%output_lines==0){
	  for(unsigned line=0;line<output_lines;++line)
		out<<luminosity[current+line]<<"\n";
	  current+=output_lines;
	  double t1=MPI_Wtime();
	  std::cout<<current/output_lines<<"*" << output_lines<<" turns:\t"<<t1-t0<<" seconds."<<std::endl;
	  t0=t1;
	}

	++index;
  }

  if(rank==0)
	out.close();

  MPI_Finalize();
  return 0;
}
