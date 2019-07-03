//
// Created by Yue Hao on 5/26/16.
//

#ifndef EPIC_RING_MODEL_H
#define EPIC_RING_MODEL_H
//enum Twiss {Beta, Alpha, Gamma, Tune, Chrom};
#include <vector>

class COneTurnMap{

public:
    double m11,m12,m21,m22;
    double beta;
    double alpha;
    double gamma;
    double tune;
    double chrom;
    COneTurnMap(){}
    COneTurnMap(const double& beta, const double& alpha ,const double& tune, const double& chrom);
    void Pass(double& x, double& px, const double& de) const;
    void SetMatrix(const double& m11, const double& m12,const double& m21, const double& m22);
};



class Crf{
public:
    double gammat;
	//07/02/2019:provide a inclass initializer for voltage
    double voltage=-1.0;
    double phase;
    double harm;
    double freq0;
    int keep_z;
    double tune_syn;
    Crf(){}
	//07/02/2019:one parater constructor
	Crf(double nus):tune_syn(nus){}
    Crf(const double& gmt, const double& v, const double& ph, const double& h, const double& freq, const int& keep_z=1):
            gammat(gmt),voltage(v), phase(ph), harm(h), freq0(freq),keep_z(keep_z), tune_syn(0.0){};
    void force_tune(double nus) {tune_syn=nus;}
};

class CCrabCavity{
private:
    double kcc;
    double crab_strength_IP;
    int before_bb;
    int use_transfer_map;
    double m11,m12,m21,m22,m14,m24,m31,m32;
    double im11,im12,im21,im22,im14,im24,im31,im32;
    double amp_factor;
public:
    double freq;
    double harmonic;
    double harm_ratio;
    double rms_voltage_jitter;
    double rms_phase_jitter;
    double betax;
    double alphax;
    double etax;
    double etapx;
    double phase_adv;
    
    
    CCrabCavity();
    CCrabCavity(const int& bbb, const double& frequency, const double& half_crossing_angle, const double& deviation=1.0, const int& harm=3, const double& ratio=0);
    void set_optics(const double& betax, const double& alphax, const double& etax, const double& etapx, const double& phase_adv,
                    const double& beta_star, const double& alpha_star, const double& eta_star, const double& etap_star);
    void set_error(const double& voltage_jitter, const double& phase_jitter);
    void pass(double& x, double& px, double& z, double& de) const;
    //std::vector<double> IP_deviation(std::vector<double>& zpos, const double& half_crossing_angle) const;
};

//07/02/2019
struct lattice_radiation_property{
  lattice_radiation_property &SetEmit(double, double);
  lattice_radiation_property &SetDamping(double, double, double);
  lattice_radiation_property &SetLongitudinal(double, double);
  lattice_radiation_property &SetTransverse(const COneTurnMap&, const COneTurnMap&); //Must be called after "SetEmit"

  double ex=-1.0, ey=-1.0;
  double zsize=-1.0, deltaE=-1.0;
  double damping_turns_x=-1.0, damping_turns_y=-1.0, damping_turns_z=-1.0; // 1/e damping time in unit of turns
  double damping_strength_x, damping_strength_y, damping_strength_z;
  double excitation_strength_x, excitation_strength_y, excitation_strength_z;
  double xsize=-1.0, ysize=-1.0, pxsize=-1.0, pysize=-1.0;
  bool is_damping=false, is_excitation=false;
};


#endif //EPIC_RING_MODEL_H

