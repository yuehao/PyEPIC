//
// Created by Yue Hao on 5/26/16.
//
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include "acc_model.h"

#ifndef EPIC_BEAM_H
#define EPIC_BEAM_H


class beam {
public:
  friend double beambeam_pass(beam &, beam &);
  friend double beambeam_slicepass(beam &, unsigned, std::vector<std::vector<double>> &, beam &, unsigned, std::vector<std::vector<double>> &);
  void print_pdf(const std::string &, const double &) const;
  void print_coord(const std::string &, int=-1) const;
  void set_longitudinal_slices2();
private:
  std::vector<std::vector<unsigned long>> par_inslice;
private:
    //Should be in dimension of n_macro
    std::vector<double> x_;
    std::vector<double> px_;
    std::vector<double> y_;
    std::vector<double> py_;
    std::vector<double> z_;
    std::vector<double> delta_;
    
    std::vector<int> inslice; //The main slice
    std::vector<int> adslice; //The partial slice
    std::vector<double> partition; //The part in main slice
    
    //Should be in dimension of zslice
    std::vector<double> slice_center;
    std::vector<double> slice_npar;  //Particle count
    std::vector<double> slice_xc;
    std::vector<double> slice_yc;
    std::vector<double> slice_xrms;
    std::vector<double> slice_yrms;
    
    //std::vector<double> slice_npar_gaussian;
    std::vector<double> slice_crab_deviation;
    //std::vector<double> slice_xrms_gaussian;
    //std::vector<double> slice_yrms_gaussian;
    
    //Beam boundaries
    double zmax, zmin;
    double xmax, xmin, ymax, ymin;
    
    //Beam statistics
    double xsize_ini, ysize_ini, zsize_ini, deltaE_ini;
    double xsize, ysize, zsize, deltaE;
    double pxsize, pysize;
    double xemit, yemit;
    double xmean, ymean, pxmean, pymean, zmean, demean;
    double xfourth, yfourth;
    
    //General things
    double classrad; //Classical radius
    double aperture;//Aperture

public:
    friend class beambeam;
    beam();
    
    ~beam() {}
    
    double n_par; //Bunch Intensity
    unsigned long n_macro, n_live; //Number of macro_particle, Nr still alive
    double energy, mass, gamma_e, charge;
    unsigned zslice;
    COneTurnMap mapz;
    int formfactor;
    
    
    void SetN(double Npar, unsigned long Nm);
    
    void SetSize(double xs, double ys, double zs, double delta);
    
    void SetEmit(double xemit, double yemit);
    
    void SetEnergy(double energy, double m, double charge);
    
    void SetSlice(unsigned zslice){
        this->zslice = zslice;
    }
    
    void Ini6D(const COneTurnMap &mapx, const COneTurnMap &mapy, const Crf &rfz, std::mt19937& rdgen);
    
    void OneTurn(const COneTurnMap &mapx, const COneTurnMap &mapy, const Crf &rfz);
    
    void findzBoundary();
    
    void CalEmit();
    
    void Drift(const double &length);
    
    void Fold(const double &ratio);
    
    void QuadKick(const double &dx, const double &dy);
    
    void PosKick(const double &dx, const double &dpx, const double &dy, const double &dpy);
    
    void LorentzTransform(const double &half_cross_angle, const int &forward = 1);
    
    void crab_kick(const CCrabCavity &crabcavity);
    
    void set_longitudinal_slices();
    
    void simplified_crab_deviation(const CCrabCavity &cc, const double &hca);
    
    double
    weak_strong_Gaussian(const double &charge, const double &npar, const double &cx, const double &cy, const double &sx,
                       const double &sy);
    
    friend std::vector<double> beambeam_2pass(beam &longbeam, beam &shortbeam, int output=0);
    
};

int Efield_Gaussian(double xsize, double ysize, double x, double y,
                    double& Ex,double& Ey, const int& linear=0);

#endif //EPIC_BEAM_H
