//
// Created by Yue Hao on 5/26/16.
//

#include "acc_model.h"
#include "mathfunc.h"
#include <cmath>
#include <iostream>
#include <cassert>

void COneTurnMap::SetMatrix(const double &m11, const double &m12, const double &m21, const double &m22) {
    this->m11=m11;
    this->m12=m12;
    this->m21=m21;
    this->m22=m22;
}

void COneTurnMap::Pass(double &x, double &px, const double &de) const {
    double tx=x,tpx=px;
    x=m11*tx+m12*tpx;
    px=m21*tx+m22*tpx;
    if (this->chrom!=0){
        double sintemp=sin(2.0*M_PI*this->chrom*de), costemp=cos(2.0*M_PI*this->chrom*de);
        double tx=x,tpx=px;
        x=(costemp+alpha*sintemp)*tx+(beta*sintemp)*tpx;
        px=(-gamma*sintemp)*tx+(costemp-alpha*sintemp)*tpx;
    }
    return;


}
COneTurnMap::COneTurnMap(const double &beta, const double &alpha, const double &tune, const double &chrom) :
        beta(beta), alpha(alpha), tune(tune), chrom(chrom)
{
    this->gamma=(1+alpha*alpha)/beta;
    double angle=2.0*M_PI*tune;

    m11=cos(angle)+alpha*sin(angle); m12=beta*sin(angle);
    m21=-gamma*sin(angle); m22=cos(angle)-alpha*sin(angle);
}

CCrabCavity::CCrabCavity():harmonic(1), harm_ratio(0), freq(400e6), rms_voltage_jitter(0.0), rms_phase_jitter(0.0),
                           use_transfer_map(0), amp_factor(1.0) , before_bb(1){
    this->kcc = 2.0 * M_PI * freq / clight;
    this->crab_strength_IP = 0;
    
}

CCrabCavity::CCrabCavity(const int& bbb, const double &frequency, const double &half_crossing_angle, const double& deviation, const int& harm, const double& ratio) :
harmonic(harm), harm_ratio(ratio), freq(frequency), rms_voltage_jitter(0), rms_phase_jitter(0),
use_transfer_map(0), amp_factor(1.0) , before_bb(bbb){
    this->kcc = 2.0 * M_PI * frequency / clight;
    this->crab_strength_IP = half_crossing_angle * deviation / kcc;
    //std::cout << this->crab_strength_IP<<std::endl;
 
}

void CCrabCavity::set_error(const double &voltage_jitter, const double &phase_jitter) {
    this->rms_phase_jitter = phase_jitter;
    this->rms_voltage_jitter = voltage_jitter;
}
void CCrabCavity::set_optics(const double &betax, const double &alphax, const double &etax, const double &etapx,
                             const double &phase_adv, const double& b1, const double& a1, const double& eta1, const double& etap1) {
    use_transfer_map = 1;

    this->betax = betax;
    this->alphax = alphax;
    this->etax = etax;
    this->etapx = etapx;
    this->phase_adv = phase_adv;
    double b2 = betax, a2 = alphax;
    amp_factor=sqrt(b1*b2);
    m11 = sqrt(b2 / b1) * (cos(phase_adv) + a1 * sin(phase_adv));
    m12 = amp_factor * sin(phase_adv);
    m14 = etax - m11 * eta1 - m12 * etap1;
    m21 = (-(1.0 + a1 * a2) * sin(phase_adv) + (a1 - a2) * cos(phase_adv)) / amp_factor;
    m22 = sqrt(b1 / b2) * (cos(phase_adv) - a2 * sin(phase_adv));
    m24 = etapx - m21 * eta1 - m22 * etap1;
    m31 = m14 * m21 - m24 * m11;  //checked equiv m21 * etax - m11 * etapx + etap1
    m32 = m14 * m22 - m24 * m12;    //checked equiv m22 * etax - m12 * etapx - eta1
    
    im11 = m22;
    im12 = -m12;
    im21 = -m21;
    im22 = m11;
    im14 = -m14 * m22 + m24 * m12;
    im24 = m14 * m21 - m24 * m11;
    im31 = m24;
    im32 = -m14;
    //std::cout<<m11<<'\t'<<m12<<'\t'<<m21<<'\t'<<m22<<std::endl;
    //std::cout<<im11<<'\t'<<im12<<'\t'<<im21<<'\t'<<im22<<std::endl;
}

void CCrabCavity::pass(double &x, double& px, double& z, double &de) const {
    if (use_transfer_map>0) {
        double tempx = x, temppx = px;
        x = m11 * tempx + m12 * temppx + m14 * de;
        px = m21 * tempx + m22 * temppx + m24 * de;
        z += m31 * tempx + m32 * temppx;
    
        double dpx = crab_strength_IP * sin(kcc * z) / amp_factor;
        double dde = -crab_strength_IP * kcc * cos(kcc * z) * x / amp_factor;
        if (harmonic > 0 && harm_ratio != 0) {
            dpx *= (1 - harm_ratio);
            dde *= (1 - harm_ratio);
            dpx += harm_ratio * crab_strength_IP * sin(harmonic * kcc * z) / harmonic / amp_factor;
            dde -= harm_ratio * crab_strength_IP * kcc * cos(harmonic * kcc * z) * x / amp_factor;
        }
        px += dpx;
        de += dde;

        tempx = x;
        temppx = px;
        
        x = im11 * tempx + im12 * temppx + im14 * de;
        px = im21 * tempx + im22 * temppx + im24 * de;
        z += im31 * tempx + im32 * temppx;
    
    }
    else {
        double dx = -crab_strength_IP * sin(kcc * z) * before_bb ;
        double dde = crab_strength_IP * kcc * cos(kcc * z) * px * before_bb ;
        if (harmonic > 0 && harm_ratio != 0) {
            dx *= (1 - harm_ratio);
            dde *= (1 - harm_ratio);
            dx -= harm_ratio * crab_strength_IP * sin(harmonic * kcc * z) * before_bb / harmonic;
            dde += harm_ratio * crab_strength_IP * kcc * cos(harmonic * kcc * z) * px * before_bb;
        }
        x += dx;
        de += dde;
    }
    
}


/*************************************************************************************************************/
lattice_radiation_property & lattice_radiation_property::SetEmit(double ex0, double ey0){
  ex=ex0;ey=ey0;
  return *this;
}
lattice_radiation_property & lattice_radiation_property::SetDamping(double nx, double ny, double nz){
  assert(nx>0.0 && ny>0.0 && nz>0.0);

  damping_turns_x=nx;
  damping_turns_y=ny;
  damping_turns_z=nz;
  damping_strength_x=exp(-1.0/damping_turns_x);
  damping_strength_y=exp(-1.0/damping_turns_y);
  damping_strength_z=exp(-1.0/damping_turns_z);
  excitation_strength_x=sqrt(1.0-damping_strength_x*damping_strength_x);
  excitation_strength_y=sqrt(1.0-damping_strength_y*damping_strength_y);
  excitation_strength_z=sqrt(1.0-damping_strength_z*damping_strength_z);

  is_damping=true;
  is_excitation=true;
  return *this;
}
lattice_radiation_property & lattice_radiation_property::SetLongitudinal(double z, double delta){
  zsize=z;
  deltaE=delta;
  return *this;
}
lattice_radiation_property & lattice_radiation_property::SetTransverse(const COneTurnMap &mapx, const COneTurnMap &mapy){
  assert(ex>0.0 && ey>0.0);
  xsize=sqrt(mapx.beta*ex);
  pxsize=sqrt(mapx.gamma*ex);
  ysize=sqrt(mapy.beta*ey);
  pysize=sqrt(mapy.gamma*ey);
  return *this;
}


