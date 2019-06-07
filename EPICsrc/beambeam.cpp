//
// Created by Yue Hao on 7/29/16.
//

#include "beambeam.h"
void beambeam::set_beam_param(int beamid, double np, unsigned long nm, double sx, double sy, double sz, double sdelta) {
    if (beamid==1){
        this->beam1.SetN(np,nm);
        this->beam1.SetSize(sx,sy,sz,sdelta);
    }
    else{
        this->beam2.SetN(np,nm);
        this->beam2.SetSize(sx,sy,sz,sdelta);
    }
}

void beambeam::set_beam(int beamid, double energy, double mass, double charge, int zslice) {
    if (beamid==1){
        this->beam1.SetEnergy(energy,mass,charge);
        this->beam1.SetSlice(zslice);
    }
    if (beamid==2){
        this->beam2.SetEnergy(energy,mass,charge);
        this->beam2.SetSlice(zslice);
    }
}

void beambeam::set_crab(int beamid, double freq, double deviation, int harm, double ratio) {
    if (beamid==1){
        this->cc1b=CCrabCavity(1, freq,this->crossing_angle, deviation, harm, ratio);
        this->cc1a=CCrabCavity(-1, freq,this->crossing_angle, deviation, harm, ratio);
    }
    else{
        this->cc2b=CCrabCavity(1, freq,this->crossing_angle, deviation, harm, ratio);
        this->cc2a=CCrabCavity(-1, freq,this->crossing_angle, deviation, harm, ratio);
    }
}

void beambeam::set_crab_optics(int beamid, double beta, double alpha, double eta, double etap, double phase_adv ){
    if (beamid==11) {
        this->cc1b.set_optics(beta,alpha,eta,etap,phase_adv,this->map1x.beta, this->map1x.alpha, 0.0, 0.0);
    }
    else if (beamid==12) {
        this->cc1a.set_optics(beta, alpha, eta, etap, phase_adv, this->map1x.beta, this->map1x.alpha, 0.0, 0.0);
    }
    else if (beamid==21) {
        this->cc2b.set_optics(beta, alpha, eta, etap, phase_adv, this->map2x.beta, this->map2x.alpha, 0.0, 0.0);
    }
    else if (beamid==22){
        this->cc2a.set_optics(beta,alpha,eta,etap,phase_adv,this->map2x.beta, this->map2x.alpha, 0.0, 0.0);
    }
    
}


void beambeam::set_IPmap(int beamid, double b, double a, double phi, double chrom) {
    if (beamid==11){
        this->map1x=COneTurnMap(b,a,phi,chrom);
    }
    else if (beamid==12){
        this->map1y=COneTurnMap(b,a,phi,chrom);
    }
    else if (beamid==21){
        this->map2x=COneTurnMap(b,a,phi,chrom);
    }
    else if (beamid==22){
        this->map2y=COneTurnMap(b,a,phi,chrom);
    }
}
void beambeam::set_rf(int beamid, double voltage, double freq, double phase, int h, double rt, double tunez) {
    if (beamid==1){
        this->rf1=Crf(rt,voltage,phase,h,freq);
        this->rf1.force_tune(tunez);
    }
    else{
        this->rf2=Crf(rt,voltage,phase,h,freq);
        this->rf2.force_tune(tunez);
    }
}



std::vector<double> beambeam::get_distribution(const int &ind) {
    if (ind==10) return beam1.x_;
    else if (ind==11) return beam1.px_;
    else if (ind==12) return beam1.y_;
    else if (ind==13) return beam1.py_;
    else if (ind==14) return beam1.z_;
    else if (ind==15) return beam1.delta_;
    else if (ind==20) return beam2.x_;
    else if (ind==21) return beam2.px_;
    else if (ind==22) return beam2.y_;
    else if (ind==23) return beam2.py_;
    else if (ind==24) return beam2.z_;
    else if (ind==25) return beam2.delta_;
    
    
}

std::vector<double> beambeam::get_statistics(){
    
        return std::vector<double>{beam1.xmean,beam1.pxmean,beam1.xsize,beam1.pxsize,beam1.xemit,
                                   beam1.ymean,beam1.pymean,beam1.ysize,beam1.pysize,beam1.yemit,
                                   beam2.xmean,beam2.pxmean,beam2.xsize,beam2.pxsize,beam2.xemit,
                                   beam2.ymean,beam2.pymean,beam2.ysize,beam2.pysize,beam2.yemit,};
}

std::vector<double> beambeam::get_optics(){
    return std::vector<double>{beam1.zmax, beam1.zmin, beam2.zmax, beam2.zmin};
}


std::vector<double> beambeam::beambeam_interaction_2pass() {
    this->crab_crossing(1, 1);
    this->crab_crossing(2, 1);
    std::vector<double> slice_output;
    slice_output=beambeam_2pass(beam1, beam2, 1);

    this->crab_crossing(1, -1);
    this->crab_crossing(2, -1);
    return slice_output;
}

void beambeam::crab_crossing(int beamid, int forward){
    if (beamid==1){
        if (forward>0) {
            this->beam1.crab_kick(cc1b);
            this->beam1.LorentzTransform(this->crossing_angle, forward);
        }
        else{
            this->beam1.LorentzTransform(this->crossing_angle, forward);
            this->beam1.crab_kick(cc1a);
        }
    }
    else {
        if (forward>0) {
            this->beam2.crab_kick(cc2b);
            this->beam2.LorentzTransform(this->crossing_angle, forward);
        }
        else{
            this->beam2.LorentzTransform(this->crossing_angle, forward);
            this->beam2.crab_kick(cc2a);
        }
    }
}

