//
// Created by Yue Hao on 7/29/16.
//

#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include "mpi.h"
#include "beam.h"
#include "acc_model.h"
#include "omp.h"
#ifndef EPIC_BEAMBEAM_H
#define EPIC_BEAMBEAM_H

const double sqrtpi=sqrt(M_PI);


class beambeam{
private:
    int model;
    std::mt19937 gen;
public:
    
    MPI_Comm mcomm;
    double crossing_angle;
    beam beam1;
    beam beam2;
    CCrabCavity cc1b, cc1a;
    CCrabCavity cc2b, cc2a;
    COneTurnMap map1x, map1y;
    COneTurnMap map2x, map2y;
    Crf rf1;
    Crf rf2;

    beambeam():gen(std::random_device()()),crossing_angle(0),model(1){}
    beambeam(const double& crossing_angle, const int& model=0):crossing_angle(crossing_angle), model(model),gen(std::random_device()()){
        omp_set_num_threads(4);
        #pragma omp parallel
        if (omp_get_thread_num() == 1)
            std::cout << omp_get_num_threads()<<std::endl;
    }
    
    //void set_mpi(MPI_Comm mcomm){this->mcomm=mcomm;}
    void set_crab(int beamid, double freq, double deviation=1.0, int harm=3, double ratio=0.0);
    void set_crab_optics(int beamid, double beta, double alpha, double eta, double etap, double phase_adv );
    void set_IPmap(int beamid, double b, double a, double phi, double chrom);
    void set_rf(int beamid, double voltage, double freq, double phase, int h, double rt, double tunez);
    void set_beam_param(int beamid, double np, unsigned long nm, double sx, double sy, double sz, double sdelta);
    void set_beam(int beamid, double energy, double mass, double charge, int zlice);
    void initialize(){ beam1.Ini6D(map1x, map1y, rf1, gen); beam2.Ini6D(map2x, map2y, rf2, gen);}
    void cal_emit(){beam1.CalEmit(); beam2.CalEmit();}
    std::vector<double> get_distribution(const int& ind);
    std::vector<double> get_statistics();
    std::vector<double> get_optics();
    std::vector<double> beambeam_interaction_2pass();
    void oneturn(){beam1.OneTurn(map1x, map1y, rf1); beam1.OneTurn(map2x, map2y, rf2);}
    void crab_crossing(int beamid, int forward);
};




#endif //EPIC_BEAMBEAM_H
