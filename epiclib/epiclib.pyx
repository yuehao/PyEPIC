cimport cython
from mpi4py import MPI


from libcpp.vector cimport vector
from libcpp.string cimport string
import numpy as np



cdef extern from "../../../EPIC/beambeam.h":
    cdef cppclass beambeam:
        beambeam() except +
        beambeam(const double& crossing_angle, const int& model=0) except +


        void set_crab(int beamid, double freq, double deviation, int harm, double ratio);
        void set_crab_optics(int beamid, double beta, double alpha, double eta, double etap, double phase_adv);
        void set_IPmap(int beamid, double b, double a, double phi, double chrom);
        void set_rf(int beamid, double voltage, double freq, double phase, int h, double rt, double tunez);

        void set_beam_param(int beamid, double np, unsigned long nm, double sx, double sy, double sz, double sdelta);
        void set_beam(int beamid, double energy, double mass, double charge, int zslice);
        void crab_crossing(int beamid, int forward);
        void initialize()
        void cal_emit()
        vector[double] get_distribution(const int& ind);
        vector[double] get_statistics();
        vector[double] get_optics();
        vector[double] beambeam_interaction_2pass();




cdef class epiclib:
    cdef:
        beambeam bb
    def __cinit__(self, float crossing_angle):
        self.bb = beambeam(crossing_angle)


    def initialize(self):
        self.bb.initialize()


    def set_crab(self, int beamid, float freq, float deviation, int harm, double ratio, *optics):
        self.bb.set_crab(beamid, freq, deviation, harm, ratio)
        if len(optics)>0:
            self.bb.set_crab_optics(beamid*10+1, optics[0], optics[1], optics[2], optics[3], optics[4])
            self.bb.set_crab_optics(beamid*10+2, optics[5], optics[6], optics[7], optics[8], optics[9])



    def set_IPmap(self, int beamid, float b, float a, float phi, float chrom):
        self.bb.set_IPmap(beamid, b, a, phi, chrom)

    def set_rf(self, int beamid, float voltage, float freq, float phase, int harm, float gammat, float tunez):
        self.bb.set_rf(beamid, voltage, freq, phase, harm, gammat, tunez)

    def set_beam(self, int beamid, float num_part, long num_macro, float energy, float mass, float charge, int zslice, float xsize, float ysize, float zsize, float deltaE):
        self.bb.set_beam(beamid, energy, mass, charge, zslice)
        self.bb.set_beam_param(beamid, num_part, num_macro, xsize, ysize, zsize, deltaE)

    def cal_emit(self):
        self.bb.cal_emit()

    def get_distribution(self):
        temp10 = self.bb.get_distribution(10)
        temp11 = self.bb.get_distribution(11)
        temp12 = self.bb.get_distribution(12)
        temp13 = self.bb.get_distribution(13)
        temp14 = self.bb.get_distribution(14)
        temp15 = self.bb.get_distribution(15)
        #beam1=np.vstack([temp0,temp1,temp2,temp3,temp4,temp5])
        temp20 = self.bb.get_distribution(20)
        temp21 = self.bb.get_distribution(21)
        temp22 = self.bb.get_distribution(22)
        temp23 = self.bb.get_distribution(23)
        temp24 = self.bb.get_distribution(24)
        temp25 = self.bb.get_distribution(25)
        #beam2=np.vstack([temp0,temp1,temp2,temp3,temp4,temp5])
        return (np.array(temp10), np.array(temp11), np.array(temp12), np.array(temp13), np.array(temp14), np.array(temp15),
                np.array(temp20), np.array(temp21), np.array(temp22), np.array(temp23), np.array(temp24), np.array(temp25))

    def get_statistics(self):
        return self.bb.get_statistics()

    def get_optics(self):
        return self.bb.get_optics()

    def beambeam_interaction_2pass(self):
        return np.array(self.bb.beambeam_interaction_2pass()).reshape((9,-1))


    def crab_crossing(self, int beamid, int forward):
        return self.bb.crab_crossing(beamid, forward)



