from .beam import Beam
import numpy as np
from . import epicIO
#from .epiclib import epiclib as epiclib
import scipy.constants as constants
import scipy.integrate as integrate


class PyEPIC(object):
    def __init__(self, inputfile=None):
        self.spline = 1
        self.crossing_angle = 0.0
        self.tracking_turns = 0
        self.beam1 = Beam()
        self.beam2 = Beam()
        if inputfile is not None:
            self.readfrom(inputfile)
        self.beam1.param.crossing_angle = self.crossing_angle
        self.beam2.param.corssing_angle = self.crossing_angle





    def readfrom(self, inputfile):
        with open(inputfile,'r') as fs:
            glbdict,b1dict,b2dict = epicIO.parse_file(fs)
        self.beam1.read_input(b1dict)
        self.beam2.read_input(b2dict)
        self.set_param(**glbdict)




    def set_param(self, **param_dict):

        namelist = self.__dict__.keys()
        for k, v in param_dict.items():
            if k in namelist:
                self.__dict__[k] = v
            else:
                print("Warning, {} is not a recognized key, ignored the action {}={}.".format(k, k, v))


        self._calculate()


    def get_param(self, param_name):
        namelist = self.__dict__.keys()
        if param_name.lower() in namelist:
            return self.__dict__[param_name.lower()]
        else:
            print("{} is not a recognized key.".format(param_name))

    def _calculate(self):
        self.beam1.param.crossing_angle = self.crossing_angle
        self.beam2.param.crossing_angle = self.crossing_angle

    def weak_strong_beambeam(self, strong, weak, from_distribution=False):
        '''Beam 2 is always weak.'''
        #self.lumi, self.weak_strong_result = self.beam2.weak_strong_beam_beam(self.beam1, even_z=True, make_movie=make_evolving_movie, )
        #return self.lumi,self.beam2.emittance[0],self.beam2.beamsize[0],self.beam2.geo_beta[0],self.beam2.geo_alpha[0], self.beam2.effective_rms_emittance[0],self.beam2.effective_max_emittance[0], self.beam2.max_emittance[0], np.max(self.weak_strong_result[:,7])

    #def get_distribution(self):
    #    (self.beam1._x, self.beam1._px, self.beam1._y, self.beam1._py, self.beam1._ct, self.beam1._dE,
    #     self.beam2._x, self.beam2._px, self.beam2._y, self.beam2._py, self.beam2._ct, self.beam2._dE)= self.epic_cpp.get_distribution()


    def weak_strong_plot_evolution(self):
        import matplotlib.pyplot as plt

        f1, ax1 = plt.subplots()
        p1=ax1.scatter(self.weak_strong_result[:,0], self.weak_strong_result[:,1]*1e6, color='b', label='Geometric emittance')
        p2=ax1.scatter(self.weak_strong_result[:,0], self.weak_strong_result[:,3]*1e6, color='r', label='Effective emittance')

        ax1.set_xlabel('s [m]')
        ax1.set_ylabel('rms horizontal emittance [mm mrad]')
        ax2=ax1.twinx()
        p3=ax2.scatter(self.weak_strong_result[:,0], self.weak_strong_result[:,2]*1e6, color='g', label='Beam size')
        p4=ax2.plot(self.weak_strong_result[:,0], np.sqrt(1+np.power((self.weak_strong_result[:,0]-
                       self.beam1.param.s_star[0])/self.beam1.param.beta_star[0],2.0))*np.sqrt(self.beam1.param.emittance[0]*self.beam1.param.beta_star[0])*1e6,
                    linestyle='None', c='0.85')
        ax1.set_ylim(bottom=0, top=0.02)
        ax2.set_ylabel('rms horizontal beam size [micron]')
        ax2.set_ylim(bottom=0, top=100)
        plt.legend((p1,p2,p3),(p1.get_label(),p2.get_label(),p3.get_label()))
        #plt.show()
        plt.savefig("evolution_x.eps")



        f1, ax1 = plt.subplots()
        p1 = ax1.scatter(self.weak_strong_result[:, 0], self.weak_strong_result[:, 4] * 1e6, color='b',
                         label='Geometric emittance')
        p2 = ax1.scatter(self.weak_strong_result[:, 0], self.weak_strong_result[:, 6] * 1e6, color='r',
                         label='Effective emittance')

        ax1.set_xlabel('s [m]')
        ax1.set_ylabel('rms vertical emittance [mm mrad]')
        ax2 = ax1.twinx()
        p3 = ax2.scatter(self.weak_strong_result[:, 0], self.weak_strong_result[:, 5] * 1e6, color='g',
                         label='Beam size')
        p4 = ax2.plot(self.weak_strong_result[:, 0], np.sqrt(1 + np.power((self.weak_strong_result[:, 0] -
                                                                              self.beam1.param.s_star[1]) /
                                                                             self.beam1.param.beta_star[1],
                                                                             2.0)) * np.sqrt(
            self.beam1.param.emittance[1] * self.beam1.param.beta_star[1]) * 1e6,linestyle='None',marker='o' ,c='0.85')
        ax1.set_ylim(bottom=0, top=0.02)
        ax2.set_ylabel('rms vertical beam size [micron]')
        ax2.set_ylim(bottom=0, top=100)
        plt.legend((p1, p2, p3), (p1.get_label(), p2.get_label(), p3.get_label()))
        # plt.show()
        plt.savefig("evolution_y.eps")

    def weak_strong_plot_radiation(self):
        '''
        import pylab
        pylab.clf()
        ax1 = pylab.subplot()
        pylab.hist(self.beam2.eloss/1000*self.beam2.param.energy,bins=100)
        ax1.set_xlabel("Radiation energy [KeV]")
        ax1.set_ylabel("Histogram")
        pylab.show()
        pylab.clf()
        ax1 = pylab.subplot()
        pylab.hist(self.beam2.ncl,bins=100)
        ax1.set_xlabel("Photon number")
        ax1.set_ylabel("Histogram")
        pylab.show()
        pylab.clf()
        ax1 = pylab.subplot()
        pylab.hist(self.beam2.critical_energy/1000*self.beam2.param.energy,bins=100)
        ax1.set_xlabel("Critical photon energy [KeV]")
        ax1.set_ylabel("Histogram")
        pylab.show()
        #print(np.average(eloss)*self.param.energy,np.max(eloss)*self.param.energy)
        '''


    def beambeam_parameters(self):
        pa1 = self.beam1.param.rms_beamsize[-1] / self.beam1.param.rms_beamsize[0] * self.crossing_angle
        pa2 = self.beam2.param.rms_beamsize[-1] / self.beam2.param.rms_beamsize[0] * self.crossing_angle
        print("Piwinsky Angle, beam_1: {:.2f} \t beam_2: {:.2f}".format(pa1,pa2))

        print("")

        cc_product = self.beam1.param.charge * self.beam2.param.charge
        cc_sign=cc_product / np.abs(cc_product)
        xi1_ratio = 2.0 * np.abs(
            cc_product) * self.beam2.param.n_particle * self.beam1.param.classical_radius / self.beam1.param.gamma
        sum_beamsize = self.beam2.param.rms_beamsize[0] + self.beam2.param.rms_beamsize[1]
        xix = self.beam1.param.beta_IP[0] * xi1_ratio / sum_beamsize / self.beam2.param.rms_beamsize[0] / 4 / np.pi
        xiy = self.beam1.param.beta_IP[1] * xi1_ratio / sum_beamsize / self.beam2.param.rms_beamsize[1] / 4 / np.pi
        print("Beam-Beam parameter for beam1, Horizontal:{:.2e}, Vertical:{:.2e}".format(xix, xiy))
        self.beambeam_parameter_1 = xix, xiy
        cos_phi = np.cos(2 * np.pi * self.beam1.param.tune[0:2])
        sin_phi = np.sin(2 * np.pi * self.beam1.param.tune[0:2])
        trace = 2 * cos_phi + sin_phi * 4 * np.pi * np.array(self.beambeam_parameter_1) * cc_sign
        newphase = np.arccos(trace / 2.0)
        dtune = newphase / 2 / np.pi - np.array(self.beam1.param.tune[0:2])
        sin_newphase = np.sin(newphase)
        print("Tune change for beam1, Horizontal:{:.2e}, Vertical:{:.2e}".format(dtune[0], dtune[1]))
        new_beta = np.array(self.beam1.param.beta_IP) * np.sin(2 * np.pi * self.beam1.param.tune[0:2]) / sin_newphase
        new_alpha = (2 * self.beam1.param.alpha_IP + 4 * np.pi * np.array(
            self.beambeam_parameter_1) * cc_sign) * np.sin(2 * np.pi * self.beam1.param.tune[0:2]) / 2 / sin_newphase
        print("Dynamic optics for beam1, Horizontal:({:.2f},{:.2f}), Vertical:({:.2f},{:.2f})".format(new_beta[0],
                                                                                                      new_alpha[0],
                                                                                                      new_beta[1],
                                                                                                      new_alpha[1]))
        print("")


        xi2_ratio = 2.0 * np.abs(cc_product) * self.beam1.param.n_particle * self.beam2.param.classical_radius / self.beam2.param.gamma
        sum_beamsize = self.beam1.param.rms_beamsize[0] + self.beam1.param.rms_beamsize[1]
        xix = self.beam2.param.beta_IP[0] * xi2_ratio / sum_beamsize / self.beam1.param.rms_beamsize[0] / 4 / np.pi
        xiy = self.beam2.param.beta_IP[1] * xi2_ratio / sum_beamsize / self.beam1.param.rms_beamsize[1] / 4 / np.pi
        print("Beam-Beam parameter for beam2, Horizontal:{:.2e}, Vertical:{:.2e}".format(xix, xiy))
        self.beambeam_parameter_2 = xix, xiy
        cos_phi = np.cos(2 * np.pi * self.beam2.param.tune[0:2])
        sin_phi = np.sin(2 * np.pi * self.beam2.param.tune[0:2])
        trace = 2 * cos_phi + sin_phi * 4 * np.pi * np.array(self.beambeam_parameter_2) * cc_sign
        newphase = np.arccos(trace / 2.0)
        dtune = newphase / 2 / np.pi - np.array(self.beam2.param.tune[0:2])
        sin_newphase = np.sin(newphase)
        print("Tune change for beam2, Horizontal:{:.2e}, Vertical:{:.2e}".format(dtune[0], dtune[1]))
        new_beta = np.array(self.beam2.param.beta_IP) * np.sin(2 * np.pi * self.beam2.param.tune[0:2]) / sin_newphase
        new_alpha = (2 * self.beam2.param.alpha_IP + 4 * np.pi * np.array(
            self.beambeam_parameter_2) * cc_sign) * np.sin(2 * np.pi * self.beam2.param.tune[0:2]) / 2 / sin_newphase
        print("Dynamic optics for beam2, Horizontal:({:.2f},{:.2f}), Vertical:({:.2f},{:.2f})".format(new_beta[0],
                                                                                                      new_alpha[0],
                                                                                                      new_beta[1],
                                                                                                      new_alpha[1]))


    def sample_inputfile(self, filename):
        with open(filename, 'w') as f:
            f.write('''
                    &beam1
                    n_particle=3e11
                    n_macro=100000
                    energy=275e9
                    species=proton
                    charge=1
                    rms_energy_spread=1e-4
                    emittance=3.41e-9,3.41e-9,0
                    rms_beamsize=0,0,0.16
                    slice_z=101
                    beta_star=0.26,0.13
                    s_star= 0.0, 0.0
                    tune=0.675,0.685,0.0043
                    chrom=0.0,0.0
                    damping_decrement=0.0,0.0,0.0
                    &end
                    
                    &beam2
                    n_particle=0.33e11
                    n_macro=100000
                    energy=12.0e9
                    species=electron
                    charge=-1
                    rms_energy_spread=1e-4
                    emittance=3.41e-9,3.41e-9,0
                    rms_beamsize=0,0,0.004
                    slice_z=1
                    beta_star=0.26,0.13
                    offset=0e-6,0,0,0,0,0
                    s_star=0.0, 0.0
                    tune=0.675,0.685,0.0043
                    chrom=0.0,0.0
                    damping_decrement=0.0,0.0,0.0
                    &end
                    
                    &global
                    spline=3
                    &end
                    ''')

    def lumi_integration(self):
        csa= np.cos(self.crossing_angle)
        sna= np.sin(self.crossing_angle)
        zrange = 15 * (self.beam1.param.rms_beamsize[-1] + self.beam2.param.rms_beamsize[-1])

        def sx1sqr(s):
            return self.beam1.param.rms_beamsize[0] ** 2 * (1 + s*s / self.beam1.param.beta_IP[0]**2)

        def sx2sqr(s):
            return self.beam2.param.rms_beamsize[0] ** 2 * (1 + s*s / self.beam2.param.beta_IP[0]**2)


        def t0(s, ct):
            return (- (-ct + csa * s) ** 2 / (2 * self.beam1.param.rms_beamsize[-1] ** 2)
                    - (ct + csa * s ) ** 2 / (2 * self.beam2.param.rms_beamsize[-1] ** 2)
                    - (-s * sna + self.beam1.param.crabbing_scale * sna * (self.beam1.crab_deviation(-ct+s))) ** 2 / (2 * sx1sqr(s))
                    - (s * sna - self.beam2.param.crabbing_scale * sna * (self.beam2.crab_deviation(ct+s))) ** 2 / (2 * sx2sqr(s)))

        def t1(s, ct):
            return (-(((-ct + csa * s) * sna) / self.beam1.param.rms_beamsize[-1] ** 2)
                    + ((ct + csa * s) * sna) /self.beam2.param.rms_beamsize[-1] ** 2
                    - (csa * (-s*sna + self.beam1.param.crabbing_scale * sna * self.beam1.crab_deviation(-ct+s))) / sx1sqr(s)
                    - (csa * (s * sna - self.beam2.param.crabbing_scale * sna * self.beam2.crab_deviation(ct+s))) / sx2sqr(s)
                    )

        def t2(s, ct):
            return (1.0/ 2*(csa**2 * (-(1 / sx1sqr(s)) - 1 / sx2sqr(s)) +
                            sna ** 2 * (-(1 / self.beam1.param.rms_beamsize[-1] ** 2)
                                        - 1 / self.beam1.param.rms_beamsize[-1] ** 2 )))

        def lumi(s, ct):
            return (np.exp(t0(s, ct) - t1(s, ct)**2/(4 * t2(s, ct)))
                    /(np.sqrt(-t2(s,ct)) * np.sqrt(sx1sqr(s) * sx2sqr(s)))
                    /self.beam1.param.rms_beamsize[-1]/self.beam2.param.rms_beamsize[-1]
                    )


        return integrate.dblquad(lumi, -zrange, zrange, -zrange, zrange)


    def kink_instability_matrix(self, mnp1, mnp2):
        dim = (mnp1+mnp2)*2
        mat = np.zeros((dim,dim))
        def phi_ij (us, Nump, i, j):
            return (us - 2 * np.pi * (Nump - i + j) / Nump) / 2
        us_temp = self.beam1.param.tune[-1] * 2 * np.pi
        imn1 = np.array([np.sin(mnp1*phi_ij(us_temp, mnp1, i ,j))/mnp1/phi_ij(us_temp, mnp1, i ,j) for j in range(mnp1) for i in range(mnp1)])
        imn1.reshape((mnp1,mnp1))
        return imn1

    def beambeam_interaction_2pass(self):
        return self.epic_cpp.beambeam_interaction_2pass()

    def to_beambeam3d(self, ncpu=(32,2), grids=(128,128,2048), ):

        with open("beam1.in", 'w') as f1, open("beam2.in", 'w') as f2:
            np.set_printoptions(precision=16)
            # number of cpus
            f1.write('{:d}\t{:d}\n'.format(*ncpu))
            f2.write('{:d}\t{:d}\n'.format(*ncpu))

            # dim_phase_space and num macro particle
            f1.write('{:d}\t{:d}\n'.format(6, int(self.beam1.param.n_macro)))
            f2.write('{:d}\t{:d}\n'.format(6, int(self.beam2.param.n_macro)))

            # grids
            f1.write('{:d}\t{:d}\t{:d}\n'.format(*grids))
            f2.write('{:d}\t{:d}\t{:d}\n'.format(*grids))

            #longitudinal slices
            f1.write('{:d}\n'.format(int(self.beam1.param.slice_z)))
            f2.write('{:d}\n'.format(int(self.beam1.param.slice_z)))

            # nturns
            f1.write('{:d}\n'.format(int(self.tracking_turns)))
            f2.write('{:d}\n'.format(int(self.tracking_turns)))

            # initial distributions, 2 for Gaussian
            f1.write('2\n')
            f2.write('2\n')

            # distributions parameters x : beamsize, beta, alpha, x-sc, px-sc, cn.x, cn.px
            f1.write('{:e}\t{:f}\t{:f}\t0.0\t0.0\t0.0\t0.0\n'.format(self.beam1.param.rms_beamsize[0],
                                                                        self.beam1.param.beta_IP[0],
                                                                        self.beam1.param.alpha_IP[0],
                                                                        ))
            f2.write('{:e}\t{:f}\t{:f}\t0.0\t0.0\t0.0\t0.0\n'.format(self.beam2.param.rms_beamsize[0],
                                                                        self.beam2.param.beta_IP[0],
                                                                        self.beam2.param.alpha_IP[0],
                                                                        ))

            # distributions parameters y : beamsize, beta, alpha, y-sc, py-sc, cn.y, cn.py
            f1.write('{:e}\t{:f}\t{:f}\t0.0\t0.0\t0.0\t0.0\n'.format(self.beam1.param.rms_beamsize[1],
                                                                        self.beam1.param.beta_IP[1],
                                                                        self.beam1.param.alpha_IP[1],
                                                                        ))
            f2.write('{:e}\t{:f}\t{:f}\t0.0\t0.0\t0.0\t0.0\n'.format(self.beam2.param.rms_beamsize[1],
                                                                        self.beam2.param.beta_IP[1],
                                                                        self.beam2.param.alpha_IP[1],
                                                                        ))

            # distributions parameters z : bunch_length, energy_spread, al, c,c,s,s
            f1.write('{:f}\t{:e}\t0.0\t1.0\t1.0\t0.0\t0.0\n'.format(self.beam1.param.rms_beamsize[2],
                                                                     self.beam1.param.rms_energy_spread
                                                                     ))
            f2.write('{:f}\t{:e}\t0.0\t1.0\t1.0\t0.0\t0.0\n'.format(self.beam2.param.rms_beamsize[2],
                                                                     self.beam2.param.rms_energy_spread
                                                                     ))

            # num_particle_in_bunch, beam_energy, beam_rest_mass, beam_charge
            f1.write('{:e}\t{:e}\t{:e}\t{:f}\n'.format(self.beam1.param.n_particle,
                                                       self.beam1.param.energy,
                                                       self.beam1.param.mass,
                                                       self.beam1.param.charge
                                                       )
                     )
            f2.write('{:e}\t{:e}\t{:e}\t{:f}\n'.format(self.beam2.param.n_particle,
                                                       self.beam2.param.energy,
                                                       self.beam2.param.mass,
                                                       self.beam2.param.charge
                                                       )
                     )

            #tunes: x,y,z
            f1.write('{:f}\t{:f}\t{:f}\n'.format(*self.beam1.param.tune))
            f2.write('{:f}\t{:f}\t{:f}\n'.format(*self.beam2.param.tune))

            #optical alpha: x,y
            f1.write('{:f}\t{:f}\n'.format(*(self.beam1.param.alpha_IP[0:2])))
            f2.write('{:f}\t{:f}\n'.format(*(self.beam2.param.alpha_IP[0:2])))

            # optical beta: x,y
            f1.write('{:f}\t{:f}\n'.format(*(self.beam1.param.beta_IP[0:2])))
            f2.write('{:f}\t{:f}\n'.format(*(self.beam2.param.beta_IP[0:2])))

            # transverse emittance: x,y
            f1.write('{:e}\t{:e}\n'.format(*(self.beam1.param.emittance[0:2])))
            f2.write('{:e}\t{:e}\n'.format(*(self.beam2.param.emittance[0:2])))

            # damping turns in: x , y , z, switch
            switch=0
            if np.all(self.beam1.param.damping_turns):
                switch=1
            f1.write('{:e}\t{:e}\t{:e}\t{:d}\n'.format(*(self.beam1.param.damping_turns), switch))
            if np.all(self.beam2.param.damping_turns):
                switch=1
            f2.write('{:e}\t{:e}\t{:e}\t{:d}\n'.format(*(self.beam2.param.damping_turns), switch))

            # longitudinal: sigmaz, energy spread
            f1.write('{:e}\t{:e}\n'.format(self.beam1.param.rms_beamsize[2],self.beam1.param.rms_energy_spread))
            f2.write('{:e}\t{:e}\n'.format(self.beam2.param.rms_beamsize[2],self.beam2.param.rms_energy_spread))

            # longitudinal: sigmaz, energy spread
            f1.write('{:e}\t{:e}\n'.format(self.beam1.param.rms_beamsize[2], self.beam1.param.rms_energy_spread))
            f2.write('{:e}\t{:e}\n'.format(self.beam2.param.rms_beamsize[2], self.beam2.param.rms_energy_spread))

            # sweeping: switch, ??
            f1.write('{:d}\t{:f}\t{:f}\n'.format(0, 0.0, 0.01))
            f2.write('{:d}\t{:f}\t{:f}\n'.format(0, 0.0, 0.01))

            # orbit feed back: switch, ????
            f1.write('{:d}\t{:f}\t{:f}\t{:f}\t{:f}\n'.format(0, 0.0, 0.0, 6.28, 6.28))
            f2.write('{:d}\t{:f}\t{:f}\t{:f}\t{:f}\n'.format(0, 0.0, 0.0, 6.28, 6.28))

            # closed orbit squeeze: switch, ?
            f1.write('{:d}\t{:d}\n'.format(0, 20))
            f2.write('{:d}\t{:d}\n'.format(0, 20))

            # closed orbit x/y coordinate: x, y
            f1.write('{:e}\t{:e}\n'.format(0.0, 0.0))
            f2.write('{:e}\t{:e}\n'.format(0.0, 0.0))

            # longitudinal again: sigmaz, energy spread
            f1.write('{:e}\t{:e}\n'.format(self.beam1.param.rms_beamsize[2], self.beam1.param.rms_energy_spread))
            f2.write('{:e}\t{:e}\n'.format(self.beam2.param.rms_beamsize[2], self.beam2.param.rms_energy_spread))

            pass





'''
        
        
        
        
        
        0.0
        d0
        12.5
        d - 3 / cross
        ang
        alpha(x in x - y), phi(s
        s - x)
        0
        0
        0
        0
        0 / fix.cmp.dm., var.slc., ext.mp, ext
        .2
        nd
        mp, ext
        .3
        rd
        map, ext.
        4
        th
        ma
        p
        1 /  # of collisions per turn
        -5.0
        5.0 - 5.0
        5.0
        2.0
        2.0 / comp.dom: xmin, xmax, ymin, ymax, zmin, zmax
        2
        1 / PE
        group, type
        of
        s - w
        interaction
        1.
        d0
        1.
        d0 / x and y
        linear
        chromaticity
        0.
        d0 / curvature
        of
        ring

'''