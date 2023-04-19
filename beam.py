import numpy as np
import scipy as sp
import scipy.constants as const
import scipy.optimize as optimize
from . import listfunc as lf
from . import element
import copy
import multiprocessing as mp


class BeamParameter(object):
    def __init__(self,):

        self.n_particle = 0  # number of particle in bunch
        self.n_macro = 0
        self.energy = 1e9  # eV
        self.species = 'electron'
        self.charge = -1
        
        self.rms_energy_spread = 0
        self.rms_emittance = np.array([0.0,0.0,0.0])  # m-rad, mrad, ev-s
        self.rms_beamsize = np.array([0.0,0.0,0.0])  # sigma_x,sigma_y, sigma_z

        self.offset = np.zeros(6) #imperfection of [x,px,y,py,z,dE]
        self.damping_decrement = np.array([0.0, 0.0, 0.0])
        self.chrom=np.array([0.0,0.0])

        self.slice_z = 0
        self.gaussian_cutoff = 5
        self.transverse_dist = 'Gaussian'
        self.longitudinal_dist = 'Gaussian'

        self.beta_star = np.array([0, 0])
        self.beta_IP = np.array([1, 1])
        self.alpha_IP = np.array([0, 0])
        self.s_star = np.array([0, 0])

        self.beam_beta = None
        self.beam_alpha = None

        self.tune = np.array([0, 0, 0])


        self.rf_voltage = 0
        self.rf_harm = 1
        self.rf_freq = 28.15e6
        self.rf_phase = 0
        self.momentum_compaction = 1/23.57/23.57

        #self.crossing_angle=5e-3   # half of the crossing angle
        self.crab_cavity_freq=197e6
        self.crabbing_scale = 1.0
        self.crab_cavity_beta=[1000.0,1000.0]
        self.crab_cavity_alpha = [0, 0]
        self.crab_cavity_eta = [0, 0]
        self.crab_cavity_etap = [0, 0]
        self.crab_cavity_phaseadv_IP = [np.pi/2.0, -np.pi/2.0]

        self.crab_cavity_harmonic_number=np.array([3, ])
        self.crab_cavity_harmonic_ratio=np.array([0.0, ])

        self.mass = 0
        self.classical_radius = 0
        self.gamma = 0
        self.beta = 0

        self.crossing_angle=0  #half angle



    def _calculate(self):
        #beam species determines mass and r and gamma beta
        self.species = self.species.upper()
        if self.species == 'ELECTRON':
            self.mass = const.physical_constants['electron mass energy equivalent in MeV'][0]*1e6
            self.classical_radius = 2.8179403267e-15
            self.cgamma = 8.846e-5
        elif self.species == 'PROTON':
            self.mass = const.physical_constants['proton mass energy equivalent in MeV'][0]*1e6
            self.classical_radius = 2.8179403267e-15*const.physical_constants['electron mass energy equivalent in MeV'][0]/const.physical_constants['proton mass energy equivalent in MeV'][0]
            self.cgamma = 7.783e-18
        self.gamma = self.energy/self.mass
        self.beta = np.sqrt(1-1/self.gamma/self.gamma)
        self.pc = np.sqrt(self.energy * self.energy - self.mass * self.mass)

        #crab cavity
        self.crab_cavity_harmonic_freq= self.crab_cavity_freq * self.crab_cavity_harmonic_number

        #s* shift optics function
        
        if self.beta_star.any():
            self.alpha_IP = self.s_star/self.beta_star
            self.beta_IP = self.beta_star+self.s_star*self.s_star/self.beta_star
            self.gamma_IP = (1+self.alpha_IP*self.alpha_IP)/self.beta_IP

        else:
            self.gamma_IP = (1+self.alpha_IP*self.alpha_IP)/self.beta_IP
            self.beta_star = 1/self.gamma_IP
            self.s_star = self.beta_star*self.alpha_IP

        if self.beam_beta is None:
            self.beam_beta = self.beta_IP*1.0
            self.beam_alpha = self.alpha_IP*1.0


        if self.rms_beamsize[:2].any():
            self.rms_emittance[:2] = self.rms_beamsize[:2]*self.rms_beamsize[:2]/self.beta_IP
        else:
            self.rms_beamsize[:2] = np.sqrt(self.rms_emittance[:2] * self.beta_IP)

        self.initial_beamsize = copy.deepcopy(self.rms_beamsize)
        self.initial_energy_spread = self.rms_energy_spread





        self.IPmap=np.eye(6,6)
        phase = self.tune[0] * 2 * np.pi
        beta = self.beta_IP[0]
        alpha = self.alpha_IP[0]
        gamma = self.gamma_IP[0]
        self.IPmap[0, 0] = np.cos(phase) + alpha * np.sin(phase)
        self.IPmap[0, 1] = beta * np.sin(phase)
        self.IPmap[1, 0] = -gamma * np.sin(phase)
        self.IPmap[1, 1] = np.cos(phase) - alpha * np.sin(phase)

        phase = self.tune[1] * 2 * np.pi
        beta = self.beta_IP[1]
        alpha = self.alpha_IP[1]
        gamma = self.gamma_IP[1]
        self.IPmap[2, 2] = np.cos(phase) + alpha * np.sin(phase)
        self.IPmap[2, 3] = beta * np.sin(phase)
        self.IPmap[3, 2] = -gamma * np.sin(phase)
        self.IPmap[3, 3] = np.cos(phase) - alpha * np.sin(phase)

        phase = self.tune[2] * 2 * np.pi
        beta = self.rms_beamsize[-1]/self.rms_energy_spread
        alpha = 0
        gamma = 1.0/beta
        self.IPmap[4, 4] = np.cos(phase) + alpha * np.sin(phase)
        self.IPmap[4, 5] = beta * np.sin(phase)
        self.IPmap[5, 4] = -gamma * np.sin(phase)
        self.IPmap[5, 5] = np.cos(phase) - alpha * np.sin(phase)

        if self.rf_voltage != 0:
            self.rf_eta = self.momentum_compaction-1.0/self.gamma**2
            etacosphi=self.rf_eta*np.cos(self.rf_phase)

            if etacosphi>0:
                self.rf_phase=np.pi-self.rf_phase
                etacosphi*=-1

            self.rf_unstable_phase = np.pi - self.rf_phase
            tune = np.sqrt(-self.rf_voltage*self.rf_harm*etacosphi/2/np.pi/self.beta**2/self.energy)
            ratio = self.rf_freq*2*np.pi/const.c*tune/np.abs(self.rf_eta)/self.rf_harm
            self.tune[-1]=tune
            self.rms_energy_spread=ratio*self.rms_beamsize[-1]


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



class Beam(object):
    def __init__(self, nproc=1):
        self.param = BeamParameter()
        self._x = np.empty(0)
        self._px = np.empty(0)
        self._y = np.empty(0)
        self._py = np.empty(0)
        self._ct = np.empty(0)
        self._dE = np.empty(0)
        self._ndim = 0
        self.nproc=nproc

    def read_input(self, input_dict):
        from multiprocessing.sharedctypes import RawArray
        self.param.set_param(**input_dict)

        self.param.n_macro=(self.param.n_macro//self.nproc)*self.nproc
        nmp = int(self.param.n_macro)

        self._XRAW = RawArray('d', nmp)
        self._PXRAW = RawArray('d', nmp)
        self._YRAW = RawArray('d', nmp)
        self._PYRAW = RawArray('d', nmp)
        self._CTRAW = RawArray('d', nmp)
        self._DERAW = RawArray('d', nmp)

        self._x=np.frombuffer(self._XRAW)
        self._px=np.frombuffer(self._PXRAW)
        self._y = np.frombuffer(self._YRAW)
        self._py = np.frombuffer(self._PYRAW)
        self._ct = np.frombuffer(self._CTRAW)
        self._dE = np.frombuffer(self._DERAW)


    def get_optics(self, loc):
        beta=self.param.beta_star+np.power(loc-self.param.s_star,2.0)/self.param.beta_star
        alpha=(loc-self.param.s_star)/self.param.beta_star
        return beta[0], alpha[0], beta[1], alpha[1]

    def get_action(self, direction='x', optics=(1,0)):
        if direction == 'x':
            return (self._x * self._x + np.power(optics[1] * self._x + optics[0] * self._px, 2.0)) / optics[0] / 2.0
        elif direction == 'y':
            return (self._y * self._y + np.power(optics[1] * self._y + optics[0] * self._py, 2.0)) / optics[0] / 2.0

    def get_emittance(self):

        self.average = np.array([np.average(self._x), np.average(self._y), np.average(self._ct)])
        self.average_p = np.array([np.average(self._px), np.average(self._py), np.average(self._dE)])
        self.beamsize = np.array([np.std(self._x), np.std(self._y), np.std(self._ct)])
        self.beamsize_p = np.array([np.std(self._px), np.std(self._py), np.std(self._dE)])
        self.correlate = np.array([np.correlate(self._x, self._px)[0], np.correlate(self._y, self._py)[0],
                                   np.correlate(self._ct, self._dE)[0]]) / len(self._x) - self.average * self.average_p
        self.emittance = np.sqrt(
            self.beamsize * self.beamsize * self.beamsize_p * self.beamsize_p - self.correlate * self.correlate)
        self.beam_beta = self.beamsize * self.beamsize / self.emittance
        self.beam_alpha = -self.correlate / self.emittance

    def get_effect_emittance(self,loc=0, get_max=True):
        optics_at_loc = self.get_optics(loc)
        self.action_x = self.get_action(direction='x', optics=(optics_at_loc[0], optics_at_loc[1]))
        self.action_y = self.get_action(direction='y', optics=(optics_at_loc[2], optics_at_loc[3]))
        self.effective_rms_emittance = (np.average(self.action_x), np.average(self.action_y))
        self.effective_max_emittance = (np.max(self.action_x), np.max(self.action_y))
        if get_max:
            geo_action_x = self.get_action(direction='x', optics=(self.beam_beta[0], self.beam_alpha[0]))
            geo_action_y = self.get_action(direction='y', optics=(self.beam_beta[1], self.beam_alpha[1]))
            self.max_emittance = np.array([np.max(geo_action_x), np.max(geo_action_y)])
        return optics_at_loc

    def fold(self, ratio=1):
        self._x += self._ct * self._px * ratio
        self._y += self._ct * self._py * ratio

    def one_turn_4D(self):
        self._x, self._px, self._y, self._py=self.param.IPmap[0:4,0:4].dot(np.array([self._x, self._px, self._y, self._py]))

    def one_turn_x(self):
        self._x, self._px=self.param.IPmap[0:2,0:2].dot(np.array([self._x, self._px]))
    
    def one_turn_x_chromatic(self, chrom):
        dphi = self._dE*chrom*2*np.pi
        beta = self.param.beta_IP[0]
        alpha = self.param.alpha_IP[0]
        gamma = self.param.gamma_IP[0]
        cosdphi=np.cos(dphi)
        sindphi=np.sin(dphi)
        m11 = cosdphi + alpha * sindphi
        m12 = beta * sindphi
        m21 = -gamma * sindphi
        m22 = cosdphi - alpha * sindphi
        self._x, self._px = m11*self._x+m12*self._px , m21*self._x+m22*self._px

        

    def one_turn_y(self):
        self._y, self._py=self.param.IPmap[2:4,2:4].dot(np.array([self._y, self._py]))

    def one_turn_y_chromatic(self, chrom):
        dphi = self._dE*chrom*2*np.pi
        beta = self.param.beta_IP[1]
        alpha = self.param.alpha_IP[1]
        gamma = self.param.gamma_IP[1]
        cosdphi=np.cos(dphi)
        sindphi=np.sin(dphi)
        m11 = cosdphi + alpha * sindphi
        m12 = beta * sindphi
        m21 = -gamma * sindphi
        m22 = cosdphi - alpha * sindphi
        self._y, self._py = m11*self._y+m12*self._py , m21*self._y+m22*self._py

    def one_turn_z(self):
        self._ct, self._dE=self.param.IPmap[4:,4:].dot(np.array([self._ct, self._dE]))

    def one_turn_z_RF(self):
        ztphi=2*np.pi*self.param.rf_freq/const.c
        sinphi=np.sin(self.param.rf_phase)
        c1=self.param.rf_voltage/self.param.energy/self.param.beta**2
        c2=2*np.pi*self.param.rf_harm*self.param.rf_eta/ztphi

        self._dE+=c1*(np.sin(ztphi*self._ct+self.param.rf_phase)-sinphi)
        self._ct+=c2*self._dE

    def plot_distribution(self, axis, xy=('x','px'), units=('[mm]', '[mrad]'), scale=(1.0, 1.0), bins=400, range=None, zrange=None,
                          draw_emit_ellipse=False, threshold=0.001, location=0, alpha=0.9, cmap=None, norm=None):
        from matplotlib import colors
        from matplotlib.patches import Ellipse
        from matplotlib.colors import LogNorm

        if location != 0:
            #self.drift(location)
            self._ct += location
        h, xedge, yedge, = np.histogram2d(self.__dict__['_' + xy[1]]*scale[1], self.__dict__['_' + xy[0]]*scale[0],
                                          bins=bins, range=range, normed=False)
        np.ma.masked_where(h>=threshold, h, copy=False)

        if zrange is None:
            maxz=np.max(h)
            minz=np.min(h)
            print(maxz, minz)
        else:
            maxz = zrange[0]
            minz = zrange[1]


        if cmap is None:
            cmap = colors.LinearSegmentedColormap.from_list('mycolor',
                                                        [(0, 'white'),
                                                         (1e-4, 'magenta'),
                                                         (0.1, 'blue'),
                                                         (0.2, 'cyan'),
                                                         (0.4, 'green'),
                                                         (0.6, 'yellow'),
                                                         (0.8, 'orange'),
                                                         (1.0, 'red')])
        else:
            cmap.set_bad(color='white')
            h=np.ma.masked_where(h <= threshold, h, copy=False)

        image = axis.imshow(h, interpolation='nearest', origin='lower', extent=[yedge[0], yedge[-1], xedge[0], xedge[-1]],
                            aspect='auto', cmap=cmap, vmin=minz, vmax=maxz,alpha=alpha, norm=norm)
        # axis.scatter(self.__dict__['_'+xy[0]]*1.0e3,self.__dict__['_'+xy[1]]*1.0e3)
        axis.set_xlabel(xy[0] + ' ' + units[0])
        axis.set_ylabel(xy[1] + ' ' + units[1])
        #axis.set_xlim(xrange)
        #axis.set_ylim(yrange)

        if draw_emit_ellipse:
            optics_at_loc = self.get_effect_emittance(loc=0, get_max=True)
            ang = np.arange(500) / (500 / 2 / np.pi)
            emitlist = [self.emittance[0], self.max_emittance[0], self.effective_rms_emittance[0],
                        self.effective_max_emittance[0]]
            betalist = [self.beam_beta[0], self.beam_beta[0], optics_at_loc[0], optics_at_loc[0]]
            alphalist = [self.beam_alpha[0], self.beam_alpha[0], optics_at_loc[1], optics_at_loc[1]]
            colorlist = ['green', 'yellow', 'red', 'blue']
            legend=['rms geometric emittance', '100% geometric emittance', 'rms effective emittance', '100% effective emittance']
            for i in range(len(emitlist)):
                r = np.sqrt(2 * betalist[i] * emitlist[i] * 1.0e6)
                x = r * np.cos(ang)
                xp = (-r * np.sin(ang) - alphalist[i] * x) / betalist[i]
                axis.plot(x, xp, color=colorlist[i], label=legend[i])
            axis.legend(loc=1)

        '''if draw_position:
            extra_height = yrange[1] / 3
            ep = Ellipse((zshift * xrange[1] / 4.0, yrange[1] + extra_height / 2.0), width=xrange[1],
                         height=extra_height / 2.0, facecolor='Blue')
            ee = Ellipse((-zshift * xrange[1] / 4.0, yrange[1] + extra_height / 2.0), width=xrange[1] / 20.0,
                         height=extra_height / 2.0, facecolor='red')
            axis.set_ylim([yrange[0], yrange[1] + extra_height])
            ep.set_alpha(0.5)
            ee.set_alpha(0.5)
            axis.add_patch(ep)
            axis.add_patch(ee)
            '''
        if location != 0:
            #self.drift(-location)
            self._ct -= location
        return image

    def synchrotron_radiation(self, cur=0, l=0, bx=0, by=0):
        u = 1 + self._dE
        rl = (self._x * cur + 1) * l
        absn = np.sqrt(np.power(1 + self._x * cur, 2.0) + self._px * self._px / u / u + self._py * self._py / u / u)
        ex = self._px / u / absn
        ey = self._py / u / absn
        ez = (1 + self._x * cur) / absn
        b2 = by * ez * by * ez + bx * ez * bx * ez + np.power(bx * ey - by * ex, 2.0)
        # omegac = 3*np.power(self.param.gamma,3.0)*np.sqrt(b2)*3e8*6.52e-16
        # print(np.max(omegac))

        eloss = self.param.cgamma / 2.0 / np.pi * np.power(self.param.energy / 1.0e9, 3.0) * u * u * b2 * rl
        return eloss

    def put_error(self):
        self._x += self.param.offset[0]
        self._px += self.param.offset[1]
        self._y += self.param.offset[2]
        self._py += self.param.offset[3]
        self._ct += self.param.offset[4]
        self._dE += self.param.offset[5]

    def SR_damping_excitation(self):
        if self.param.damping_decrement.any():
            size= len(self._x)
            exp_damp= np.exp(-self.param.damping_decrement)
            sqrt_exp_damp=np.sqrt(1.0-exp_damp*exp_damp)
            self._x *= exp_damp[0]
            self._px *= exp_damp[0]
            self._y *= exp_damp[1]
            self._py *= exp_damp[1]
            self._ct *= exp_damp[2]
            self._dE *= exp_damp[2]

            self._x += np.random.randn(size) * (sqrt_exp_damp[0] * self.param.rms_beamsize[0])
            self._px += np.random.randn(size) * (sqrt_exp_damp[0] * self.param.rms_beamsize[0]/self.param.beta_IP[0])

            self._y += np.random.randn(size) * (sqrt_exp_damp[1] * self.param.rms_beamsize[1])
            self._py += np.random.randn(size) * (sqrt_exp_damp[1] * self.param.rms_beamsize[1]/self.param.beta_IP[1])

            self._ct += np.random.randn(size) * (sqrt_exp_damp[2] * self.param.rms_beamsize[2])
            self._dE += np.random.randn(size) * (sqrt_exp_damp[2] * self.param.rms_energy_spread)

    def zslicing(self, from_distribution=False, even_z=True):



        import scipy.special as spefunc

        self.nslice = int(self.param.slice_z) / 2 * 2 + 1  # Always make it an odd number
        halfn = (self.nslice - 1) / 2

        nump = self.param.n_particle
        size_z = self.param.rms_beamsize[2]

        cutoff = self.param.gaussian_cutoff
        if even_z:
            self.zlist = (np.arange(self.nslice) - halfn) * size_z ** 2.0 / self.nslice  # where pos is located
            self.zsep = np.concatenate(
                (np.array([-cutoff * size_z]), self.zlist + cutoff * size_z / self.nslice))  # where momentum is located
            # print(zsep,zlist)

        else:
            portion_in_cutoff = spefunc.erf(cutoff / np.sqrt(2.0))
            partition = (np.arange(halfn) + 1) / (halfn + 0.5) * portion_in_cutoff
            halfzlist = spefunc.erfinv(partition) * np.sqrt(2.0)
            partition = (np.arange(halfn + 1) + 0.5) / (halfn + 0.5) * portion_in_cutoff
            halfzsep = spefunc.erfinv(partition) * np.sqrt(2.0)
            self.zlist = np.concatenate((-halfzlist[::-1], np.array([0]), halfzlist)) * size_z
            self.zsep = np.concatenate((-halfzsep[::-1], halfzsep)) * size_z


        if from_distribution:
            self.np_in_zslices = np.histogram(self._ct, bins=self.zsep)
        else:
            self.np_in_zslices = np.diff(nump * (0.5 + 0.5 * spefunc.erf(self.zsep / np.sqrt(2.0) / size_z)))

        np_sum = np.sum(self.np_in_zslices)
        self.np_in_zslices *= (nump / np_sum)






    def bb_field_2d(self, nbin=1000, method='GAUSS', plt_axis = None): # normalized kick=E*r0/gamma default is Gauss law, later can extend to PIC solver
        if method == 'GAUSS':
            r = np.sqrt(self._x * self._x + self._y * self._y)
            r2 = np.average(r * r) / 2.0
            # r2=self.param.rms_beamsize[0]*self.param.rms_beamsize[0]
            nc, rc = np.histogram(r, bins=nbin)

            efield = 2.0 * np.cumsum(nc) / rc[1:] / self.param.n_macro * self.param.n_particle
            efield2 = 2.0 * self.param.n_particle * (1 - np.exp(-rc[1:] * rc[1:] / 2.0 / r2)) / rc[1:]


            plt_axis.plot(rc[1:] * 1e6, efield, label='From beam distribution')
            plt_axis.plot(rc[1:] * 1e6, efield2, label='1 Gaussian Fit')
            plt_axis.set_xlabel("Radius [micron]")
            plt_axis.set_ylabel("Beam-beam field [Arb. Unit]")

            return rc[1:]



    def transverse_moment(self, max_order=3):
        self.x_moment = np.zeros(max_order + 1)
        self.y_moment = np.zeros(max_order + 1)
        self.r_moment = np.zeros(max_order + 1)
        average_x = 0
        average_y = 0
        std_x = 0
        std_y = 0

        for i in range(max_order + 1):
            if i == 0:
                self.x_moment[0] = 1
                self.y_moment[0] = 1
                self.r_moment[0] = 1
            elif i == 1:
                average_x = np.average(self._x)
                self.x_moment[i] = average_x
                average_y = np.average(self._y)
                self.y_moment[i] = average_y
                self.r_moment[i] = 0
            elif i == 2:
                std_x = np.std(self._x)
                std_y = np.std(self._y)
                self.x_moment[i] = std_x
                self.y_moment[i] = std_y
                self.r_moment[i] = np.sqrt(
                    np.average(np.power(self._x - average_x, 2.0) + np.power(self._y - average_y, 2.0)) / 2.0)
            else:
                self.x_moment[i] = np.average(np.power(self._x - average_x, 1.0 * i)) / np.power(std_x, 1.0 * i)
                self.y_moment[i] = np.average(np.power(self._y - average_y, 1.0 * i)) / np.power(std_y, 1.0 * i)
                self.r_moment[i] = np.average(
                    np.power((np.power(self._x - average_x, 2.0) + np.power(self._x - average_x, 2.0)) / 2.0,
                             0.5 * i)) / np.power(self.r_moment[2], 1.0 * i)


    def lorentz_boost(self):
        '''
        :param angle: the crossing angle in x direction
        :return: Change the coordinates in funciton
        '''
        angle = self.param.crossing_angle
        cos_ang = np.cos(angle)
        sin_ang = np.sin(angle)
        tan_ang = sin_ang / cos_ang

        ps_list = np.sqrt((1.0 + self._dE) * (1.0 + self._dE) - self._px * self._px - self._py * self._py)
        h_list = 1.0 + self._dE - ps_list

        self._py = self._py / cos_ang
        h_list = h_list / (cos_ang*cos_ang)
        self._px = self._px / cos_ang - h_list * sin_ang;
        self._dE = self._dE - self._px * sin_ang

        ps_list = 1 + self._dE - h_list
        ds = self._x * sin_ang
        self._x += tan_ang * self._ct + ds * self._px / ps_list
        self._y += ds * self._py / ps_list
        self._ct = self._ct / cos_ang - ds * h_list / ps_list



    def inv_lorentz_boost(self):
        angle = self.param.crossing_angle
        cos_ang = np.cos(angle)
        sin_ang = np.sin(angle)
        tan_ang = sin_ang / cos_ang

        ps_list = np.sqrt((1.0 + self._dE) * (1.0 + self._dE) - self._px * self._px - self._py * self._py)
        h_list = 1.0 + self._dE - ps_list

        self._x=(self._x-self._ct*sin_ang)/(1+(self._px+h_list*sin_ang)*sin_ang/ps_list)
        self._y-=self._py*self._x*sin_ang/ps_list
        self._ct=(self._ct+h_list*self._x*sin_ang/ps_list)*cos_ang


        self._py*=cos_ang
        self._dE+=self._px*sin_ang
        self._px=(self._px+h_list*sin_ang)*cos_ang

    def crab_kick(self, direction=1, scale=1):
        angle=self.param.crossing_angle
        krf = 2 * np.pi * self.param.crab_cavity_freq / const.c

        rsum=np.sum(self.param.crab_cavity_harmonic_ratio)

        dx = -(1 - rsum) * angle / krf * np.sin(krf * self._ct)
        de = - (1-rsum) * angle * np.cos(krf*self._ct) * self._px

        if rsum==0.0:
            self._x += dx * direction * scale
            self._dE += de * direction * scale
            return

        for m,r in zip(self.param.crab_cavity_harmonic_number, self.param.crab_cavity_harmonic_ratio):
            dx -= r * angle / krf / m * np.sin(m * krf * self._ct)
            de -= r * angle * np.sin(m * krf * self._ct) * self._px

        self._x += dx * direction * scale
        self._dE += de * direction * scale


    def crab_kick_ver(self, direction=1, scale=0):
        angle=self.param.crossing_angle
        krf = 2 * np.pi * self.param.crab_cavity_freq / const.c

        rsum=np.sum(self.param.crab_cavity_harmonic_ratio)

        dx = -(1 - rsum) * angle / krf * np.sin(krf * self._ct)
        de = - (1-rsum) * angle * np.cos(krf*self._ct) * self._px

        if rsum==0.0:
            self._y += dx * direction * scale
            self._dE += de * direction * scale
            return

        for m,r in zip(self.param.crab_cavity_harmonic_number, self.param.crab_cavity_harmonic_ratio):
            dx -= r * angle / krf / m * np.sin(m * krf * self._ct)
            de -= r * angle * np.sin(m * krf * self._ct) * self._px

        self._y += dx * direction * scale
        self._dE += de * direction * scale

    def solenoid_lab(self, bfield, lstart, lend):
        g = bfield/2.0/self.param.pc * const.speed_of_light
        theta= g*(lend-lstart)
        print(theta)

        element.drift_2D(self._x, self._px, lstart)
        element.drift_2D(self._y, self._py, lstart)
        element.thick_quad_2d(self._x, self._px, g * g, lend-lstart)
        element.thick_quad_2d(self._y, self._py, g * g, lend-lstart)
        element.rotation(self._x, self._y, theta)
        element.rotation(self._px, self._py, theta)
        element.drift_2D(self._x, self._px, -lend)
        element.drift_2D(self._y, self._py, -lend)

    def xy_rotation(self, theta):
        element.rotation(self._x, self._y, theta)
        element.rotation(self._px, self._py, theta)


    def crab_deviation(self, z):
        krf = 2 * np.pi * self.param.crab_cavity_freq / const.c
        rsum = np.sum(self.param.crab_cavity_harmonic_ratio)


        dz=(self.param.crabbing_scale-rsum)*np.sin(krf*z)/krf
        if rsum==0.0:
            return (z-dz)*self.param.crossing_angle
        else:
            for m, r in zip(self.param.crab_cavity_harmonic_number, self.param.crab_cavity_harmonic_ratio):
                dz+=r / krf / m * np.sin(m * krf * z)*self.param.crabbing_scale
            return (z-dz)*self.param.crossing_angle

    def one_pass_linear_BB_map(self, oppo_beam, even_z=True, kick_IP=None, initial_position=None):
        import scipy.special as spefunc
        oppo_np = oppo_beam.param.n_particle
        nslice = int(oppo_beam.param.slice_z / 2) * 2 + 1  # Always make it an odd number
        halfn = (nslice - 1) / 2
        cutoff = oppo_beam.param.gaussian_cutoff
        oppo_size_z=oppo_beam.param.rms_beamsize[-1]
        total_mat = np.eye(3)
        if kick_IP is not None:
            total_mat[0, 2] = kick_IP[0]
            total_mat[1, 2] = kick_IP[1]
            if initial_position is not None:
                initial_position[0] += kick_IP[0]
                initial_position[1] += kick_IP[1]

        if even_z:
            zlist = (np.arange(nslice)-halfn)*oppo_size_z*cutoff*2.0/nslice   #where pos is located
            zsep = np.concatenate((np.array([-cutoff*oppo_size_z]),zlist+cutoff*oppo_size_z/nslice)) # where momentum is located
            #print(zsep,zlist)
            oppo_np_in_z=np.diff(oppo_np*(0.5+0.5*spefunc.erf(zsep/np.sqrt(2.0)/oppo_size_z)))
        else:
            portion_in_cutoff = spefunc.erf(cutoff/np.sqrt(2.0))
            partition = (np.arange(halfn)+1)/(halfn+0.5)*portion_in_cutoff

            halfzlist = spefunc.erfinv(partition)*np.sqrt(2.0)
            #print(halfzlist)
            partition = (np.arange(halfn+1)+0.5)/(halfn+0.5)*portion_in_cutoff

            halfzsep = spefunc.erfinv(partition)*np.sqrt(2.0)
            zlist = np.concatenate((-halfzlist[::-1],np.array([0]),halfzlist))*oppo_size_z
            zsep = np.concatenate((-halfzsep[::-1],halfzsep))*oppo_size_z

            oppo_np_in_z = np.diff(oppo_np*(0.5+0.5*spefunc.erf(zsep/np.sqrt(2.0)/oppo_size_z)))
        oppo_np_sum = np.sum(oppo_np_in_z)

        if oppo_np_sum > 0:
            oppo_np_in_z *= (oppo_np/oppo_np_sum)

        if initial_position is not None:
            temp_x=np.zeros(nslice)
            temp_px=np.zeros(nslice)


        for i in reversed(range(nslice)):
            this_mat=np.eye(3)
            zpos = -(zlist[i])/2.0
            #lf.drift_2D(self._x, self._px, zpos)
            #lf.drift_2D(self._y, self._py, zpos)
            if initial_position is not None:
                initial_position[0] += initial_position[1] * zpos
                temp_x[i] = initial_position[0]
                temp_px[i] = initial_position[1]

            this_mat[0,1] = zpos
            total_mat = np.dot(this_mat , total_mat)

            sigma_x = np.sqrt(1 + np.power((zpos - oppo_beam.param.s_star[0]) / oppo_beam.param.beta_star[0], 2.0)) * oppo_beam.param.rms_beamsize[0]
            sigma_y = np.sqrt(1 + np.power((zpos - oppo_beam.param.s_star[1]) / oppo_beam.param.beta_star[1], 2.0)) * oppo_beam.param.rms_beamsize[1]

            sigma_2x2my2 = np.sqrt(2*(sigma_x*sigma_x-sigma_y*sigma_y))
            ccdev = -oppo_beam.crab_deviation(zlist[i])

            expterm=np.exp(-ccdev*ccdev/2/sigma_x/sigma_x)
            argx1 = ccdev/sigma_2x2my2
            w1r = np.exp(-argx1 * argx1)
            w1i = 2.0 * spefunc.dawsn(argx1) / np.sqrt(np.pi)

            argx2 = ccdev*sigma_y/sigma_x / sigma_2x2my2
            w2r = np.exp(-argx2 * argx2)
            w2i = 2.0 * spefunc.dawsn(argx2) / np.sqrt(np.pi)


            factor = 2 * self.param.charge * oppo_beam.param.charge * oppo_np_in_z[i] * self.param.classical_radius * np.sqrt(np.pi) / self.param.gamma / sigma_2x2my2
            #invf = invf / (sigma_x+sigma_y)/sigma_x/self.param.gamma

            zero_order = (w1i - expterm * w2i) * factor
            first_order = 2/np.sqrt(np.pi)/sigma_2x2my2-expterm*2/np.sqrt(np.pi)*sigma_y/sigma_x/sigma_2x2my2
            first_order -= 2* argx1 * w1i / sigma_2x2my2
            first_order += (ccdev/sigma_x/sigma_x + 2*argx2*sigma_y/sigma_x/sigma_2x2my2) *expterm * w2i
            first_order *= factor
            this_mat = np.eye(3)

            this_mat[1,0]= first_order

            this_mat[1,2]= zero_order
            total_mat = np.dot(this_mat, total_mat)

            if initial_position is not None:
                #initial_position[1] += (initial_position[0]-ccdev) * invf
                #initial_position[1] += zero_order
                initial_position[1] += zero_order + first_order * initial_position[0]



            this_mat = np.eye(3)
            # lf.drift_2D(self._x, self._px, -zpos)
            # lf.drift_2D(self._y, self._py, -zpos)
            this_mat[0, 1] = -zpos
            total_mat = np.dot(this_mat, total_mat)
            if initial_position is not None:
                initial_position[0] -= initial_position[1] * zpos


        if kick_IP is not None:
            this_mat = np.eye(3)

            this_mat[0, 2] = kick_IP[0]
            this_mat[1, 2] = -kick_IP[1]
            total_mat = np.dot(this_mat, total_mat)
            if initial_position is not None:
                initial_position[0] += kick_IP[0]
                initial_position[1] -= kick_IP[1]

        oneturn = np.eye(3)
        oneturn[0:2, 0:2] = self.param.IPmap[0:2, 0:2]
        oneturn = np.dot(oneturn, total_mat)

        if initial_position is not None:
            return oneturn, (zlist, temp_x, temp_px)

        closed_orbit_IP=np.dot(np.linalg.inv(np.eye(2)-oneturn[0:2,0:2]), oneturn[0:2, 2])

        #print(np.dot(np.linalg.inv(np.eye(2)-oneturn[0:2,0:2]),oneturn[0:2,0:2]-self.param.IPmap[0:2,0:2]))
        return self.one_pass_linear_BB_map(oppo_beam, even_z, kick_IP, closed_orbit_IP)



    def strong_beam_setup(self, even_z=True):
        import scipy.special as spefunc
        oppo_np = self.param.n_particle

        size_z = self.param.rms_beamsize[-1]
        #oppo_s_star_x, oppo_s_star_y = self.param.s_star
        #oppo_beta_star_x, oppo_beta_star_y = self.param.beta_star
        self.ss_nslice = int(self.param.slice_z) // 2 * 2 + 1  # Always make it an odd number
        halfn = (self.ss_nslice - 1) / 2
        cutoff = self.param.gaussian_cutoff
        if even_z:
            self.ss_zlist = (np.arange(self.ss_nslice)-halfn)*size_z*cutoff*2.0/self.ss_nslice   #where pos is located
            self.ss_zsep = np.concatenate((np.array([-cutoff*size_z]),self.ss_zlist+cutoff*size_z/self.ss_nslice)) # where momentum is located
            #print(zsep,zlist)
            self.np_in_z=np.diff(oppo_np*(0.5+0.5*spefunc.erf(self.ss_zsep/np.sqrt(2.0)/size_z)))
        else:
            portion_in_cutoff = spefunc.erf(cutoff/np.sqrt(2.0))
            partition = (np.arange(halfn)+1)/(halfn+0.5)*portion_in_cutoff

            halfzlist = spefunc.erfinv(partition)*np.sqrt(2.0)
            #print(halfzlist)
            partition = (np.arange(halfn+1)+0.5)/(halfn+0.5)*portion_in_cutoff

            halfzsep = spefunc.erfinv(partition)*np.sqrt(2.0)
            self.ss_zlist = np.concatenate((-halfzlist[::-1],np.array([0]),halfzlist))*size_z
            self.ss_zsep = np.concatenate((-halfzsep[::-1],halfzsep))*size_z

            self.np_in_z = np.diff(oppo_np * (0.5 + 0.5 * spefunc.erf(self.ss_zsep / np.sqrt(2.0) / size_z)))
        oppo_np_sum = np.sum(self.np_in_z)
        #print(len(zsep))
        delta_z = np.diff(self.ss_zlist)
        if oppo_np_sum > 0:
            self.np_in_z *= (oppo_np / oppo_np_sum)
        self.ss_crabdev=self.crab_deviation(self.ss_zlist)





    def weak_strong_beam_beam(self, oppo_beam, freq=7.827e4, nbunch=1190, solenoid_str=None, solenoid_tilting_angle=0):


        import scipy.special as spefunc

        oppo_size_x, oppo_size_y, oppo_size_z = oppo_beam.param.rms_beamsize
        oppo_s_star_x, oppo_s_star_y = oppo_beam.param.s_star
        oppo_beta_star_x, oppo_beta_star_y = oppo_beam.param.beta_star

        wsresult=np.empty(0)

        lumi = 0
        charge = 1.0*self.param.charge*oppo_beam.param.charge

        sol_g=0
        bbrho = 0
        if solenoid_str is not None:

            sol_g = solenoid_str / 2.0 / self.param.pc * const.speed_of_light
            self._px+=solenoid_tilting_angle
            zedge1 = (self._ct + oppo_beam.ss_zsep[0]) / 2.0
            sol_alpha = sol_g * zedge1
            element.rotation(self._x, self._y, sol_alpha)
            element.rotation(self._px, self._py, sol_alpha)
            self._px -= solenoid_tilting_angle


        for i in range(oppo_beam.ss_nslice):
            #this_mat=np.eye(3)
            #if i==0:
            self._x -= oppo_beam.ss_crabdev[i]
            zpos = (self._ct + oppo_beam.ss_zlist[i])/2.0
            zedge1= (self._ct + oppo_beam.ss_zsep[i])/2.0
            zedge2 = (self._ct + oppo_beam.ss_zsep[i+1]) / 2.0
            dz=(oppo_beam.ss_zsep[i+1]-oppo_beam.ss_zsep[i])/2.0

            element.drift_2D(self._x, self._px, zpos)
            element.drift_2D(self._y, self._py, zpos)
            #this_mat[0,1] = zpos


            sigma_x = np.sqrt(1+np.power((zpos-oppo_s_star_x)/oppo_beta_star_x,2.0))*oppo_size_x
            sigma_y = np.sqrt(1+np.power((zpos-oppo_s_star_y)/oppo_beta_star_y,2.0))*oppo_size_y
            #print(sigma_x.shape, sigma_y.shape)
            lumi+=np.sum(oppo_beam.np_in_z[i] / 2.0 / np.pi / sigma_x / sigma_y * np.exp(-np.power(self._x / sigma_x, 2.0) / 2.0 - np.power(self._y / sigma_y, 2.0) / 2.0))

            sigma_x2_y2 = 2.0*(sigma_x*sigma_x-sigma_y*sigma_y)

            r2 = self._x*self._x+self._y*self._y
            epsilon=(oppo_size_x+oppo_size_y)/2.0e2*0
            #sigma_flag_0 = np.logical_and(sigma_x2_y2<epsilon, sigma_x2_y2>-epsilon )

            sigma_flag_pos = (sigma_x2_y2 > epsilon )
            sigma_flag_neg = (sigma_x2_y2 < -epsilon)
            dpx = np.zeros_like(sigma_x)
            dpy = np.zeros_like(sigma_y)
            #dp2 = oppo_np_in_z[i]*self.param.classical_radius/self.param.gamma/r2[sigma_flag_0]*(1-np.exp(-r2[sigma_flag_0]/2.0/np.power(sigma_x[sigma_flag_0],2.0)))

            #dpx[sigma_flag_0] = 2.0*charge*dp2*self._x[sigma_flag_0]
            #dpy[sigma_flag_0] = 2.0*charge*dp2*self._y[sigma_flag_0]
            w1pos = sp.special.wofz((np.abs(self._x[sigma_flag_pos])+1.0j*np.abs(self._y[sigma_flag_pos]))/np.sqrt(sigma_x2_y2[sigma_flag_pos]))

            w2pos = sp.special.wofz((np.abs(self._x[sigma_flag_pos])*sigma_y[sigma_flag_pos]/sigma_x[sigma_flag_pos]
                                      + 1.0j*np.abs(self._y[sigma_flag_pos])*sigma_x[sigma_flag_pos]/sigma_y[sigma_flag_pos])/np.sqrt(sigma_x2_y2[sigma_flag_pos]))

            dp2 = -2.0j*oppo_beam.np_in_z[i]*self.param.classical_radius/self.param.gamma*np.sqrt(np.pi/sigma_x2_y2[sigma_flag_pos])*\
                                (w1pos-np.exp(-np.power(self._x[sigma_flag_pos]/sigma_x[sigma_flag_pos],2.0)/2.0-np.power(self._y[sigma_flag_pos]/sigma_y[sigma_flag_pos],2.0)/2.0)*w2pos)
            dpx[sigma_flag_pos] = charge*dp2.real*np.sign(self._x[sigma_flag_pos])
            dpy[sigma_flag_pos] = -charge*dp2.imag*np.sign(self._y[sigma_flag_pos])

            w1pos = sp.special.wofz((np.abs(self._y[sigma_flag_neg])+1.0j*np.abs(self._x[sigma_flag_neg]))/np.sqrt(-sigma_x2_y2[sigma_flag_neg]))
            w2pos = sp.special.wofz((np.abs(self._y[sigma_flag_neg])*sigma_x[sigma_flag_neg]/sigma_y[sigma_flag_neg]
                                      + 1.0j*np.abs(self._x[sigma_flag_neg])*sigma_y[sigma_flag_neg]/sigma_x[sigma_flag_neg])/np.sqrt(-sigma_x2_y2[sigma_flag_neg]))

            dp2 = -2.0j*oppo_beam.np_in_z[i]*self.param.classical_radius/self.param.gamma*np.sqrt(-np.pi/sigma_x2_y2[sigma_flag_neg])*\
                                (w1pos-np.exp(-np.power(self._x[sigma_flag_neg]/sigma_x[sigma_flag_neg],2.0)/2.0
                                              -np.power(self._y[sigma_flag_neg]/sigma_y[sigma_flag_neg],2.0)/2.0)*w2pos)
            dpy[sigma_flag_neg] = charge*dp2.real*np.sign(self._y[sigma_flag_neg])
            dpx[sigma_flag_neg] = -charge*dp2.imag*np.sign(self._x[sigma_flag_neg])
            self._px += dpx
            self._py += dpy
            #dpxy = np.sqrt(dpx*dpx+dpy*dpy)
            #dz = (zsep[i+1]-zsep[i])/2.0
            #self.eloss+=self.synchrotron_radiation(l=dz,bx=dpx/dz,by=dpy/dz)
            #self.ncl += 2.5/np.sqrt(3.0)*7.297e-3*self.param.gamma*dpxy
            #self.critical_energy=15*np.sqrt(3.0)/8*self.eloss/self.ncl
            #dpola += 5*7*9/4.0/np.sqrt(3.0)/54.0*7.297e-3*2.426e-12*2.426e-12*np.power(self.param.gamma,5.0)*np.power(dpxy,3.0)/np.power(dz,2.0)
            #self.get_emittance()
            #wsresult = np.append(wsresult, [zlist[i]/2.0, self.emittance[0],self.beamsize[0],self.effective_rms_emittance[0],
            #                                self.emittance[1], self.beamsize[1], self.effective_rms_emittance[1], np.max(self._px), np.max(self._py)])
            #if i < nslice-1:
            #    lf.drift_2D(self._x,self._px, delta_z[i]/2.0)
            #    lf.drift_2D(self._y,self._py, delta_z[i]/2.0)
            #    zpos += delta_z[i]/2.0

            #else:
            element.drift_2D(self._x, self._px, -zpos)
            element.drift_2D(self._y, self._py, -zpos)
            self._x += oppo_beam.ss_crabdev[i]
            if solenoid_str is not None:
                sol_alpha=dz*sol_g
                self._px += solenoid_tilting_angle
                element.rotation(self._x, self._y, sol_alpha)
                element.rotation(self._px, self._py, sol_alpha)
                self._px -= solenoid_tilting_angle




            #if make_movie:
            #    f1,a1=plt.subplots()
            #    image=self.plot_distribution(a1, zshift=zlist[i]/zlist[-1], draw_ellipse=True, draw_position=True)
            #    plt.colorbar(image)
            #    plt.savefig('fig{:03d}.png'.format(i))
            #    plt.clf()

        if solenoid_str is not None:
            zedge2 = (self._ct + oppo_beam.ss_zsep[-1]) / 2.0
            sol_alpha = -sol_g * zedge2
            self._px += solenoid_tilting_angle
            element.rotation(self._x, self._y, sol_alpha)
            element.rotation(self._px, self._py, sol_alpha)
            self._px -= solenoid_tilting_angle



        self.get_emittance()

        #wsresult = wsresult.reshape((oppo_beam.ss_nslice, -1))
        lumi = lumi*self.param.n_particle/self.param.n_macro*freq*nbunch/1.0e4
        return lumi
   
    def mirror_copy_distribution(self, dimlim=6):
        if dimlim <= 2:
            self._x = np.append(self._x, -self._x)
            self._y = np.append(self._y, self._y)
            self._px = np.append(self._px, self._px)
            self._py = np.append(self._py, self._py)
            self._ct = np.append(self._ct, self._ct)
            self._dE = np.append(self._dE, self._dE)
            self.param.n_macro = self.param.n_macro*2

            self._x = np.append(self._x, self._x)
            self._y = np.append(self._y, self._y)
            self._px = np.append(self._px, -self._px)
            self._py = np.append(self._py, self._py)
            self._ct = np.append(self._ct, self._ct)
            self._dE = np.append(self._dE, self._dE)
            self.param.n_macro = self.param.n_macro*2
        if dimlim<=4:
            self._x = np.append(self._x, self._x)
            self._y = np.append(self._y, -self._y)
            self._px = np.append(self._px, self._px)
            self._py = np.append(self._py, self._py)
            self._ct = np.append(self._ct, self._ct)
            self._dE = np.append(self._dE, self._dE)
            self.param.n_macro = self.param.n_macro*2

            self._x = np.append(self._x, self._x)
            self._y = np.append(self._y, self._y)
            self._px = np.append(self._px, self._px)
            self._py = np.append(self._py, -self._py)
            self._ct = np.append(self._ct, self._ct)
            self._dE = np.append(self._dE, self._dE)
            self.param.n_macro = self.param.n_macro*2

        if dimlim<=6:
            self._x = np.append(self._x, self._x)
            self._y = np.append(self._y, self._y)
            self._px = np.append(self._px, self._px)
            self._py = np.append(self._py, self._py)
            self._ct = np.append(self._ct, -self._ct)
            self._dE = np.append(self._dE, self._dE)
            self.param.n_macro = self.param.n_macro*2

            self._x = np.append(self._x, self._x)
            self._y = np.append(self._y, self._y)
            self._px = np.append(self._px, self._px)
            self._py = np.append(self._py, self._py)
            self._ct = np.append(self._ct, self._ct)
            self._dE = np.append(self._dE, -self._dE)
            self.param.n_macro = self.param.n_macro*2

       
        

    def initialize_distribution(self, dim=6, seed=None, mirror=False):

        if seed is not None:
            np.random.seed(seed)

        if self._ndim<dim:
            self._ndim=dim
        if self.param.n_macro==0:
            print("Warning, zero macro particle requested")


        moment_mat=np.eye(dim)
        dis=[]
        if dim >= 2:
            self._x[:] = lf.randn_cutoff(self.param.gaussian_cutoff, int(self.param.n_macro)) #*self.param.initial_beamsize[0]
            self._px[:] = lf.randn_cutoff(self.param.gaussian_cutoff, int(self.param.n_macro)) #*self.param.initial_beamsize[0]/self.param.beam_beta[0]
            #element.thin_quad_2d(self._x,self._px,self.param.beam_alpha[0]/self.param.beam_beta[0])
            moment_mat[0,1]=np.mean(self._x*self._px)
            dis.append(self._x)
            dis.append(self._px)
        if dim >= 4:
            self._y[:] = lf.randn_cutoff(self.param.gaussian_cutoff, int(self.param.n_macro))   #*self.param.initial_beamsize[1]
            self._py[:] = lf.randn_cutoff(self.param.gaussian_cutoff, int(self.param.n_macro))   #*self.param.initial_beamsize[1]/self.param.beam_beta[1]
            #element.thin_quad_2d(self._y,self._py,self.param.beam_alpha[1]/self.param.beam_beta[1])
            moment_mat[0, 2] = np.mean(self._x * self._y)
            moment_mat[0, 3] = np.mean(self._x * self._py)
            moment_mat[1, 2] = np.mean(self._px * self._y)
            moment_mat[1, 3] = np.mean(self._px * self._py)
            moment_mat[2, 3] = np.mean(self._y * self._py)
            dis.append(self._y)
            dis.append(self._py)
        if dim >= 6:
            self._ct[:] = lf.randn_cutoff(self.param.gaussian_cutoff, int(self.param.n_macro))    #*self.param.initial_beamsize[2]
            self._dE[:] = lf.randn_cutoff(self.param.gaussian_cutoff, int(self.param.n_macro))      #*self.param.initial_energy_spread
            moment_mat[0, 4] = np.mean(self._x * self._ct)
            moment_mat[0, 5] = np.mean(self._x * self._dE)
            moment_mat[1, 4] = np.mean(self._px * self._ct)
            moment_mat[1, 5] = np.mean(self._px * self._dE)
            moment_mat[2, 4] = np.mean(self._y * self._ct)
            moment_mat[2, 5] = np.mean(self._y * self._dE)
            moment_mat[3, 4] = np.mean(self._py * self._ct)
            moment_mat[3, 5] = np.mean(self._py * self._dE)
            moment_mat[4, 5] = np.mean(self._ct * self._dE)
            dis.append(self._ct)
            dis.append(self._dE)

        moment_mat=moment_mat+moment_mat.transpose()-np.eye(dim)
        dis = [self._x, self._px, self._y, self._py, self._ct, self._dE]
        evl, ev = np.linalg.eig(moment_mat)
        vec = ev @ np.diag(1 / np.sqrt(evl)) @ ev.transpose()
        if dim==2:
            self._x[:], self._px[:] = vec.dot(np.array(dis))
            self._x *= self.param.initial_beamsize[0]
            self._px *= self.param.initial_beamsize[0] / self.param.beam_beta[0]
            element.thin_quad_2d(self._x,self._px,self.param.beam_alpha[0]/self.param.beam_beta[0])
        elif dim==4:
            self._x[:], self._px[:], self._y[:], self._py[:]=vec.dot(np.array(dis))
            self._x *= self.param.initial_beamsize[0]
            self._px *= self.param.initial_beamsize[0] / self.param.beam_beta[0]
            self._y *= self.param.initial_beamsize[1]
            self._py *= self.param.initial_beamsize[1] / self.param.beam_beta[1]
            element.thin_quad_2d(self._x, self._px, self.param.beam_alpha[0] / self.param.beam_beta[0])
            element.thin_quad_2d(self._y, self._py, self.param.beam_alpha[1] / self.param.beam_beta[1])
        elif dim==6:
            self._x[:], self._px[:], self._y[:], self._py[:], self._ct[:], self._dE[:] = vec.dot(np.array(dis))
            self._x *= self.param.initial_beamsize[0]
            self._px *= self.param.initial_beamsize[0] / self.param.beam_beta[0]
            self._y *= self.param.initial_beamsize[1]
            self._py *= self.param.initial_beamsize[1] / self.param.beam_beta[1]
            self._ct *= self.param.initial_beamsize[2]
            self._dE *= self.param.initial_energy_spread
            element.thin_quad_2d(self._x, self._px, self.param.beam_alpha[0] / self.param.beam_beta[0])
            element.thin_quad_2d(self._y, self._py, self.param.beam_alpha[1] / self.param.beam_beta[1])
            if self.param.rf_voltage!=0:
                for i in range(int(3/self.param.tune[-1])):
                    self.one_turn_z_RF()
                ztphi = 2 * np.pi * self.param.rf_freq / const.c
                maxphase=np.abs(self.param.rf_unstable_phase-self.param.rf_phase)
                mask=np.abs(self._ct*ztphi)>np.abs(self.param.rf_unstable_phase-self.param.rf_phase)
                self._ct[mask]=np.mod(np.abs(self._ct[mask]*ztphi), maxphase)/ztphi
                if self.param.rf_phase > np.pi/2:
                    self._ct[mask]*=-1
                self._dE[mask]=0.0
                for i in range(int(5/self.param.tune[-1])):
                    self.one_turn_z_RF()

        if mirror is True:
            self.mirror_copy_distribution()



if __name__=='__main__':
    test=Beam()
    pass

