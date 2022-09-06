#!/usr/bin/env python
# coding: utf-8

# In[1]:


import json
import numpy as np
import matplotlib.pyplot as plt


import xobjects as xo
import xpart as xp
import xtrack as xt
import xfields as xf
from cpymad.madx import Madx

from scipy.constants import e, m_p, c
from scipy.constants import physical_constants

#Imports needed in part 2. creating the python_beam
from PyHEADTAIL.particles import generators
from PyHEADTAIL.particles.generators import RFBucketMatcher
from PyHEADTAIL.trackers.rf_bucket import RFBucket
from PyHEADTAIL.particles.rfbucket_matching import ThermalDistribution
#----

from xpart.pyheadtail_interface.pyhtxtparticles import PyHtXtParticles as Pht_XT
from PyHEADTAIL.particles.slicing import UniformBinSlicer
from xfields import SpaceCharge3D

print('done')


# In[2]:


import pyopencl as cl

platform = cl.get_platforms()
for i in range(3):
    device = platform[i].get_devices()
    print(platform)
    print("---")
    print(device)
    print('\n','------------------------------------------------------------------------------------','\n')
    


# In[3]:


class SpaceChargeWrapper(object):
    def __init__(self, spcharge_node, length):
        self.length = length
        self.spcharge_node = spcharge_node
        
    def track(self, beam):
        save = self.spcharge_node.length
        self.spcharge_node.length = self.length
        self.spcharge_node.track(beam)
        self.spcharge_node.length = save


# In[4]:


class XsuiteSimulation(object):
    # All the variables involved
    n_macroparticles = int(1e6) #1e6
    n_turns = 1

    with_SC = True
    losses = False
    N_factor = 1.0 # times nominal intensity

    A = 238 # mass number
    Q = 28 # elementary charges per particle
    Ekin_per_nucleon = 0.2e9 # in eV

    epsx_rms_fin = 35e-6 / 4 # geometrical emittances
    epsy_rms_fin = 15e-6 / 4

    limit_n_rms_x = 2
    limit_n_rms_y = 2
    limit_n_rms_z = 4

    sig_z_rms = 58. / 4. * 0.3/0.33
    sig_dp_rms = 0.5e-3

    sig_z = sig_z_rms
    sig_dp = sig_dp_rms

    sc_mode = "Bunched"
    intensity = 0.625e11 * N_factor

    # fixed grid for space charge
    n_slices_sc = 64
    n_mesh_nodes = 128
    n_mesh_sigma = 12
    n_scnodes = int(1e3)

    # monitor
    n_stored_particles = 1 #10000

    # RF:
    n_cavities = 14
    lag_h10 = 0
    f_h10 = 1.5705898344588984e6
    v_h10 = 58.2e3

    nmass = physical_constants[
        'atomic mass constant energy equivalent in MeV'][0] * 1e-3

    mass = A * nmass * 1e9 * e / c**2 # in kg
    charge = Q * e # in Coul

    Ekin = Ekin_per_nucleon * A
    p0c = np.sqrt(Ekin**2 + 2*Ekin*mass/e * c**2) # in eV

    Etot = np.sqrt(p0c**2 + (mass/e)**2 * c**4) * 1e-9 # in GeV
    p0 = p0c / c * e # in SI units
    gamma = np.sqrt(1 + (p0 / (mass * c))**2)
    beta = np.sqrt(1 - gamma**-2)


    epsx_gauss = epsx_rms_fin * 1.43
    epsy_gauss = epsy_rms_fin * 1.41

    epsn_x = epsx_gauss * beta * gamma
    epsn_y = epsy_gauss * beta * gamma
    
    use_SC_3D = False
    

    def run(self):
        '''
        Main function: it runs the full simulation calling the different
        funtions of the class.
        '''
        #self.context = xo.ContextCpu()
        self.context = xo.ContextPyopencl('0.0')
        print(self.context)
        
        qx, qy, e_seed = 18.9, 18.88, 1 
        qx = float(qx)
        qy = float(qy)
        e_seed = int(e_seed)
        qqx, qqy = int(np.round((qx%1) * 1e4)), int(np.round((qy%1) * 1e4))
        
        outputpath=str('./')
        filename_error_table = (
                    outputpath +
                    "/error_tables/errors_{qqx}_{qqy}_{eseed:d}".format(
                        qqx=qqx, qqy=qqy, eseed=e_seed))
        
        madx = Madx() # initializes a mad x object
        madx.options.echo = False
        madx.options.warn = False
        madx.options.info = False
        # reading files and importing info to madx obj
        twiss_dict = self.setup_madx(madx, qx, qy, filename_error_table)
        
        #Importing the madx lattice to xsuite
        xtrack_elements =self.lattice_madx_to_xtrack(madx)
        xtrack_elements = xtrack_elements.remove_zero_length_drifts()
        
        #creating particles at pyHEADTAIL
        pyht_beam = self.setup_pyheadtail_particles(twiss_dict,self.n_macroparticles)
        
        #importing the particles to a xpart obj
        xparticles = self.particles_to_xpart(pyht_beam)
        
        #creates the base-sc_node used later in PIC algorithm
        spcharge_node = self.sc_node_creator(pyht_beam)
        
        #locating the scn nodes - not used- positions are copyed from the ther program
        #sc_locations, sc_lengths = self.get_sc_locations(twiss_dict['length'], self.n_scnodes, xtrack_elements)
        s_spacecharge = np.load('/home/pablo/scripts/PIC/s_scnodes.npy')
        
        #setting the sc_nodes
        xtrack_elements = self.set_the_sc_nodes(s_spacecharge, twiss_dict, xtrack_elements, spcharge_node )
        
        #tracking
        tracker = self.trackABC(xtrack_elements, xparticles)
        
        return tracker
    
        
    
    def setup_madx(self,madx, qx, qy, filename_error_table):
        '''
        Reading a file assigns the information to a madx object
        '''
        madx.call('./SIS100RING.seq')
        k1nl_s52qd11_factor = 1.0
        k1nl_s52qd12_factor = 1.0
        seqname = 'sis100cold'

        madx.command.beam(particle='ion', mass=self.A*self.nmass,
                          charge=self.Q, energy=self.Etot)

        assert madx.input('''
                select, flag=seqedit, class=collimator;

                seqedit, sequence={sn};
                    remove, element=selected;
                    flatten;
                endedit;

                select, flag=seqedit, class=marker;
                seqedit, sequence={sn};
                    remove, element=selected;
                    install, element={sn}$START, s=0;
                    flatten;
                endedit;
            '''.format(sn=seqname))

        madx.use(sequence=seqname)

        assert madx.command.select(
            flag='MAKETHIN',
            class_='QUADRUPOLE',
            slice_='9',)

        assert madx.command.select(
            flag='MAKETHIN',
            class_='SBEND',
            slice_='9',)

        assert madx.command.makethin(
            makedipedge=True,
            style='teapot',
            sequence=seqname,)

        madx.call('OpticsYEH_BeamParameters.str')

        madx.use(sequence=seqname)

        madx.input(f'''
                kqd := -0.2798446835;
                kqf := 0.2809756135;
                kqd := -2.81308e-01;
                kqf := 2.81990e-01;
                K1NL_S00QD1D :=  kqd ;
                K1NL_S00QD1F :=  kqf ;
                K1NL_S00QD2F :=  kqf ;
                K1NL_S52QD11 :=  {k1nl_s52qd11_factor}   *   kqd ;
                K1NL_S52QD12 :=  {k1nl_s52qd11_factor}   *   kqf ;
            ''')

        # match normally using MAD-X
        madx.input('''
                    match, sequence={sn};
                    global, sequence={sn}, q1={qx}, q2={qy};
                    vary, name=kqf, step=0.00001;
                    vary, name=kqd, step=0.00001;
                    lmdif, calls=500, tolerance=1.0e-10;
                    endmatch;
            '''.format(qx=qx, qy=qy, sn=seqname)
            )

        twiss = madx.twiss()

        madx.input('cavity_voltage = {}/{}'.format(
            self.v_h10 * 1e-6, self.n_cavities))

        twiss_dict = {
            'x': twiss['x'][0],
            'px': twiss['px'][0],
            'y': twiss['y'][0],
            'py': twiss['py'][0],
            'betx': twiss['betx'][0],
            'alfx': twiss['alfx'][0],
            'bety': twiss['bety'][0],
            'alfy': twiss['alfy'][0],
            'dx': twiss['dx'][0],
            'dpx': twiss['dpx'][0],
            'dy': twiss['dy'][0],
            'dpy': twiss['dpy'][0],
            'length': twiss.summary.length,
            'gammatr': twiss.summary.gammatr, }
        return twiss_dict

    
    
    
    def lattice_madx_to_xtrack(self, madx):
        '''
        this function imports the madX lattice to a xtrack element
        '''
        seqname = 'sis100cold'
        sis100 = getattr(madx.sequence, seqname)
        #lattice transfer and preparation!
        xtrack_elements = xt.Line.from_madx_sequence(sis100, exact_drift=False, install_apertures=False)
        
        return xtrack_elements 
    

    def setup_pyheadtail_particles(self, twiss_table, npart=n_macroparticles):
        '''
        this function initializes the particles in pht
        '''

        x_co = 0
        xp_co = 0
        y_co = 0
        yp_co = 0

        D_x_0 = twiss_table['dx'] * self.beta
        D_y_0 = twiss_table['dy'] * self.beta

        Dp_x_0 = twiss_table['dpx'] * self.beta
        Dp_y_0 = twiss_table['dpy'] * self.beta

        np.random.seed(0)

        pyht_beam = generators.generate_Gaussian6DTwiss(
            npart, self.intensity, self.charge, self.mass,
            twiss_table['length'], self.gamma,
            twiss_table['alfx'], twiss_table['alfy'],
            twiss_table['betx'], twiss_table['bety'],
            1, self.epsn_x, self.epsn_y, 1,
            dispersion_x=None, #D_x_0 if D_x_0 else None,
            dispersion_y=None, #D_y_0 if D_y_0 else None,
            limit_n_rms_x=self.limit_n_rms_x**2,
            limit_n_rms_y=self.limit_n_rms_y**2,
            limit_n_rms_z=self.limit_n_rms_z**2,
        )

        if self.v_h10:
            rfb = RFBucket(
                twiss_table['length'], self.gamma, self.mass,
                self.charge, [twiss_table['gammatr']**-2],
                p_increment=0., harmonic_list=[10],
                voltage_list=[self.v_h10],
                phi_offset_list=[np.pi + self.lag_h10 * np.pi * 2])

            rfb_matcher = RFBucketMatcher(
                rfb, ThermalDistribution, sigma_z=self.sig_z_rms)
            z, dp, _, _ = rfb_matcher.generate(npart)
            pyht_beam.z, pyht_beam.dp = z, dp

        # recentre on 0 to avoid dipolar motion:
        pyht_beam.x -= pyht_beam.mean_x()
        pyht_beam.xp -= pyht_beam.mean_xp()
        pyht_beam.y -= pyht_beam.mean_y()
        pyht_beam.yp -= pyht_beam.mean_yp()
        pyht_beam.z -= pyht_beam.mean_z()
        pyht_beam.dp -= pyht_beam.mean_dp()

        # add dispersive contribution to coordinates:
        pyht_beam.x += D_x_0 * pyht_beam.dp
        pyht_beam.y += D_y_0 * pyht_beam.dp
        # also need to add D'_{x,y} to momenta:
        pyht_beam.xp += Dp_x_0 * pyht_beam.dp
        pyht_beam.yp += Dp_y_0 * pyht_beam.dp

        # PyHT generates around 0, closed orbit offset:
        pyht_beam.x += x_co
        pyht_beam.xp += xp_co
        pyht_beam.y += y_co
        pyht_beam.yp += yp_co

        return pyht_beam

    
    def particles_to_xpart(self,pyht_beam):
        '''
        It imports the pyHEADTAIL particles into a xpart
        object 
        important:new atribute in xsuite weight necessary 
        por the correct behavior of the sc_nodes
        '''

        particles=Pht_XT.from_pyheadtail(pyht_beam)
        x=particles.x
        py=particles.py
        px=particles.px
        y=particles.y
        z=particles.z
        delta=particles.dp

        X_particles = xp.Particles(_context=self.context,
                                   mass0= self.A*self.nmass*1e9,
                                   q0=pyht_beam.charge/e,
                                   p0c=self.p0c, #eV
                                   x=x,
                                   px=px,
                                   y=y,
                                   py=py,
                                   zeta=z,
                                   delta=delta,
                                   weight=pyht_beam.charge_per_mp / pyht_beam.charge,
                                  )
        return(X_particles)

    def sc_node_creator(self, pyht_beam):
        '''
        It creates the sc node object that will be used
        on the PIC algorithm later for calculating the 
        effect of the space charge interactions
        '''
        
        sig_x = pyht_beam.sigma_x()
        sig_y = pyht_beam.sigma_y()
        slicer_sc = UniformBinSlicer(self.n_slices_sc, n_sigma_z=1.1 *self.limit_n_rms_z)
        slices = pyht_beam.get_slices(slicer_sc)
        
        x_lim = self.n_mesh_sigma*sig_x
        y_lim = self.n_mesh_sigma*sig_y
        z_lim_tup = [slices.z_cut_tail, slices.z_cut_head]
        #mesh_3d.z0 / pyht_beam.gamma, (mesh_3d.z0 + mesh_3d.nz*mesh_3d.dz) / pyht_beam.gamma

        delta_x = 2 * x_lim / self.n_mesh_nodes
        delta_y = 2 * y_lim / self.n_mesh_nodes
        delta_z = np.diff(z_lim_tup)[0] / self.n_slices_sc

        z_lim_tup[1] -= delta_z

        spcharge_node = SpaceCharge3D(
                _context=self.context,
                length=1.0, update_on_track=True, apply_z_kick=False,
                x_range=(-x_lim, x_lim - delta_x),
                y_range=(-y_lim, y_lim - delta_y),
                z_range=z_lim_tup,
                nx=self.n_mesh_nodes, ny=self.n_mesh_nodes, nz=self.n_slices_sc,
                solver='FFTSolver3D' if self.use_SC_3D else 'FFTSolver2p5D',
                gamma0=pyht_beam.gamma)
        
        return spcharge_node
    
    def get_sc_locations(self, circumference, n_scnodes, line):
        '''
        It calculates the appropriate location for the sc nodes
        '''
        l_target = circumference / n_scnodes
        length_fuzzy = l_target / 2
        s_elements=np.array(line.get_s_elements())
        length_target=s_elements[-1]/float(n_scnodes)
        s_targets = np.arange(0, s_elements[-1], length_target)
        sc_locations = []
        for s in s_targets:
            idx_closest = (np.abs(s_elements - s)).argmin()
            s_closest = s_elements[idx_closest]
            if abs(s - s_closest) < length_fuzzy / 2.0:
                sc_locations.append(s_closest)
            else:
                sc_locations.append(s)
        sc_lengths = np.diff(sc_locations).tolist() + [s_elements[-1] - sc_locations[-1]]

        return sc_locations, sc_lengths
    
    def set_the_sc_nodes(self, s_spacecharge, twiss_dict, xtrack_elements, spcharge_node ):
        '''
        It sets the nodes on the correct place
        A especial objet called SpaceChargeWrapper is used to allow the 
        use of the same scnode just changin the length what is
        much faster 
        '''
        
        sc_lengths = np.diff(s_spacecharge).tolist() + [twiss_dict['length'] - s_spacecharge[-1]]
        sc_elements = []
        sc_names = []
        for ii, (scl, ss) in enumerate(zip(sc_lengths, s_spacecharge)):
            sc_elements.append(SpaceChargeWrapper(spcharge_node, scl))
            sc_names.append(f'spacecharge_{ii}')        
            xtrack_elements.insert_element(name=sc_names[-1], element=sc_elements[-1],
                                    at_s=ss) #s_tol=tol_spacecharge_position
        return xtrack_elements

    def trackABC(self, xtrack_elements, xparticles):
        tracker = xt.Tracker(_context=self.context, line=xtrack_elements)
        tracker.track(xparticles, num_turns=self.n_turns, turn_by_turn_monitor= None )
        #tracker.track(xparticles, num_turns=self.n_turns, turn_by_turn_monitor= None )
        return tracker


# In[5]:


simulation=XsuiteSimulation()


# In[6]:


import cProfile


# In[7]:


cProfile.run('simulation.run()')






