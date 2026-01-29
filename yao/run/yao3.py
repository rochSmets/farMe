## 03/12/2025 
## what if there's no applied initial velocity, but only thermal expansion?

#!/usr/bin/env python3


import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv
# from pyphare.pharesee.hierarchy import get_times_from_h5
# from tests.diagnostic import all_timestamps
# from pyphare.pharesee.run import Run
# import os
# from pyphare.pharein.global_vars import sim



# import matplotlib.pyplot as plt
# import matplotlib as mpl
import numpy as np
# mpl.use('Agg')


center = [[20., 20.], [60., 20.]]
r_density = [15.0, 15.0]
r_velocity = [15.0, 15.0]
w_velocity = [0.5, 0.5]
nMain = 0.9
nBack = 0.1

Epsilon = 1e-8


def slot(x):
    """ return a square slot for -1 < x < +1 """
    return np.heaviside(x+1, 0.5)-np.heaviside(x-1, 0.5)


def polynom(x):
    """ return a bell shaped curve for -1 < x < +1 """
    return -6*np.abs(x)**5+15*x**4-10*np.abs(x)**3+1


def profile(x):
    return slot(x)*polynom(x)


def density(x, y):
    assert x.shape == y.shape
    funcs = np.zeros((x.shape[0], 2))
    for i in range(2):  # because 2 plasma plumes
        x_ = (x-center[i][0])/r_density[0]
        y_ = (y-center[i][1])/r_density[1]
        r_ = np.sqrt(x_**2+y_**2)
        funcs[:, i] = nMain*profile(r_)
    return nBack+funcs.sum(axis=1)


def bx(x, y):
    return 0.0


def by(x, y):
    return 0.0


def bz(x, y):
    return 1.0


def v0x(x, y):
    assert x.shape == y.shape
    funcs = np.zeros((x.shape[0], 2))
    for i in range(2):
        x_ = (x-center[i][0])/r_velocity[0]
        y_ = (y-center[i][1])/r_velocity[1]
        r_ = np.sqrt(x_**2+y_**2)
        s_ = (r_-1)*w_velocity[i]
        funcs[:, i] = profile(s_)*x_/r_
    return funcs.sum(axis=1)


def v0y(x, y):
    assert x.shape == y.shape
    funcs = np.zeros((x.shape[0], 2))
    for i in range(2):
        x_ = (x-center[i][0])/r_velocity[0]
        y_ = (y-center[i][1])/r_velocity[1]
        r_ = np.sqrt(x_**2+y_**2)
        s_ = (r_-1)*w_velocity[i]
        funcs[:, i] = profile(s_)*y_/r_
    return funcs.sum(axis=1)


def v0z(x, y):
    assert x.shape == y.shape
    return 0.0


def vth(x, y):
    return np.sqrt(0.2)

def density_amb(x, y):
    return 1.0

def v0_amb(x, y):
    return 0.

def vth_amb(x, y):
    return np.sqrt(0.1)

x1_ = 5.
x2_ = 30.

def density_s1(x, y):
    return np.piecewise(x, [x < x1_, ((x1_ < x) & (x < x2_)), x > x2_], [lambda x:0.0, lambda x:2.0, lambda x:0.0])

def v_s1_x(x, y):
    return np.piecewise(x, [x < x1_, ((x1_ < x) & (x < x2_)), x > x2_], [lambda x:0.0, lambda x:10.0, lambda x:0.0])

def v_s1_y(x, y):
    return 0.

def v_s1_z(x, y):
    return 0.

x3_ = 35.
x4_ = 60.
def density_s2(x, y):
    return np.piecewise(x, [x < x3_, ((x3_ < x) & (x < x4_)), x > x4_], [lambda x:0.0, lambda x:2.0, lambda x:0.0])

def v_s2_x(x, y):
    return np.piecewise(x, [x < x3_, ((x3_ < x) & (x < x4_)), x > x4_], [lambda x:0.0, lambda x:-10.0, lambda x:0.0])

def v_s2_y(x, y):
    return 0.

def v_s2_z(x, y):
    return 0.


def config():

    Simulation(
        time_step = 0.001,
        final_time = 1.0,
        boundary_types = ("periodic", "periodic"),
        cells = (650, 300),
        dl = (0.1, 0.1),
        refinement_boxes = {},
        tag_buffer="10",
        resistivity = 0.001,
        hyper_resistivity = 0.001,
        diag_options = {"format" : "phareh5",
                        "options" : {"dir" : "ds_t1",
                                     "mode" : "overwrite"}},
        #restart_options={"dir" : "checks",
        #                 "mode" : "overwrite",
        #                 "timestamps" : [10., 20., 30., 40., 50., 60., 70., 80., 90., 100.],
        #                 #"restart_time":30.},
        #interp_order = 1,
    )

    MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        s1={"charge": 1,
              "mass": 1,
              "density": density_s1,
              "vbulkx": v_s1_x, 
              "vbulky": v_s1_y, 
              "vbulkz": v_s1_z,
              "vthx": vth,
              "vthy": vth,
              "vthz": vth,
              "nbr_part_per_cell": 1,
              "init": {"seed": 12}},
        s2={"charge": 1,
              "mass": 1,
              "density": density_s2,
              "vbulkx": v_s2_x,
              "vbulky": v_s2_y,
              "vbulkz": v_s2_z,
              "vthx": vth,
              "vthy": vth,
              "vthz": vth,
              "nbr_part_per_cell": 1,
              "init": {"seed": 12}},
        amb={"charge": 1,
              "mass": 1,
              "density": density_amb,
              "vbulkx": v0_amb,
              "vbulky": v0_amb,
              "vbulkz": v0_amb,
              "vthx": vth_amb,
              "vthy": vth_amb,
              "vthz": vth_amb,
              "nbr_part_per_cell": 1,
              "init": {"seed": 12}},
    )

    ElectronModel(closure="isothermal", Te=0.02)

    sim = gv.sim
    start_time = sim.start_time()

    dt = 100.*sim.time_step
    nt = (sim.final_time-start_time)/dt+1
    timestamps_fine = start_time+dt*np.arange(nt)

    dt = 1000.*sim.time_step
    nt = (sim.final_time-start_time)/dt+1
    timestamps_coarse = start_time+dt*np.arange(nt)

    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps_fine,
        )

    for quantity in ["charge_density", "bulkVelocity"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps_fine,
        )
    
    poplist = ["s1","s2", "amb"]
    for pop in poplist:
        for quantity in ["density", "flux"]:
            FluidDiagnostics(quantity=quantity,
                         write_timestamps=timestamps_fine,
                         population_name=pop)

        for quantity in ['domain']: #, 'levelGhost', 'patchGhost']:
            ParticleDiagnostics(quantity=quantity,
                            write_timestamps=timestamps_coarse,
                            population_name=pop)





def main():

    config()
    sim = Simulator(gv.sim)
    sim.initialize()
    sim.run()


if __name__=="__main__":
    main()

