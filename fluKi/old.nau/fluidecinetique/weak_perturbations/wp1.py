#!/usr/bin/env python3


import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv
from pyphare.pharesee.run import Run


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.use('Agg')


from tests.diagnostic import all_timestamps

def density(x):
    return 1.


def bx(x):
    return 1.


def by(x):
    return 0.


def bz(x):
    return 0.0



def T(x):
    return 0.01


def vx(x):
    from pyphare.pharein.global_vars import sim
    L = sim.simulation_domain()[0]
    x0 = 0.5*L
    lamb = 6.
    bubble = 0.8*np.exp(-(x-x0)**2/lamb)
    return bubble


def vy(x):
    return 0.


def vz(x):
    return 0.


def vthx(x):
    return T(x)


def vthy(x):
    return T(x)


def vthz(x):
    return T(x)

vvv = {"vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz }


def config(**kwargs):

    Simulation(
        time_step=0.005,
        final_time=100.,
        boundary_types="periodic",
        hyper_resistivity=0.001,
        cells=512,
        dl=0.25,
        diag_options={"format": "phareh5",
                      "options": {"dir": kwargs["diagdir"],"mode":"overwrite"}}
    )


    MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1,
                 "density": density,
                 "nbr_part_per_cell":200,
                 **vvv}
    )

    ElectronModel(closure="isothermal", Te=0.05)


    sim = ph.global_vars.sim
    dt = sim.time_step*100
    timestamps = np.arange(0,sim.final_time+dt, dt)


    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )


    for quantity in ["density", "bulkVelocity"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            )





def main():
    from pyphare.cpp import cpp_lib
    cpp = cpp_lib()

    config(diagdir="wp3")
    Simulator(gv.sim).run()
    gv.sim = None




if __name__=="__main__":
    main()
