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



def config():

    Ti = 0.01
    Te = 0.1
    n1 = 0.06
    L = 40.0
    gamma_e = 1
    gamma_i = 3

    cs = np.sqrt(gamma_e*Te+gamma_i*Ti)

    Simulation(
        time_step=0.005,
        final_time=30.0,
        boundary_types="periodic",
        hyper_resistivity=0.001,
        cells=200,
        dl=0.2,
        diag_options={"format": "phareh5",
                      "options": {"dir": "test_iaw",
                                  "mode":"overwrite"}
                     }
    )

    def density(x):
        return 1.0+np.sin(2*np.pi*x/L)*n1

    def bx(x):
        return 1.0

    def by(x):
        return 0.0

    def bz(x):
        return 0.0

    def v1(x):
        return np.sin(2*np.pi*x/L)*n1*np.sqrt(gamma_e*Te+gamma_i*Ti)

    def v0(x):
        return 0.

    def vth(x):
        return np.sqrt(Ti)

    vvv = {"vbulkx": v1,
           "vbulky": v0,
           "vbulkz": v0,
           "vthx": vth,
           "vthy": vth,
           "vthz": vth }

    MaxwellianFluidModel(bx=bx,
                         by=by,
                         bz=bz,
                         protons={"charge": 1,
                                  "density": density,
                                  "nbr_part_per_cell": 1000,
                                  **vvv}
                        )

    ElectronModel(closure="isothermal", Te=Te)

    sim = ph.global_vars.sim
    dt = sim.time_step*400
    timestamps = np.arange(0,sim.final_time+dt, dt)


    # for quantity in ["E", "B"]:
    #     ElectromagDiagnostics(
    #         quantity=quantity,
    #         write_timestamps=timestamps,
    #     )

    # for quantity in ["density", "charge_density", "mass_density", "flux", "bulkVelocity", "momentum_tensor"]:
    for quantity in ["charge_density", "bulkVelocity"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            )

    # for quantity in ["density", "flux"]:
    #     FluidDiagnostics(
    #         quantity=quantity,
    #         write_timestamps=timestamps,
    #         population_name="protons",
    #         )

    # for quantity in ['domain']:  # , 'levelGhost', 'patchGhost']:
    #     ParticleDiagnostics(quantity=quantity,
    #                         write_timestamps=timestamps,
    #                         population_name="protons")


def sine_func(x, A, B, C):
    return A * np.sin(k * x + B) + C


def main():
    # from pyphare.cpp import cpp_lib
    # import sys

    # cpp = cpp_lib()

    config()
    Simulator(gv.sim).run()
    gv.sim = None



if __name__=="__main__":
    main()
