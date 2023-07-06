#!/usr/bin/env python3

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics,FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')



def config():
    L=0.5
    Simulation(
        time_step=0.005,
        final_time=20.,
        cells=(400, 100),
        dl=(0.40, 0.40),
        refinement="tagging",
        loadbalancing="nppc",
        #refinement_boxes=None,
        max_nbr_levels=2,
        nesting_buffer=1,
        clustering="tile",
        tag_buffer="10",
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={"format": "phareh5",
                      "options": {"dir": "02a",
                                  "mode":"overwrite"}
                     },
        # restart_options={"dir":"checks", "mode":"overwrite",
        #                  "timestamps":[20.,],
        #                  #"restart_time":80.
        #                 }
    )

    x0_ = 0.5
    y1_ = 0.3
    y2_ = 0.7

    def density(x, y):
        from pyphare.pharein.global_vars import sim
        Ly = sim.simulation_domain()[1]
        return 0.4 + 1./np.cosh((y-Ly*y1_)/L)**2 + 1./np.cosh((y-Ly*y2_)/L)**2


    def S(y, y0, l):
        return 0.5*(1. + np.tanh((y-y0)/l))


    def by(x, y):
        from pyphare.pharein.global_vars import sim
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.
        dB = 0.1

        x0 = (x - x0_ * Lx)
        y1 = (y - y1_ * Ly)
        y2 = (y - y2_ * Ly)

        dBy1 =  2*dB*x0 * np.exp(-(x0**2 + y1**2)/(sigma)**2)
        dBy2 = -2*dB*x0 * np.exp(-(x0**2 + y2**2)/(sigma)**2)

        return dBy1 + dBy2


    def bx(x, y):
        from pyphare.pharein.global_vars import sim
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.
        dB = 0.1

        x0 = (x - x0_ * Lx)
        y1 = (y - y1_ * Ly)
        y2 = (y - y2_ * Ly)

        dBx1 = -2*dB*y1 * np.exp(-(x0**2 + y1**2)/(sigma)**2)
        dBx2 =  2*dB*y2 * np.exp(-(x0**2 + y2**2)/(sigma)**2)

        v1=-1
        v2=1.
        return v1 + (v2-v1)*(S(y,Ly*y1_,L) -S(y, Ly*y2_,L)) + dBx1 + dBx2


    def bz(x, y):
        return 0.


    def b2(x, y):
        return bx(x,y)**2 + by(x, y)**2 + bz(x, y)**2


    def T(x, y):
        K = 0.7
        temp = 1./density(x, y)*(K - b2(x, y)*0.5)
        assert np.all(temp >0)
        return temp

    def v0(x, y):
        return 0.

    def vth(x, y):
        return np.sqrt(T(x, y))

    def vthz(x, y):
        return np.sqrt(T(x, y))


    MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1,
                 "mass": 1.,
                 "density": density,
                 "vbulkx":v0,
                 "vbulky":v0,
                 "vbulkz":v0,
                 "vthx":vth,
                 "vthy":vth,
                 "vthz":vth,
                 "nbr_part_per_cell": 500,
                 "init": {"seed": 12}},
    )

    ElectronModel(closure="isothermal", Te=0.0)



    sim = ph.global_vars.sim

    dt = 200.*sim.time_step
    nt = (sim.final_time)/dt+1
    timestamps_fine = dt * np.arange(nt)

    dt = 2000.*sim.time_step
    nt = (sim.final_time)/dt+1
    timestamps_coarse = dt * np.arange(nt)



    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps_fine,
            compute_timestamps=timestamps_fine,
        )


    for quantity in ["density", "bulkVelocity"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps_fine,
            compute_timestamps=timestamps_fine,
            )


    for quantity in ["density", "flux"]:
        FluidDiagnostics(quantity=quantity,
                         write_timestamps=timestamps_fine,
                         compute_timestamps=timestamps_fine,
                         population_name="protons")


    for quantity in ['domain']: #, 'levelGhost', 'patchGhost']:
        ParticleDiagnostics(quantity=quantity,
                            compute_timestamps=timestamps_coarse,
                            write_timestamps=timestamps_coarse,
                            population_name="protons")



def main():

    config()
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()


if __name__=="__main__":
    main()
