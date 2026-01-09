#!/usr/bin/env python3

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv


import matplotlib as mpl
import numpy as np
mpl.use('Agg')


cells = (200, 10)
dl = (0.2, 0.2)
L = tuple([n*d for n, d in zip(cells, dl)])
xs1_ = (0.2*L[0], 0.4*L[0])
xs2_ = (0.6*L[0], 0.8*L[0])

nB_ = 0.2
nLeft_ = 1.0
nRight_ = 1.0

vB_ = 0.2

TBack_ = 0.01
TLeft_ = 0.08
TRight_ = 0.08

Te = 0.02


def config():

    Simulation(
        time_step = 0.001,
        final_time = 10.,
        boundary_types = ("periodic", "periodic"),
        cells = cells,
        dl = dl,
        refinement_boxes = {},
        tag_buffer = "10",
        resistivity = 0.001,
        hyper_resistivity = 0.001,
        diag_options = {"format" : "phareh5",
                        "options" : {"dir" : "run",
                                     "mode" : "overwrite"}},
        #restart_options={"dir" : "checks",
        #                 "mode" : "overwrite",
        #                 "timestamps" : [10., 20., 30., 40., 50., 60., 70., 80., 90., 100.],
        #                 #"restart_time":30.},
        interp_order = 1,
    )

    def densityBack(x, y):
        return np.piecewise(x, [x < xs1_[0], ((xs1_[0] < x) & (x < xs1_[1])), ((x > xs1_[1]) & (x < xs2_[0])), ((xs2_[0] < x) & (x < xs2_[1])), x > xs2_[1]] , [lambda x:nB_, lambda x:0.0 , lambda x:nB_ , lambda x:0.0 , lambda x:nB_])

    def densityLeft(x, y):
        return np.piecewise(x, [x < xs1_[0], ((xs1_[0] < x) & (x < xs1_[1])), x > xs1_[1]], [lambda x:0.0, lambda x:nLeft_, lambda x:0.0])

    def densityRight(x, y):
        return np.piecewise(x, [x < xs2_[0], ((xs2_[0] < x) & (x < xs2_[1])), x > xs2_[1]], [lambda x:0.0, lambda x:nRight_, lambda x:0.0])

    def bx(x, y):
        return 0.

    def by(x, y):
        return 0.

    def bz(x, y):
        return 1.

    def vBLeft(x, y):
        return np.piecewise(x, [x < xs1_[0], ((xs1_[0] < x) & (x < xs1_[1])), x > xs1_[1]], [lambda x:0.0, lambda x:+vB_, lambda x:0.0])

    def vBRight(x, y):
        return np.piecewise(x, [x < xs2_[0], ((xs2_[0] < x) & (x < xs2_[1])), x > xs2_[1]], [lambda x:0.0, lambda x:-vB_, lambda x:0.0])

    def v0(x, y):
        return 0.

    def vthBack(x, y):
        return np.sqrt(TBack_)

    def vthLeft(x, y):
        return np.sqrt(TLeft_)

    def vthRight(x, y):
        return np.sqrt(TRight_)

    MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        back={"charge": 1,
              "mass": 1,
              "density": densityBack,
              "vbulkx": v0,
              "vbulky": v0,
              "vbulkz": v0,
              "vthx": vthBack,
              "vthy": vthBack,
              "vthz": vthBack,
              "nbr_part_per_cell": 100,
              "init": {"seed": 12}},
        left={"charge": 1,
              "mass": 1,
              "density": densityLeft,
              "vbulkx": vBLeft,
              "vbulky": v0,
              "vbulkz": v0,
              "vthx": vthLeft,
              "vthy": vthLeft,
              "vthz": vthLeft,
              "nbr_part_per_cell": 100,
              "init": {"seed": 12}},
        right={"charge": 1,
              "mass": 1,
              "density": densityRight,
              "vbulkx": vBRight,
              "vbulky": v0,
              "vbulkz": v0,
              "vthx": vthRight,
              "vthy": vthRight,
              "vthz": vthRight,
              "nbr_part_per_cell": 100,
              "init": {"seed": 12}},
    )

    ElectronModel(closure="isothermal", Te=Te)

    sim = ph.global_vars.sim
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
    for quantity in ["density", "charge_density", "bulkVelocity"]:
             FluidDiagnostics(
                 quantity=quantity,
                 write_timestamps=timestamps_fine,
             )

    poplist = ["back", "left", "right"]
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

