#!/usr/bin/env python3

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv
from pyphare.pharesee.hierarchy import get_times_from_h5
from tests.diagnostic import all_timestamps
from pyphare.pharesee.run import Run
import os
from pyphare.pharein.global_vars import sim



import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.use('Agg')



def config():

    Simulation(
        time_step = 0.001,
        final_time = 80.,
        boundary_types = ("periodic", "periodic"),
        cells = (1000, 200),
        dl = (1., 1.),
        refinement_boxes = {},
        resistivity = 0.001,
        hyper_resistivity = 0.001,
        diag_options = {"format" : "phareh5",
                        "options" : {"dir" : "01b",
                                     "mode" : "overwrite"}},
        restart_options={"dir" : "checks",
                         "mode" : "overwrite",
                         "timestamps" : [40., 50., 60., 70., 80.],
                         "restart_time":30.},
    )

    def densityMain(x, y):
        return 1.

    def densityBeam(x, y):
        return .01

    def bx(x, y):
        return 1.

    def by(x, y):
        return 0.

    def bz(x, y):
        return 0.

    def vB(x, y):
        return 50.

    def v0(x, y):
        return 0.

    def vthMain(x, y):
        return np.sqrt(5.0)

    def vthBeam(x, y):
        return np.sqrt(1.0)

    MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        main={"charge": 1,
              "mass": 1,
              "density": densityMain,
              "vbulkx": v0,
              "vbulky": v0,
              "vbulkz": v0,
              "vthx": vthMain,
              "vthy": vthMain,
              "vthz": vthMain,
              "nbr_part_per_cell": 500,
              "init": {"seed": 12}},
        beam={"charge": 1,
              "mass": 1,
              "density": densityBeam,
              "vbulkx": vB,
              "vbulky": v0,
              "vbulkz": v0,
              "vthx": vthBeam,
              "vthy": vthBeam,
              "vthz": vthBeam,
              "nbr_part_per_cell": 100,
              "init": {"seed": 12}},
    )

    ElectronModel(closure="isothermal", Te=1.0)

    sim = ph.global_vars.sim

    dt = 100.*sim.time_step
    nt = (sim.final_time)/dt+1
    timestamps_fine = dt * np.arange(nt)

    dt = 1000.*sim.time_step
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
                 compute_timestamps=timestamps_fine
             )

    poplist = ["main", "beam"]
    for pop in poplist:
        for quantity in ["density", "flux"]:
            FluidDiagnostics(quantity=quantity,
                             write_timestamps=timestamps_fine,
                             compute_timestamps=timestamps_fine,
                             population_name=pop)

        for quantity in ['domain']: #, 'levelGhost', 'patchGhost']:
            ParticleDiagnostics(quantity=quantity,
                                compute_timestamps=timestamps_coarse,
                                write_timestamps=timestamps_coarse,
                                population_name=pop)





def main():

    config()
    sim = Simulator(gv.sim)
    sim.initialize()
    sim.run()


if __name__=="__main__":
    main()

