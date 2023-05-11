#!/usr/bin/env python3

import pyphare.pharein as ph
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
        final_time = 1,
        boundary_types = ("periodic", "periodic"),
        cells = (400,200),
        dl = (0.2,0.2),
        refinement_boxes = {},
        resistivity = 0.001,
        hyper_resistivity = 0.001,
        diag_options = {"format": "phareh5",
                        "options": {"dir": "run", "mode": "overwrite"}}
    )

    sim = ph.global_vars.sim
    L = sim.simulation_domain()[0]

    def densityToRight(x, y):
        center = 0.25*L
        width = 0.1*L
        return 0.8*np.exp((x-center)/width)

    def densityToLeft(x, y):
        center = 0.75*L
        width = 0.1*L
        return 0.8*np.exp((x-center)/width)

    def densityBackground(x, y):
        return 0.2

    def bx(x, y):
        return 1.0

    def by(x, y):
        return 0.0

    def bz(x, y):
        return 0.0

    def vToRight(x, y):
        return 1.0

    def vToLeft(x, y):
        return -1.0

    def v0(x, y):
        return 0.0

    def vth(x, y):
        return np.sqrt(1.0)

    MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        background={"charge": 1,
                    "mass": 1,
                    "density": densityBackground,
                    "vbulkx": v0,
                    "vbulky": v0,
                    "vbulkz": v0,
                    "vthx": vth,
                    "vthy": vth,
                    "vthz": vth,
                    "nbr_part_per_cell": 40,
                    "init": {"seed": 12}},
        toright={"charge": 1,
                 "mass": 1,
                 "density": densityToRight,
                 "vbulkx": vToRight,
                 "vbulky": v0,
                 "vbulkz": v0,
                 "vthx": vth,
                 "vthy": vth,
                 "vthz": vth,
                 "nbr_part_per_cell": 40,
                 "init": {"seed": 12}},
        toleft={"charge": 1,
                "mass": 1,
                "density": densityToRight,
                "vbulkx": vToRight,
                "vbulky": v0,
                "vbulkz": v0,
                "vthx": vth,
                "vthy": vth,
                "vthz": vth,
                "nbr_part_per_cell": 40,
                "init": {"seed": 12}},
    )

    ElectronModel(closure="isothermal", Te=1.0)

    timestamps_particles = np.linspace(0, 1, 11)
    timestamps_fields =  np.linspace(0, 1, 11)

    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps_fields,
            compute_timestamps=timestamps_fields,
        )
    for quantity in ["density", "bulkVelocity"]:
             FluidDiagnostics(
                 quantity=quantity,
                 write_timestamps=timestamps_fields,
                 compute_timestamps=timestamps_fields
             )

    poplist = ["background", "toright", "toleft"]
    for pop in poplist:
        for quantity in ["density", "flux"]:
            FluidDiagnostics(quantity=quantity,
                             write_timestamps=timestamps_fields,
                             compute_timestamps=timestamps_fields,
                             population_name=pop)

        for quantity in ['domain']: #, 'levelGhost', 'patchGhost']:
            ParticleDiagnostics(quantity=quantity,
                                compute_timestamps=timestamps_particles,
                                write_timestamps=timestamps_particles,
                                population_name=pop)




def main():

    config()
    sim = Simulator(gv.sim)
    sim.initialize()
    sim.run()


if __name__=="__main__":
    main()

