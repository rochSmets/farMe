#!/usr/bin/env python3

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.use('Agg')



def config():

    Simulation(
        smallest_patch_size=20,
        largest_patch_size=60,
        final_time=200,
        time_step=0.001,
        boundary_types="periodic",
        cells=400,
        dl=0.2,
        hyper_resistivity=0.001,
        refinement_boxes={},
        diag_options={"format": "phareh5",
                      "options": {"dir": "bernsteinA",
                                  "mode": "overwrite"}}
    )


    def density(x):
        return 1.

    def bx(x):
        return 0.08715574274765814

    def by(x):
        return 0.9961946980917455

    def bz(x):
        return 0.

    def vb(x):
        return 0.0

    def vth(x):
        return np.sqrt(0.1)


    vvv = {
        "vbulkx": vb, "vbulky": vb, "vbulkz": vb,
        "vthx": vth, "vthy": vth, "vthz": vth
    }



    MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density, **vvv}
    )


    ElectronModel(closure="isothermal", Te=1.0)


    sim = ph.global_vars.sim


    timestamps = np.arange(0, sim.final_time, 0.5)


    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )


    FluidDiagnostics(
        quantity="density",
        write_timestamps=timestamps,
        compute_timestamps=timestamps,
        )


def main():
    config()
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()

if __name__=="__main__":
    main()
