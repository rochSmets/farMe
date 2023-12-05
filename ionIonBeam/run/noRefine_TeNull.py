
import pyphare.pharein as ph
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

    # configure the simulation
    # most unstable mode at k=0.19, that is lambda = 33
    # hence the length of the box is 33

    Simulation(
        smallest_patch_size=20,
        largest_patch_size=60,
        final_time=60,
        time_step=0.0005,
        boundary_types="periodic",
        cells=165,
        dl=0.2,
        hyper_resistivity=0.01,
        refinement_boxes={},
        diag_options={"format": "phareh5",
                      "options": {"dir": "ionIonBeam_noRefine_TeNull",
                                  "mode": "overwrite"}}
    )


    def densityMain(x):
        return 1.

    def densityBeam(x):
        return .01

    def bx(x):
        return 1.

    def by(x):
        return 0.

    def bz(x):
        return 0.

    def vB(x):
        return 5.

    def v0(x):
        return 0.

    def vth(x):
        return np.sqrt(0.1)


    vMain = {
        "vbulkx": v0, "vbulky": v0, "vbulkz": v0,
        "vthx": vth, "vthy": vth, "vthz": vth
    }


    vBulk = {
        "vbulkx": vB, "vbulky": v0, "vbulkz": v0,
        "vthx": vth, "vthy": vth, "vthz": vth
    }


    MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        main={"charge": 1, "density": densityMain, **vMain},
        beam={"charge": 1, "density": densityBeam, **vBulk}
    )


    ElectronModel(closure="isothermal", Te=np.sqrt(0.0))


    sim = ph.global_vars.sim


    timestamps = np.arange(0, sim.final_time, 0.05)


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


    for popname in ("main", "beam"):
             for name in ["domain", "levelGhost", "patchGhost"]:
                 ParticleDiagnostics(quantity=name,
                                     compute_timestamps=timestamps,
                                     write_timestamps=timestamps,
                                     population_name=popname)



def main():
    config()
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()

if __name__=="__main__":
    main()
