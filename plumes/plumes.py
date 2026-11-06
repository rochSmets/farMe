#!/usr/bin/env python3

# import pyphare.pharein as ph  # lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, \
    FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator

import numpy as np


Centers = [[10., 10.], ]
Radius_B = [5., ]  # radius of the B shell (centered on r=Radius_B)
Radius_N = [8., ]  # radius of the main density profile vanishing at Radius_N
assert(len(Radius_B) == len(Centers))
assert(len(Radius_N) == len(Centers))
Nplumes = len(Radius_B)
dump_dt = 0.2  # time between two dumps
Ni = [0.8, 0.2]  # particle density for Main and Back
Ti = [1.0, 1.0]
Te = 0.2
qi = [1.0, 1.0]
mi = [1.0, 1.0]
widthRatio = 0.2  # the B bubble has a width equal Radius_B*2*widthRatio
Epsilon = 1e-8


def slot(x):
    return np.heaviside(x+1, 0.5)-np.heaviside(x-1, 0.5)  # 1 for x in [-1, 1], 0 elsewhere


def polynom(x):
    return -6*np.abs(x)**5+15*x**4-10*np.abs(x)**3+1


def profile(x):
    return slot(x)*polynom(x)


def nMain(x, y):
    assert x.shape == y.shape
    funcs = np.zeros((x.shape[0], Nplumes))
    for i in range(Nplumes):
        x_ = (x-Centers[i][0])/Radius_N[i]
        y_ = (y-Centers[i][1])/Radius_N[i]
        r_ = np.sqrt(x_**2+y_**2)
        funcs[:, i] = Ni[0]*profile(r_)
    return funcs.sum(axis=1)


def nBack(x, y):
    assert x.shape == y.shape
    return Ni[1]


def bx(x, y):
    assert x.shape == y.shape
    funcs = np.zeros((x.shape[0], Nplumes))
    for i in range(Nplumes):
        x_ = (x-Centers[i][0])/Radius_B[i]
        y_ = (y-Centers[i][1])/Radius_B[i]
        r_ = np.sqrt(x_**2+y_**2)
        s_ = (r_-1)/widthRatio
        X_ = +y_
        Y_ = -x_
        R_ = np.clip(np.sqrt(X_**2+Y_**2), Epsilon, None)
        funcs[:, i] = profile(s_)*X_/R_
    return funcs.sum(axis=1)


def by(x, y):
    assert x.shape == y.shape
    funcs = np.zeros((x.shape[0], Nplumes))
    for i in range(Nplumes):
        x_ = (x-Centers[i][0])/Radius_B[i]
        y_ = (y-Centers[i][1])/Radius_B[i]
        r_ = np.sqrt(x_**2+y_**2)
        s_ = (r_-1)/widthRatio
        X_ = +y_
        Y_ = -x_
        R_ = np.clip(np.sqrt(X_**2+Y_**2), Epsilon, None)
        funcs[:, i] = profile(s_)*Y_/R_
    return funcs.sum(axis=1)


def bz(x, y):
    return 0.


def v0(x, y):
    return 0.


def vtMain(x, y):
    return np.sqrt(Ti[0]/mi[0])


def vtBack(x, y):
    return np.sqrt(Ti[1]/mi[1])


# def Te(x, y):
#     return 0.2
#
# def Pe(x, y):
#     return density(x, y)*Te(x, y)



def config():
    sim = Simulation(
        time_step=0.002,  # {"mode": "adaptive", "cfl_wave": 0.5} is not available for hybrid yet
        final_time=20.0,
        cells=(100, 100),
        dl=(0.2, 0.2),
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={"format": "phareh5",
                      "options": {"dir": "./plumes.run",
                                  "mode": "overwrite"}, },
        # restart_options={"dir": "checks",
        #                  "mode": "overwrite",
        #                  "timestamps": [.2], }
    )

    MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        main={"charge": qi[0],
              "mass": mi[0],
              "density": nMain,
              "vbulkx": v0,
              "vbulky": v0,
              "vbulkz": v0,
              "vthx": vtMain,
              "vthy": vtMain,
              "vthz": vtMain,
              "nbr_part_per_cell": 100,
              "init": {"seed": 12}},
        back={"charge": qi[1],
              "mass": mi[1],
              "density": nBack,
              "vbulkx": v0,
              "vbulky": v0,
              "vbulkz": v0,
              "vthx": vtBack,
              "vthy": vtBack,
              "vthz": vtBack,
              "nbr_part_per_cell": 100,
              "init": {"seed": 12}},
    )

    ElectronModel(closure="isothermal", Te=Te)
    # ElectronModel(closure="polytropic", Pe=Pe, gamma=1.66)

    nt = sim.time_stepper.final_time/dump_dt
    timestamps = dump_dt * np.arange(int(round(nt))+1)
    # timestamps = [.0, .1, .2]

    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
        )

    for quantity in ["charge_density", "bulkVelocity"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            )

    for pop in ["main", "back"]:
        for quantity in ["density", "flux"]:
            FluidDiagnostics(quantity=quantity,
                             write_timestamps=timestamps,
                             population_name=pop)

            # for quantity in ['domain']:  # , 'levelGhost', 'patchGhost']:
            #     ParticleDiagnostics(quantity=quantity,
            #                         write_timestamps=timestamps,
            #                         population_name=pop)

    return sim


def main():

    simulator = Simulator(config())
    simulator.initialize()
    simulator.run()
    simulator.reset()


if __name__ == "__main__":
    main()
