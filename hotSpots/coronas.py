#!/usr/bin/env python3

import pyphare.pharein as ph  # lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, \
    FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv

import numpy as np
import matplotlib as mpl
mpl.use('Agg')


Centers = [[10., 10.], [30., 10.]]
Radius = [6., 8.]
nMain = 0.8
nBack = 0.2
widthRatio = 0.3  # the B bubble has a width equal 2*Radius*widthRatio
Epsilon = 1e-8


def slot(x):
    return np.heaviside(x+1, 0.5)-np.heaviside(x-1, 0.5)


def polynom(x):
    return -6*np.abs(x)**5+15*x**4-10*np.abs(x)**3+1


def profile(x):
    return slot(x)*polynom(x)


def density(x, y):
    assert x.shape == y.shape
    funcs = np.zeros((x.shape[0], 2))
    for i in range(2):
        x_ = (x-Centers[i][0])/Radius[i]
        y_ = (y-Centers[i][1])/Radius[i]
        r_ = np.sqrt(x_**2+y_**2)
        funcs[:, i] = nMain*profile(r_)
    return nBack+funcs.sum(axis=1)


def bx(x, y):
    assert x.shape == y.shape
    funcs = np.zeros((x.shape[0], 2))
    for i in range(2):
        x_ = (x-Centers[i][0])/Radius[i]
        y_ = (y-Centers[i][1])/Radius[i]
        r_ = np.sqrt(x_**2+y_**2)
        s_ = (r_-1)/widthRatio

        X_ = +y_
        Y_ = -x_
        R_ = np.clip(np.sqrt(X_**2+Y_**2), Epsilon, None)

        funcs[:, i] = profile(s_)*X_/R_
    return funcs.sum(axis=1)


def by(x, y):
    assert x.shape == y.shape
    funcs = np.zeros((x.shape[0], 2))
    for i in range(2):
        x_ = (x-Centers[i][0])/Radius[0]
        y_ = (y-Centers[i][1])/Radius[1]
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


def vth(x, y):
    return 1.


def Te(x, y):
    return 0.2

def Pe(x, y):
    return density(x, y)*Te(x, y)



def config():
    Simulation(
        time_step=0.004,
        final_time=.2,
        cells=(200, 100),
        dl=(0.2, 0.2),
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={"format": "phareh5",
                      "options": {"dir": "./coronas.run",
                                  "mode": "overwrite"}, },
        # restart_options={"dir": "checks",
        #                  "mode": "overwrite",
        #                  "timestamps": [.2], }
    )

    MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1,
                 "mass": 1.,
                 "density": density,
                 "vbulkx": v0,
                 "vbulky": v0,
                 "vbulkz": v0,
                 "vthx": vth,
                 "vthy": vth,
                 "vthz": vth,
                 "nbr_part_per_cell": 100,
                 "init": {"seed": 12}},
    )

    ElectronModel(closure="isothermal", Te=0.2)
    # ElectronModel(closure="polytropic", Pe=Pe, gamma=1.66)

    # sim = ph.global_vars.sim
    # dt = 100.*sim.time_step
    # nt = (sim.final_time)/dt
    # timestamps = dt * np.arange(nt+1)
    timestamps = [.0, .1, .2]

    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
        )

    for quantity in ["density", "charge_density", "bulkVelocity"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            )

    for quantity in ["density", "flux"]:
        FluidDiagnostics(quantity=quantity,
                         write_timestamps=timestamps,
                         population_name="protons")

    for quantity in ['domain']:  # , 'levelGhost', 'patchGhost']:
        ParticleDiagnostics(quantity=quantity,
                            write_timestamps=timestamps,
                            population_name="protons")


def main():

    config()
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()


if __name__ == "__main__":
    main()
