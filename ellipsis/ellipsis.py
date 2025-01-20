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


ellipsisCenter = [[40., 0.], [40., 20.]]
ellipsisAxis = [32., 8.]
nMain = 0.8
nBack = 0.2
widthRatio = 0.3
Epsilon = 1e-8


def slot(x):
    return np.heaviside(x+1, 0.5)-np.heaviside(x-1, 0.5)


def polynom(x):
    return -6*np.abs(x)**5+15*x**4-10*np.abs(x)**3+1


def profile(x):
    return slot(x)*polynom(x)


def density(x, y):
    assert x.shape == y.shape
    # funcs = np.zeros((x.shape[0], x.shape[1], 2))
    funcs = np.zeros((x.shape[0], 2))
    for i in range(2):
        x_ = (x-ellipsisCenter[i][0])/ellipsisAxis[0]
        y_ = (y-ellipsisCenter[i][1])/ellipsisAxis[1]
        r_ = np.sqrt(x_**2+y_**2)
        # funcs[:, :, i] = nMain*profile(r_)
        funcs[:, i] = nMain*profile(r_)
    # return funcs.sum(axis=2)
    return nBack+funcs.sum(axis=1)


def bx(x, y):
    assert x.shape == y.shape
    # funcs = np.zeros((x.shape[0], x.shape[1], 2))
    funcs = np.zeros((x.shape[0], 2))
    for i in range(2):
        x_ = (x-ellipsisCenter[i][0])/ellipsisAxis[0]
        y_ = (y-ellipsisCenter[i][1])/ellipsisAxis[1]
        r_ = np.sqrt(x_**2+y_**2)
        s_ = (r_-1)/widthRatio

        X_ = +y_/ellipsisAxis[1]
        Y_ = -x_/ellipsisAxis[0]
        R_ = np.clip(np.sqrt(X_**2+Y_**2), Epsilon, None)

        # normalisation to ensiure divB == 0
        Z_ = np.clip(ellipsisAxis[1]*R_/r_, Epsilon, None)

        # funcs[:, :, i] = profile(s_)*Z_*X_/R_
        funcs[:, i] = profile(s_)*Z_*X_/R_
    # return funcs.sum(axis=2)
    return funcs.sum(axis=1)


def by(x, y):
    assert x.shape == y.shape
    # funcs = np.zeros((x.shape[0], x.shape[1], 2))
    funcs = np.zeros((x.shape[0], 2))
    for i in range(2):
        x_ = (x-ellipsisCenter[i][0])/ellipsisAxis[0]
        y_ = (y-ellipsisCenter[i][1])/ellipsisAxis[1]
        r_ = np.sqrt(x_**2+y_**2)
        s_ = (r_-1)/widthRatio

        X_ = +y_/ellipsisAxis[1]
        Y_ = -x_/ellipsisAxis[0]
        R_ = np.clip(np.sqrt(X_**2+Y_**2), Epsilon, None)

        # normalisation to ensiure divB == 0
        Z_ = np.clip(ellipsisAxis[1]*R_/r_, Epsilon, None)

        # funcs[:, :, i] = profile(s_)*Z_*Y_/R_
        funcs[:, i] = profile(s_)*Z_*Y_/R_
    # return funcs.sum(axis=2)
    return funcs.sum(axis=1)


def bz(x, y):
    return 0.


def v0(x, y):
    return 0.


def vth(x, y):
    return 1.


def Te(x, y):
    return 0.2



def config():
    Simulation(
        time_step=0.004,
        final_time=.2,
        cells=(400, 100),
        dl=(0.2, 0.2),
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={"format": "phareh5",
                      "options": {"dir": "./ellipsis",
                                  "mode": "overwrite"}, },
        restart_options={"dir": "checks",
                         "mode": "overwrite",
                         "timestamps": [.2], }
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

    ElectronModel(closure="polytropic", Te=0.2, gamma=1.66)

    # sim = ph.global_vars.sim
    # dt = 100.*sim.time_step
    # nt = (sim.final_time)/dt
    # timestamps = dt * np.arange(nt+1)
    timestamps = [0., .1, .2]

    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
        )

    for quantity in ["density", "mass_density", "bulkVelocity"]:
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
