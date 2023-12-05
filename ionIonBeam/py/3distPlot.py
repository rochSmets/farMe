
import os
import h5py
import numpy as np
from numpy.fft import fft
from scipy.optimize import curve_fit

from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import finest_part_data
from pyphare.pharesee.hierarchy import get_times_from_h5

import matplotlib.pyplot as plt


run_path = "/home/smets/sherpa/onGoing/ionIonBeam/ionIonBeam_refine"
pop_name = ["main", "beam"]



allParts1 = Run(run_path).GetParticles(time=0, pop_name=pop_name)
allParts2 = Run(run_path).GetParticles(time=64, pop_name=pop_name)
allParts3 = Run(run_path).GetParticles(time=85.4, pop_name=pop_name)

fig,(ax1,ax2,ax3) = plt.subplots(nrows=3, sharex='row', figsize=(10,8))

xlim=(0, 33)
ylim=(-1.5, 5.5)

allParts1.dist_plot(axis=("x", "Vx"),
                   ax=ax1,
                   finest=True,
                   color_scale="log",
                   color_max=0.1,
                   xlim=xlim,
                   ylim=ylim,
                   vmin=ylim[0],
                   vmax=ylim[1],
                   xlabel="",
                   ylabel= ""
                  )
ax1.vlines(10, ylim[0], ylim[1], color="k")
ax1.vlines(14, ylim[0], ylim[1], color="r")
ax1.vlines(18, ylim[0], ylim[1], color="r")
ax1.vlines(22, ylim[0], ylim[1], color="k")

allParts2.dist_plot(axis=("x", "Vx"),
                   ax=ax2,
                   finest=True,
                   color_scale="log",
                   color_max=0.1,
                   xlim=xlim,
                   ylim=ylim,
                   vmin=ylim[0],
                   vmax=ylim[1],
                   xlabel="",
                   ylabel= "Vx - Velocity"
                  )
ax2.vlines(10, ylim[0], ylim[1], color="k")
ax2.vlines(14, ylim[0], ylim[1], color="r")
ax2.vlines(18, ylim[0], ylim[1], color="r")
ax2.vlines(22, ylim[0], ylim[1], color="k")

allParts3.dist_plot(axis=("x", "Vx"),
                   ax=ax3,
                   finest=True,
                   color_scale="log",
                   color_max=0.1,
                   xlim=xlim,
                   ylim=ylim,
                   vmin=ylim[0],
                   vmax=ylim[1],
                   xlabel="X - Position",
                   ylabel= ""
                  )
ax3.vlines(10, ylim[0], ylim[1], color="k")
ax3.vlines(14, ylim[0], ylim[1], color="r")
ax3.vlines(18, ylim[0], ylim[1], color="r")
ax3.vlines(22, ylim[0], ylim[1], color="k")

for ax in (ax1,ax2):
   ax.tick_params(labelbottom=False)

fig.savefig("3distPlot.pdf", bbox_inches='tight')

