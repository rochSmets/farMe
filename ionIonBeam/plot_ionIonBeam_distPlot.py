
import os
import h5py
import numpy as np
from numpy.fft import fft
from scipy.optimize import curve_fit

from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import finest_part_data
from pyphare.pharesee.hierarchy import get_times_from_h5

import matplotlib.pyplot as plt


run_path = "/home/smets/shErpA/onGoing/ionIonBeam_refine"
pop_name = ["main", "beam"]


# ___ 1 plot in phase space
time = 84.0

allParts = Run(run_path).GetParticles(time=time, pop_name=pop_name)

particles = finest_part_data(allParts)

fig, ax = plt.subplots(figsize=(10,4))

xlim=(0, 33)
ylim=(-1.5, 5.5)

allParts.dist_plot(axis=("x", "Vx"),
                   ax=ax,
                   finest=True,
                   #gaussian_filter_sigma=(1,1),
                   median_filter_size=(2,2),
                   color_scale="log",
                   cmap="viridis",
                   color_max=0.1,
                   xlim=xlim,
                   ylim=ylim,
                   vmin=ylim[0],
                   vmax=ylim[1],
                   title="2 stream [ {} ]".format(run_path),
                   xlabel="X - Position",
                   ylabel= "Vx - Velocity"
                  )

ax.vlines(10, ylim[0], ylim[1], color="k")
ax.vlines(14, ylim[0], ylim[1], color="r")
ax.vlines(18, ylim[0], ylim[1], color="r")
ax.vlines(22, ylim[0], ylim[1], color="k")

#plt.show()
fig.savefig("distPlot.pdf", bbox_inches='tight')

