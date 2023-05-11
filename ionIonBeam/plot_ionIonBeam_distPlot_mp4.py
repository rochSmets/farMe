
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

file = os.path.join(run_path, "EM_B.h5")
all_times = get_times_from_h5(file)

xlim=(0, 33)
ylim=(-1.5, 5.5)

for it,t in enumerate(all_times):

    print("{0:6d} / {1:6d}".format(it, all_times.size))
    fig, ax = plt.subplots(figsize=(10,4))

    allParts = Run(run_path).GetParticles(time=t, pop_name=pop_name)

    allParts.dist_plot(axis=("x", "Vx"),
                       ax=ax,
                       finest=True,
                       sigma=(1,1),
                       color_scale="log",
                       cmap='viridis',
                       color_max=0.1,
                       xlim=xlim,
                       ylim=ylim,
                       vmin=ylim[0],
                       vmax=xlim[1],
                       title="2 stream [ RhR ] - Time : {:8.4f}".format(t),
                       xlabel="X - Position",
                       ylabel= "Vx - Velocity"
                      )

    ax.vlines(10, ylim[0], ylim[1], color="k")
    ax.vlines(14, ylim[0], ylim[1], color="r")
    ax.vlines(18, ylim[0], ylim[1], color="r")
    ax.vlines(22, ylim[0], ylim[1], color="k")

    fig.savefig("mp4/t_{:04d}.png".format(it))
    fig.clf()

    plt.close(fig)

# ffmpeg -r 24 -y -i t_%04d.png -b:v 18000k -pix_fmt yuv420p dist.mp4

