
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

xlim=(0, 33)
ylim=(-1.5, 5.5)

# ___ then scatter plot
plt.figure(figsize=(10,4))
plt.scatter(particles["main_domain"].x,
            particles["main_domain"].v[:,0],
            s=particles["main_domain"].weights*100)
plt.scatter(particles["beam_domain"].x,
            particles["beam_domain"].v[:,0],
            s=particles["beam_domain"].weights*10000)
plt.axvspan( 0, 10,facecolor="0.8", alpha=0.2, linewidth=0.0)
plt.axvspan(14,18,facecolor="0.8", alpha=0.2, linewidth=0.0)
plt.axvspan(22,33,facecolor="0.8", alpha=0.2, linewidth=0.0)
plt.xlim((0,33))
plt.ylim((-2,+7))
plt.xlabel("X")
plt.ylabel("Vx")
plt.title("Time : {:8.4f}".format(time))

plt.show()
plt.savefig("scatterPlot.pdf")

