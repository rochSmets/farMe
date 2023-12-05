
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
fig.savefig("dist1.pdf", bbox_inches='tight')

plt.clf()

# ___ now 3 plots also in phase space
allParts1 = Run(run_path).GetParticles(time=0, pop_name=pop_name)
allParts2 = Run(run_path).GetParticles(time=64, pop_name=pop_name)
allParts3 = Run(run_path).GetParticles(time=85.4, pop_name=pop_name)

fig,(ax1,ax2,ax3) = plt.subplots(nrows=3, sharex='row', figsize=(10,8))

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

fig.savefig("dist3.pdf", bbox_inches='tight')



plt.clf()

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

plt.savefig("scatter.pdf")


plt.clf()

# ___produce several dist_plot to make a mp4
file = os.path.join(run_path, "EM_B.h5")
all_times = get_times_from_h5(file)
sparse_times = all_times[:800:20]


def getMode(t, m):
    return np.absolute(fft(t)[m])


r = Run(run_path)
first_mode = np.array([])


for time in sparse_times:
    B_hier = r.GetB(time, merged=True, interp="linear")

    by_interpolator, xyz_finest = B_hier["By"]
    bz_interpolator, xyz_finest = B_hier["Bz"]

    x = xyz_finest[0]

    by = by_interpolator(x)
    bz = by_interpolator(x)

    mode = getMode(by-1j*bz, 1)
    first_mode = np.append(first_mode, mode)


def croaCroa(x, a, b):
    return a*np.exp(np.multiply(b, x))


popt, pcov = curve_fit(croaCroa, sparse_times, first_mode, p0=[0.1, 0.1])

plt.clf()
plt.figure(figsize=(6,4))

plt.stem(sparse_times, first_mode, linefmt='-k', basefmt=' ', use_line_collection=True)
#plt.plot(sparse_times, first_mode, color='k', linestyle=' ', marker='o')
plt.plot(sparse_times, croaCroa(sparse_times, popt[0], popt[1]), color='red')

plt.xlabel("Time")
plt.ylabel("First mode")#FFT ($B_y - \imath B_z$) : Mode m=1")
#plt.text(0, 8, "$\gamma$ = {:6.4f}".format(popt[1]))
print("gamma : {}".format(popt[1]))

plt.savefig("gamma.png", bbox_inches='tight')
plt.show()















