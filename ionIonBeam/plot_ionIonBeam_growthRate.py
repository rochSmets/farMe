
import os
import h5py
import numpy as np
from numpy.fft import fft
from scipy.optimize import curve_fit

from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import finest_part_data
from pyphare.pharesee.hierarchy import get_times_from_h5

import matplotlib.pyplot as plt


run_path = "/home/smets/shErpA/onGoing/ionIonBeam_refine_TeNull"
pop_name = ["main", "beam"]


xlim=(0, 33)
ylim=(-1.5, 5.5)


file = os.path.join(run_path, "EM_B.h5")
all_times = get_times_from_h5(file)
sparse_times = all_times[:1000:20]


r = Run(run_path)
first_mode = np.array([])


for time in sparse_times:
    B_hier = r.GetB(time, merged=True, interp="linear")

    by_interpolator, xyz_finest = B_hier["By"]
    bz_interpolator, xyz_finest = B_hier["Bz"]

    x = xyz_finest[0][:-1]

    by = by_interpolator(x)
    bz = by_interpolator(x)

    mode1 = np.absolute(fft(by-1j*bz)[1])
    first_mode = np.append(first_mode, mode1)


def croaCroa(x, a, b):
    return a*np.exp(np.multiply(b, x))


popt, pcov = curve_fit(croaCroa, sparse_times, first_mode, p0=[0.08, 0.09])
p1 = popt


damped_mode=first_mode*croaCroa(sparse_times, 1/p1[0], -p1[1])

omegas = np.fabs(fft(damped_mode).real)
omega = 0.5*(omegas[1:omegas.size//2].argmax()+1)*2*np.pi/sparse_times[-1]

fig, (ax1, ax2) = plt.subplots(2, 1)

ax1.set_title("Right Hand Resonant mode (Beam instability)")
ax1.stem(sparse_times, first_mode, linefmt='-k', basefmt=' ', use_line_collection=True)
ax1.plot(sparse_times, croaCroa(sparse_times, popt[0], popt[1]), color='r', linestyle='-', marker='')
ax1.text(0.04, 0.80, "From Gary et al., 1985 (ApJ : 10.1086/162797)", transform=ax1.transAxes)
ax1.set_ylabel("Most unstable mode")
ax1.set_title("Right Hand Resonant mode (Beam instability)")
ax1.text(0.30, 0.50, "gamma = {:5.3f}... expected 0.09".format(p1[1]), transform=ax1.transAxes)

ax2.plot(sparse_times, damped_mode, color='g', linestyle='', marker='o')
ax2.set_xlabel("Time")
ax2.set_ylabel("Real mode")
ax2.text(0.40, 0.20, "omega (real) = {:5.3f}... expected 0.19".format(omega), transform=ax2.transAxes)

fig.savefig("growthRate.pdf", bbox_inches='tight')
fig.show()

