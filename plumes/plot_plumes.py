#!/usr/bin/env python3
# %% [markdown]
# # Plume run : 2D maps and radial profiles
#
# Written in "percent" format : each `# %%` line starts a cell, so the file can
# be run cell by cell (VS Code, Spyder, ipython `%run`) or converted into a
# notebook later (`jupytext --to ipynb plot_plumes.py`).
#
# Run with PHARE in the path, for instance
#   export PYTHONPATH=$B:~/far/PHARE:~/far/PHARE/pyphare:$B/pybindlibs
# with B=~/far/builds/release/master

# %% imports
import numpy as np
import matplotlib.pyplot as plt
from pyphare.pharesee.run import Run

# %% settings : edit here
RUN_DIR = "plumes.run"
CENTER = (10., 10.)      # plume center, used for radial profiles
TIMES = [0., 2., 5., 10., 20.]   # times to display (nearest available dump is used)
CMAPS = {"B": "viridis", "Bz": "RdBu_r", "n": "magma", "E": "viridis"}

run = Run(RUN_DIR)


# %% helpers
def available_times():
    """dump times available in the run"""
    import h5py
    with h5py.File(f"{RUN_DIR}/EM_B.h5") as f:
        return np.array(sorted(float(k) for k in f["t"].keys()))


def interior(pdata):
    """interior (no ghost) values and coordinates of a single-patch field"""
    gx, gy = pdata.ghosts_nbr
    sx = slice(gx, -gx if gx else None)
    sy = slice(gy, -gy if gy else None)
    return pdata.x[sx], pdata.y[sy], np.asarray(pdata.dataset)[sx, sy]


def scalar(hier, t):
    """the (unique) scalar field of a hierarchy at time t -> x, y, 2D array"""
    pds = hier.levels(t)[0].patches[0].patch_datas
    return interior(pds[list(pds.keys())[0]])


def component(vhier, comp, t):
    """one component ('x', 'y' or 'z') of a vector field hierarchy"""
    return scalar(getattr(vhier, comp), t)


def load(t):
    """all the fields used below at (available) time t, as a dict name -> (x, y, 2D)"""
    ts = available_times()
    t = ts[np.argmin(abs(ts - t))]
    B, E = run.GetB(t), run.GetE(t)
    f = {}
    x, y, bx = component(B, "x", t)
    _, _, by = component(B, "y", t)
    _, _, bz = component(B, "z", t)
    f["B"] = (x, y, np.sqrt(bx**2+by**2+bz**2))
    f["Bz"] = (x, y, bz)
    ex, ey, ez = (component(E, c, t) for c in "xyz")
    f["E"] = (ex[0], ex[1], np.sqrt(ex[2]**2+ey[2]**2+ez[2]**2))
    f["n"] = scalar(run.GetNi(t), t)
    f["n_main"] = scalar(run.GetN(t, pop_name="main"), t)
    f["n_back"] = scalar(run.GetN(t, pop_name="back"), t)
    return t, f


def show(ax, field, title, cmap, vmin=None, vmax=None):
    x, y, a = field
    # pcolormesh: data is [ix, iy] -> transpose for imshow-like orientation
    im = ax.pcolormesh(x, y, a.T, cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="auto")
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


# %% 2D maps : one row per time, one column per field
# common color scales across times so that the evolution is readable
SCALES = {"B": (0, 1), "Bz": (-0.1, 0.1), "n": (0, 1), "n_main": (0, 0.8),
          "n_back": (0.1, 0.3), "E": (0, 1.5)}
FIELDS = ["B", "n", "n_main", "n_back"]   # edit the list : see SCALES

fig, axes = plt.subplots(len(TIMES), len(FIELDS),
                         figsize=(4.2*len(FIELDS), 3.8*len(TIMES)),
                         squeeze=False)
for row, t in enumerate(TIMES):
    t, f = load(t)
    for col, name in enumerate(FIELDS):
        vmin, vmax = SCALES.get(name, (None, None))
        cmap = CMAPS.get(name, CMAPS.get(name.split("_")[0], "viridis"))
        show(axes[row, col], f[name], f"{name}  t={t:g}", cmap, vmin, vmax)
fig.tight_layout()
plt.show()


# %% radial profiles around CENTER (azimuthal average)
def radial_profile(field, center=CENTER, dr=0.4, rmax=None):
    x, y, a = field
    X, Y = np.meshgrid(x, y, indexing="ij")
    r = np.hypot(X-center[0], Y-center[1])
    rmax = rmax or r.max()
    edges = np.arange(0, rmax+dr, dr)
    idx = np.digitize(r.ravel(), edges) - 1
    prof = np.array([a.ravel()[idx == i].mean() if np.any(idx == i) else np.nan
                     for i in range(len(edges)-1)])
    return 0.5*(edges[1:]+edges[:-1]), prof


fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for t in TIMES:
    t, f = load(t)
    for ax, name in zip(axes, ["B", "n_main", "n_back"]):
        r, p = radial_profile(f[name], rmax=10.)
        ax.plot(r, p, label=f"t={t:g}")
        ax.set_title(f"{name} (azimuthal mean)")
        ax.set_xlabel("r")
for ax in axes:
    ax.grid(alpha=.3)
axes[0].legend()
fig.tight_layout()
plt.show()


# %% time history : max|E|, max|B|, min density
ts = available_times()
hist = []
for t in ts[::2]:
    t, f = load(t)
    hist.append((t, f["E"][2].max(), f["B"][2].max(), f["n"][2].min()))
hist = np.array(hist)
fig, axes = plt.subplots(1, 3, figsize=(15, 3.5))
for ax, col, name in zip(axes, (1, 2, 3), ("max |E|", "max |B|", "min n")):
    ax.plot(hist[:, 0], hist[:, col])
    ax.set_title(name)
    ax.set_xlabel("t")
    ax.grid(alpha=.3)
fig.tight_layout()
plt.show()
