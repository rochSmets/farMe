import numpy as np
from pyphare.pharesee.run import Run


def available_times(RUN_DIR, file="EM_B.h5"):
    """dump times available for {file} in {RUN_DIR}"""
    import h5py
    with h5py.File(f"{RUN_DIR}/{file}") as f:
        return np.array(sorted(float(k) for k in f["t"].keys()))


def interior(pdata):
    """interior (no ghost) values and coordinates of a single-patch field"""
    gx, gy = pdata.ghosts_nbr
    sx = slice(gx, -gx if gx else None)
    sy = slice(gy, -gy if gy else None)
    return pdata.x[sx], pdata.y[sy], np.asarray(pdata.dataset)[sx, sy]


def scalar(hier, t):
    """the (unique) scalar field of a hierarchy at time t -> x, y, 2D array"""
    # level 0 may be split in several patches (mpirun) : stitch them together
    parts = []
    for patch in hier.levels(t)[0].patches:
        pds = patch.patch_datas
        parts.append(interior(pds[list(pds.keys())[0]]))
    if len(parts) == 1:
        return parts[0]
    dx = parts[0][0][1] - parts[0][0][0]
    dy = parts[0][1][1] - parts[0][1][0]
    x0 = min(p[0][0] for p in parts)
    y0 = min(p[1][0] for p in parts)
    nx = int(round((max(p[0][-1] for p in parts) - x0)/dx)) + 1
    ny = int(round((max(p[1][-1] for p in parts) - y0)/dy)) + 1
    a = np.full((nx, ny), np.nan)
    for px, py, pa in parts:   # shared primal nodes are simply overwritten
        i0 = int(round((px[0]-x0)/dx))
        j0 = int(round((py[0]-y0)/dy))
        a[i0:i0+len(px), j0:j0+len(py)] = pa
    return x0 + dx*np.arange(nx), y0 + dy*np.arange(ny), a


def component(vhier, comp, t):
    """one component ('x', 'y' or 'z') of a vector field hierarchy"""
    return scalar(getattr(vhier, comp), t)


def az_from_b(bx, by, dx, dy):
    """Az (on the primal x primal nodes) such that Bx = dAz/dy, By = -dAz/dx.

    Periodic FFT solve, consistent with the Yee staggering of PHARE:
      Bx[i,j] sits at (x_i,      y_j+dy/2) = (Az[i,j+1]-Az[i,j])/dy
      By[i,j] sits at (x_i+dx/2, y_j     ) = -(Az[i+1,j]-Az[i,j])/dx
    so that the discrete curl is inverted exactly (Coulomb gauge, <Az> = 0).
    bx, by : periodic arrays of identical shape (nx, ny), no duplicated node.
    """
    nx, ny = bx.shape
    kx = 2*np.pi*np.fft.fftfreq(nx, dx)[:, None]
    ky = 2*np.pi*np.fft.fftfreq(ny, dy)[None, :]
    Dx = (np.exp(1j*kx*dx) - 1)/dx      # forward difference operators
    Dy = (np.exp(1j*ky*dy) - 1)/dy
    den = abs(Dx)**2 + abs(Dy)**2
    den[0, 0] = 1.                      # mean mode : Az mean is arbitrary -> 0
    a_k = (np.conj(Dy)*np.fft.fft2(bx) - np.conj(Dx)*np.fft.fft2(by))/den
    a_k[0, 0] = 0.
    return np.fft.ifft2(a_k).real


def vector_potential(RUN_DIR, t, run=None):
    """[x, y, Az] at (nearest available) time t, on level 0 (primal nodes)."""
    ts = available_times(RUN_DIR)
    t = ts[np.argmin(abs(ts - t))]
    run = run or Run(RUN_DIR)
    B = run.GetB(t)
    x, _, bx = component(B, "x", t)     # x primal, y dual
    _, y, by = component(B, "y", t)     # x dual,   y primal
    # primal directions carry one extra (duplicated) node on a periodic domain
    nx, ny = min(bx.shape[0], by.shape[0]), min(bx.shape[1], by.shape[1])
    dx, dy = x[1]-x[0], y[1]-y[0]
    az = az_from_b(bx[:nx, :ny], by[:nx, :ny], dx, dy)
    return [x[:nx], y[:ny], az]


def density_movie(RUN_DIR, filename="density.mp4", times=None, fps=10, dpi=150,
                  vmin=0.2, vmax=1.0, levels=8, cmap="cool"):
    """mp4 of the total ion density (pharesee plot) with the contours of Az
    (in-plane magnetic field lines), one frame per available dump time.

    times : optional subset of times (default : all the dumps of the run)
    levels: number of Az contours; the Az range of the first frame is used for
            all frames so that contours are comparable in time.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FFMpegWriter

    run = Run(RUN_DIR)
    times = available_times(RUN_DIR) if times is None else np.asarray(times)
    fig, ax = plt.subplots(figsize=(6, 6))
    writer = FFMpegWriter(fps=fps, codec="libx264",
                          extra_args=["-pix_fmt", "yuv420p"])
    az_levels = None
    with writer.saving(fig, filename, dpi):
        for t in times:
            fig.clf()
            ax = fig.add_subplot(111)
            run.GetNi(t).plot(qty="value", ax=ax, cmap=cmap,
                              vmin=vmin, vmax=vmax,
                              title=f"Total ion density, t={t:g}",
                              xlabel="X ($d_i$)", ylabel="Y ($d_i$)")
            x, y, az = vector_potential(RUN_DIR, t, run=run)
            if az_levels is None:
                az_levels = np.linspace(az.min(), az.max(), levels)
            ax.contour(x, y, az.T, levels=az_levels, colors="k",
                       linestyles="solid", linewidths=0.5)
            # ax.contour(x, y, az.T, levels=levels, colors="k",
            #            linestyles="solid", linewidths=0.5)
            writer.grab_frame()
            print(f"t={t:g} done", flush=True)
    plt.close(fig)
    return filename




