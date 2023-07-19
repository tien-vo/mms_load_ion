from astropy.convolution import convolve, Gaussian2DKernel
import astropy.constants as c
import astropy.units as u
import mpl_utils as mu
import numpy as np
import h5py as h5


tmin = np.datetime64("2017-07-26T07:28:40", "ns").astype("f8")
tmax = np.datetime64("2017-07-26T07:28:50", "ns").astype("f8")
with h5.File("data.h5") as h5f:
    t = h5f["/raw/t"][:]
    B = h5f["/raw/B"][:]
    b = B / np.linalg.norm(B, axis=1, keepdims=True)
    f3d_fpi = h5f["/raw/f3d_fpi"][:] * u.Unit("s3 cm-6")
    f1d_fpi = h5f["/raw/f1d_omni_fpi"][:] * u.Unit("cm-2 s-1 sr-1")
    f1d_feeps = h5f["/raw/f1d_omni_feeps"][:] * u.Unit("cm-2 s-1 sr-1")
    V_fpi = h5f["/raw/V_fpi"][:] * u.Unit("cm/s")
    theta_fpi = h5f["/raw/theta_fpi"][:] * u.deg
    phi_fpi = h5f["/raw/phi_fpi"][:] * u.deg
    W_fpi = h5f["/raw/W_fpi"][:] * u.eV
    W_feeps = h5f["/raw/W_feeps"][:] * u.eV * np.ones(f1d_feeps.shape)
    V_feeps = h5f["/raw/V_feeps"][:] * u.Unit("cm/s") * np.ones(f1d_feeps.shape)
    cond = np.where((tmin <= t) & (t <= tmax))
    print(t)
    t = t[cond]
    b = b[cond]
    f3d_fpi = f3d_fpi[cond]
    f1d_fpi = f1d_fpi[cond]
    f1d_feeps = f1d_feeps[cond]
    V_fpi = V_fpi[cond]
    theta_fpi = theta_fpi[cond]
    phi_fpi = phi_fpi[cond]
    W_fpi = W_fpi[cond]
    W_feeps = W_feeps[cond]
    V_feeps = V_feeps[cond]

F_thresh = 4e4 * u.Unit("cm-2 s-1 sr-1")

# Estimate background and filter f3d
f1d_fpi_sorted = np.take_along_axis(f1d_fpi, np.argsort(f1d_fpi, axis=1), axis=1)
f1d_fpi = f1d_fpi - np.nanmean(f1d_fpi_sorted[:, :5], axis=1)[:, np.newaxis]
f3d_fpi[f1d_fpi <= F_thresh, ...] = np.nan

# Create combined data
N_ext = 5
N, Nw, Nt, Np = f3d_fpi.shape
Nw_feeps = W_feeps.shape[1]

V = np.zeros((N, Nw + N_ext + Nw_feeps, Nt, Np)) * u.Unit("cm/s")
V[:, 0:Nw, :, :] = V_fpi
V[:, Nw + N_ext:Nw + N_ext + Nw_feeps, :, :] = V_feeps[..., np.newaxis, np.newaxis]
for n in range(N):
    V_ext = np.logspace(np.log10(V_fpi[n, -1, 0, 0].value), np.log10(V_feeps[n, 0].value), N_ext + 2) * u.Unit("cm/s")
    V[n, Nw:Nw + N_ext, :, :] = V_ext[1:-1, np.newaxis, np.newaxis]

f3d = np.zeros((N, Nw + N_ext + Nw_feeps, Nt, Np)) * u.Unit("s3 cm-6")
f3d[:, 0:Nw, :, :] = f3d_fpi
f3d[:, Nw + N_ext:Nw + N_ext + Nw_feeps, :, :] = ((2 / V_feeps ** 4) * f1d_feeps * u.sr)[..., np.newaxis, np.newaxis]

for n in range(N):
    _f_feeps = f1d_feeps[n, 0].value
    _f_fpi = f1d_fpi[n, -1].value
    _W_feeps = W_feeps[n, 0].value
    _W_fpi = W_fpi[n, -1, 0, 0].value
    W_ext = np.logspace(np.log10(_W_fpi), np.log10(_W_feeps), N_ext + 2) * u.eV
    V_ext = np.sqrt(2 * W_ext / c.si.m_p).to(u.Unit("cm/s"))
    slope = np.log10(_f_feeps / _f_fpi) / np.log10(_W_feeps / _W_fpi)
    f_ext = (
        np.power(10.0, (np.log10(_f_fpi) - np.log10(_W_fpi) * slope)) * (W_ext.value ** slope) *
        u.Unit("cm-2 s-1")
    )
    f3d[n, Nw:Nw + N_ext, :, :] = (2 / V_ext ** 4 * f_ext)[1:-1].to(u.Unit("s3 cm-6"))[:, np.newaxis, np.newaxis]


theta = theta_fpi[0, 0, :, 0][np.newaxis, np.newaxis, :, np.newaxis] * np.ones(V.shape)
phi = phi_fpi[:, 0, 0, :][:, np.newaxis, np.newaxis, :] * np.ones(V.shape)
Omega = np.nansum(np.sin(theta), axis=(2, 3)) * u.sr
f1d = (np.nansum((V ** 4 / 2) * f3d * np.sin(theta), axis=(2, 3)) / Omega).to(u.Unit("cm-2 s-1 sr-1"))
W = (0.5 * c.si.m_p * V ** 2).to(u.eV)
tg = t[:, np.newaxis] * np.ones((N, Nw + N_ext + Nw_feeps))

# DBCS
V = V[..., np.newaxis] * np.array([
    np.sin(theta) * np.cos(phi),
    np.sin(theta) * np.sin(phi),
    np.cos(theta),
]).transpose((1, 2, 3, 4, 0))
Vx = V[..., 0]
Vy = V[..., 1]
Vz = V[..., 2]

# FAC
e3 = b[:, np.newaxis, np.newaxis, np.newaxis, :]
zhat = np.zeros(e3.shape)
zhat[..., 2] = 1
e1 = zhat - e3 * np.einsum("...i,...i", zhat, e3)[..., np.newaxis]
e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)
e2 = np.cross(e3, e1, axis=-1)
V1 = np.einsum("...i,...i", V, e1)
V2 = np.einsum("...i,...i", V, e2)
V3 = np.einsum("...i,...i", V, e3)

Vmag = np.sqrt(V1 ** 2 + V2 ** 2 + V3 ** 2)
alpha = np.degrees(np.arccos(V3 / Vmag))
f3d[f3d == 0] = np.nan
cond = ~np.isnan(f3d)
F3d = (Vmag ** 4 / 2 * f3d).to(u.Unit("cm-2 s-1"))

W_bins = np.logspace(np.log10(5e0), np.log10(1e6), 40) * u.eV
A_bins = np.arange(0, 181, 5) * u.deg
Ag, Wg = np.meshgrid(A_bins[:-1], W_bins[:-1], indexing="ij")
bins = (A_bins.value, W_bins.value)
H = np.histogram2d(alpha[cond].value, W[cond].value, bins=bins)[0]
Hf = np.histogram2d(alpha[cond].value, W[cond].value, bins=bins, weights=np.log10(F3d[cond].value))[0]
H[H == 0] = np.nan
f = np.power(10.0, Hf / H)


fig, ax = mu.plt.subplots(1, 1)

cax = mu.add_colorbar(ax)
im = ax.pcolormesh(tg.astype("datetime64[ns]"), W[:, :, 0, 0].value, f1d.value, norm=mu.mplc.LogNorm(1e3, 1e6))
fig.colorbar(im, cax=cax)

ax.set_facecolor("silver")
ax.set_yscale("log")
mu.format_datetime_axis(ax)
fig.tight_layout()


kernel = Gaussian2DKernel(1)
fig, axes = mu.plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)
kw = dict(cmap="jet", norm=mu.mplc.LogNorm(1e5, 1e7))

cax = mu.add_colorbar(axes[1])
im = axes[0].pcolormesh(Ag.value, Wg.value, f, **kw)
f[np.isnan(f)] = np.nanmin(f)
im = axes[1].pcolormesh(Ag.value, Wg.value, convolve(f, kernel), **kw)
fig.colorbar(im, cax=cax)

axes[0].set_ylabel("eV")
for (i, ax) in enumerate(axes):
    ax.set_facecolor("silver")
    ax.set_yscale("log")
    ax.set_xlim(Ag.value.min(), Ag.value.max())
    ax.set_ylim(Wg.value.min(), Wg.value.max())
    ax.set_xlabel("deg")

fig.tight_layout(w_pad=0.05)
mu.plt.show()

