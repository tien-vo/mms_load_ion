from mms_load_ion import tvo_mms_load_ion_dist
import astropy.constants as c
import astropy.units as u
import mpl_utils as mu
import numpy as np
import h5py as h5


probe = "1"
drate = "brst"
trange = ["2017-07-26T07:28:00", "2017-07-26T07:29:30"]
data = tvo_mms_load_ion_dist(trange, probe=probe, drate=drate)

# Unpack
t = data["t"].astype("datetime64[ns]")
B = data["B"].to(u.nT)
V_fpi = data["V_fpi"].to(u.Unit("cm/s"))
V_feeps = data["V_feeps"].to(u.Unit("cm/s"))[np.newaxis, :]
W_fpi = data["W_fpi"].to(u.eV)
W_feeps = data["W_feeps"].to(u.eV)[np.newaxis, :]
theta_fpi = data["theta_fpi"].to(u.deg)
phi_fpi = data["phi_fpi"].to(u.deg)
f1d_omni_fpi = data["f1d_omni_fpi"].to(u.Unit("cm-2 s-1 sr-1"))
f1d_omni_feeps = data["f1d_omni_feeps"].to(u.Unit("cm-2 s-1 sr-1"))
f3d_fpi = data["f3d_fpi"].to(u.Unit("s3 cm-6"))

fig, axes = mu.plt.subplots(3, 1, sharex=True, figsize=(12, 10))

mu.add_colorbar(ax := axes[0]).remove()
ax.plot(t, B[:, 0], "-b")
ax.plot(t, B[:, 1], "-g")
ax.plot(t, B[:, 2], "-r")
ax.plot(t, np.linalg.norm(B, axis=1), "-k")

cax = mu.add_colorbar(ax := axes[1])
tg = (t[:, np.newaxis].astype("f8") * np.ones(f1d_omni_feeps.shape)).astype("datetime64[ns]")
Wg = W_feeps * np.ones(f1d_omni_feeps.shape)
im = ax.pcolormesh(tg, Wg.value, f1d_omni_feeps.value, cmap="jet", norm=mu.mplc.LogNorm(1e3, 1e6))
fig.colorbar(im, cax=cax)
ax.set_yscale("log")
ax.set_ylim(Wg.value.min(), Wg.value.max())

cax = mu.add_colorbar(ax := axes[2])
tg = (t[:, np.newaxis].astype("f8") * np.ones(f1d_omni_fpi.shape)).astype("datetime64[ns]")
Wg = W_fpi[:, :, 0, 0]
im = ax.pcolormesh(tg, Wg.value, f1d_omni_fpi.value, cmap="jet", norm=mu.mplc.LogNorm(1e3, 1e6))
fig.colorbar(im, cax=cax)
ax.set_yscale("log")
ax.set_ylim(Wg.value.min(), Wg.value.max())

mu.plt.show()


with h5.File("data.h5", "w") as h5f:
    h5f.create_dataset(f"/raw/t", data=t.astype("f8"))
    h5f.create_dataset(f"/raw/B", data=B.value.astype("f8"))
    h5f.create_dataset(f"/raw/V_fpi", data=V_fpi.value.astype("f8"))
    h5f.create_dataset(f"/raw/W_fpi", data=W_fpi.value.astype("f8"))
    h5f.create_dataset(f"/raw/theta_fpi", data=theta_fpi.value.astype("f8"))
    h5f.create_dataset(f"/raw/phi_fpi", data=phi_fpi.value.astype("f8"))
    h5f.create_dataset(f"/raw/V_feeps", data=V_feeps.value.astype("f8"))
    h5f.create_dataset(f"/raw/W_feeps", data=W_feeps.value.astype("f8"))
    h5f.create_dataset(f"/raw/f1d_omni_fpi", data=f1d_omni_fpi.value.astype("f8"))
    h5f.create_dataset(f"/raw/f1d_omni_feeps", data=f1d_omni_feeps.value.astype("f8"))
    h5f.create_dataset(f"/raw/f3d_fpi", data=f3d_fpi.value.astype("f8"))

