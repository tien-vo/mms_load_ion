__all__ = ["tvo_mms_load_ion_dist"]

from .utils import get_data, interpol
from .tvo_mms_load_feeps import tvo_mms_load_feeps
from pyspedas.mms import mms_load_fgm, mms_load_edp, mms_load_feeps, mms_load_fpi
from pytplot import del_data, get_data as tp_get_data
import astropy.constants as c
import astropy.units as u
import numpy as np


def tvo_mms_load_ion_dist(trange, species="ion", probe="1", drate="brst"):

    fpi_drate = "fast" if drate == "srvy" else drate
    dtype = "dis" if species == "ion" else "des"
    q = c.si.e if species == "ion" else -c.si.e
    m = c.si.m_p if species == "ion" else c.si.m_e

    # Load FPI distribution files
    kw = dict(trange=trange, probe=probe, time_clip=True, latest_version=True)
    mms_load_fgm(data_rate=drate, varnames=[f"mms{probe}_fgm_b_gse_{drate}_l2"], **kw)
    mms_load_edp(data_rate=fpi_drate, datatype="scpot", **kw)
    mms_load_fpi(data_rate=fpi_drate, datatype=f"{dtype}-dist", get_support_data=True, center_measurement=True, **kw)

    # ---- Unpack data
    # Transpose FPI DF to (time, energy, theta, phi)
    t_fpi, f3d_fpi, _, _, _ = get_data(f"mms{probe}_{dtype}_dist_{fpi_drate}")
    f3d_fpi = f3d_fpi.transpose((0, 3, 2, 1)).astype("f8")
    _, W_fpi = get_data(f"mms{probe}_{dtype}_energy_{fpi_drate}")
    W_fpi = W_fpi[..., np.newaxis, np.newaxis].astype("f8") * np.ones(f3d_fpi.shape)
    theta_fpi = (
        tp_get_data(f"mms{probe}_{dtype}_theta_brst")[np.newaxis, np.newaxis, :, np.newaxis] * u.deg
    ).astype("f8") * np.ones(f3d_fpi.shape)
    _, phi_fpi = get_data(f"mms{probe}_{dtype}_phi_{fpi_drate}")
    phi_fpi = phi_fpi[:, np.newaxis, np.newaxis, :].astype("f8") * np.ones(f3d_fpi.shape)

    # Get unit angle
    t_fgm, B = get_data(f"mms{probe}_fgm_b_gse_brst_l2")
    B = interpol(B, t_fgm, t_fpi, avg=True)[:, 0:3]
    b = B / np.linalg.norm(B, axis=-1, keepdims=True)

    # Filter spacecraft potential
    t_edp, Vsc = get_data(f"mms{probe}_edp_scpot_brst_l2")
    Vsc = (q * interpol(Vsc, t_edp, t_fpi, avg=True)[:, np.newaxis, np.newaxis, np.newaxis]).to(u.eV)
    _Vsc = Vsc * np.ones(f3d_fpi.shape)
    f3d_fpi[W_fpi <= _Vsc] = np.nan

    # Load feeps
    feeps_data = tvo_mms_load_feeps(trange, species=species, probe=probe, drate=drate)
    f1d_omni_feeps = interpol(feeps_data["f_omni_avg"], feeps_data["t"].astype("datetime64[ns]"), t_fpi)
    W_feeps = feeps_data["f_omni_energy"]
    V_feeps = np.sqrt(2 * W_feeps / m).to(u.Unit("cm/s"))

    # Calculate f1d
    V_fpi = np.sqrt(2 * W_fpi / m).to(u.Unit("cm/s"))
    Omega = np.sum(np.sin(theta_fpi), axis=(2, 3)) * u.sr
    f1d_omni_fpi = (
        np.sum((V_fpi ** 4 / 2) * f3d_fpi * np.sin(theta_fpi), axis=(2, 3)) / Omega
    ).to(u.Unit("cm-2 s-1 sr-1"))

    return dict(
        t=t_fpi, b=b, B=B,
        f1d_omni_fpi=f1d_omni_fpi, f3d_fpi=f3d_fpi,
        V_fpi=V_fpi, W_fpi=W_fpi, phi_fpi=phi_fpi, theta_fpi=theta_fpi,
        V_feeps=V_feeps, W_feeps=W_feeps, f1d_omni_feeps=f1d_omni_feeps,
    )
