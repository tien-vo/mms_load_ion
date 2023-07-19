__all__ = ["get_data"]

from pytplot import get_data as tp_get_data
from collections import namedtuple
import astropy.units as u


def get_data(tname):
    tp_data = tp_get_data(tname, xarray=True)
    if tp_data is None:
        return None

    data = dict()
    data["t"] = tp_data.time.values.astype("datetime64[ns]")
    data["y"] = tp_data.values
    if hasattr(tp_data, "CDF"):
        data["y"] *= u.Unit(tp_data.CDF["VATT"]["UNITS"])

    for attr in ["v", "v1", "v2", "v3"]:
        if hasattr(tp_data, attr):
            data[attr] = tp_data[attr].values

    return namedtuple("tvo_get_data", data)(**data)
