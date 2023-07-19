__all__ = ["sampling_period", "move_avg", "interpol"]

from astropy.convolution import (
    Gaussian1DKernel,
    Gaussian2DKernel,
    CustomKernel,
    Box1DKernel,
    convolve,
)
from tvolib.utils import sampling_period
from scipy.interpolate import interp1d
from functools import reduce
import astropy.units as u
import numpy as np


def sampling_period(t):
    r""" Calculate the sampling period of a time array. """

    if np.issubdtype(t.dtype, np.datetime64):
        t = t.astype("datetime64[ns]")
        return (np.diff(t).mean().astype("timedelta64[ns]").astype("f8") * u.ns).to(u.s)
    else:
        raise NotImplementedError


def move_avg(x, w, smooth=False):
    r""" Moving average of a signal `x` with window sizes `w` """

    if not isinstance(w, (tuple, list)):
        raise NotImplementedError("w must be a tuple or list.")
    if (ndim := x.ndim) != len(w := np.array(w)):
        raise NotImplementedError("Length of w-tuple must be the same as x dimensions.")

    if ndim == 1:
        if smooth:
            kernel = Gaussian1DKernel(*(w / 6))
        else:
            kernel = Box1DKernel(*w)
    elif ndim == 2:
        if smooth:
            kernel = Gaussian2DKernel(*(w / 6), theta=np.pi / 2)
        else:
            kernel = CustomKernel(reduce(np.outer, (Box1DKernel(_w) for _w in w)).reshape(w))
    else:
        raise NotImplementedError("Dimensionality is constrained to 2.")

    return convolve(x, kernel, fill_value=np.nan, boundary="fill")


def interpol(y, x, xout, avg=False, smooth=False, kind="linear"):

    if not (np.issubdtype(x.dtype, np.datetime64) and np.issubdtype(xout.dtype, np.datetime64)):
        raise NotImplementedError("Time arrays must be datetime64 dtype.")

    w = (sampling_period(xout) / sampling_period(x)).decompose()
    w = w if smooth else max(np.int64(w), 1)
    kw = dict(kind=kind, bounds_error=False, fill_value=np.nan)

    x = x.astype("f8")
    xout = xout.astype("f8")
    if len(y.shape) == 1:
        yout = interp1d(x, move_avg(y, (w,), smooth=smooth) if avg else y, **kw)(xout)
    elif len(y.shape) == 2:
        yout = np.empty((len(xout), N := y.shape[1]))
        for i in range(N):
            yout[:, i] = interp1d(x, move_avg(y[:, i], (w,), smooth=smooth) if avg else y[:, i], **kw)(xout)
    else:
        raise NotImplementedError("Dimensionality is constrained to 2.")

    if isinstance(y, u.Quantity):
        yout *= u.Unit(y.unit)

    return yout
