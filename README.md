# About
Contains IDL and Python routines to load MMS 3D ion distribution function. Works for both fast survey and burst data.

# Usage

Run
```
source setup_env.sh
```

to put `src/idl` into IDL path. The main routine is `src/idl/tvo_mms_load_ion_dist.pro`, which loads and preprocesses
FPI (and FEEPS) ion distribution function. A Python package (`src/mms_load_ion`) is included for testing purposes.

Find below the docstring for `tvo_mms_load_ion_dist.pro` and example crib at `scripts/crib.pro` for usage.
```
PROCEDURE: tvo_mms_load_ion_dist, trange=trange, probe=probe, drate=drate, /feeps

PURPOSE:
      Load and process MMS ion distribution function from FPI (and FEEPS if toggled).

KEYWORDS:
      trange: 2-element array of time range. (Fmt: 'YYYY-MM-DD/hh:mm:ss')
      probe: Which MMS spacecrafts? (Valid values: '1' '2' '3' '4')
      drate: Instrument data rate (Valid values: 'srvy' 'brst'; 'srvy' is converted to 'fast' for FPI)
      feeps: Toggle to combine FEEPS distribution into final product
      keep_raw: Toggle to keep raw (unprocessed) variables.
      extrapolate: Toggle to extrapolate the energy gap.
      suffix: Save processed variables with suffix.
      error: 1 = Error during processing, 0 = No error.
```

The processed 3D and omni-directional distribution functions will be saved into tplot variables
`mms?_ion_dist_3d_?` and `mms?_ion_dist_omni_?`. By default, the 3D DF is reshaped into (time, energy, theta, phi).
Toggle `/no_reshape` to keep it in FPI's original shape (time, phi, theta, energy).
