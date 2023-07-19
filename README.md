# About
Contains IDL and Python routines to load MMS 3D ion distribution function.

# Usage

Run
```
source setup_env.sh
```

to put `src/idl` into IDL path. The main routine is `src/idl/tvo_mms_load_ion_dist.pro`, which loads and preprocesses
FPI (and FEEPS) ion distribution function. A Python package (`src/mms_load_ion`) is included for testing purposes.
