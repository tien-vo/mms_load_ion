;+
; PROCEDURE: tvo_mms_load_ion_dist, trange=trange, probe=probe, drate=drate, /feeps
;
; PURPOSE:
;       Load and process MMS ion distribution function from FPI (and FEEPS if toggled).
;
; KEYWORDS:
;       trange: 2-element array of time range. (Fmt: 'YYYY-MM-DD/hh:mm:ss')
;       probe: Which MMS spacecrafts? (Valid values: '1' '2' '3' '4')
;       drate: Instrument data rate (Valid values: 'srvy' 'brst'; 'srvy' is converted to 'fast' for FPI)
;       feeps: Toggle to combine FEEPS distribution into final product
;       keep_raw: Toggle to keep raw (unprocessed) variables.
;       extrapolate: Toggle to extrapolate the energy gap.
;       suffix: Save processed variables with suffix.
;       error: 1 = Error during processing, 0 = No error.
;−
pro tvo_mms_load_ion_dist, trange=trange, probe=probe, drate=drate, feeps=feeps, keep_raw=keep_raw, $
    extrapolate=extrapolate, suffix=suffix, error=error

compile_opt idl2
tvo_mms_init, debug=-1

; ---- Defaults and constants
error = 0
SQRTEV_TO_CMS = 1.384d+6
SQCMS_TO_EV = 5.22d-13
N_ext = 5  ; Number of extrapolated data points
F_thresh = 4e4  ; cm-2 s-1 sr-1
if undefined(probe) then probe = '1'
if undefined(drate) then drate = 'brst'
if undefined(suffix) then suffix = ''
if undefined(trange) then begin
    dprint, dlevel=-1, 'Provide time range!'
    error = 1
    return
end
if drate eq 'srvy' then fpi_drate = 'fast' else fpi_drate = 'brst'

; -- Names
f3d_name = string(probe, fpi_drate, format='mms%s_dis_dist_%s_raw')
Vsc_name = string(probe, fpi_drate, format='mms%s_edp_scpot_%s_l2_raw')
f1d_feeps_name = string(probe, drate, format='mms%s_epd_feeps_%s_l2_ion_intensity_omni_spin_raw')
out_f3d_name = string(probe, fpi_drate, format='mms%s_ion_dist_3d_%s')
out_f1d_omni_name = string(probe, fpi_drate, format='mms%s_ion_dist_omni_%s')
; ----

; ---- Load data and preprocess
has_data = 1

dprint, dlevel=-1, 'Loading FPI data...'
mms_load_fpi, $
    trange=trange, $
    probes=probe, $
    data_rate=fpi_drate, $
    datatype='dis-dist', $
    varformat='*_dist_*', $
    suffix='_raw', $
    tplotnames=loaded_vars, $
    /latest_version, $
    /time_clip, $
    /center_measurement
if ~keyword_set(loaded_vars) then begin
    dprint, dlevel=-1, string(probe, fpi_drate, format='Failed loading MMS%s %s FPI data.')
    error = 1
    return
endif else undefine, loaded_vars
has_data = has_data && tdexists(f3d_name, trange[0], trange[1])

dprint, dlevel=-1, 'Loading EDP data...'
mms_load_edp, $
    trange=trange, $
    probes=probe, $
    data_rate=fpi_drate, $
    datatype='scpot', $
    varformat='*scpot*', $
    suffix='_raw', $
    tplotnames=loaded_vars, $
    /latest_version, $
    /time_clip
if ~keyword_set(loaded_vars) then begin
    dprint, dlevel=-1, string(probe, fpi_drate, format='Failed loading MMS%s %s EDP data.')
    error = 1
    return
endif else undefine, loaded_vars
has_data = has_data && tdexists(Vsc_name, trange[0], trange[1])

if defined(feeps) then begin
    mms_load_feeps, $
        trange=trange, $
        probes=probe, $
        data_rate=drate, $
        datatype='ion', $
        tplotnames=loaded_vars, $
        suffix='_raw', $
        /latest_version, $
        /time_clip
    if ~keyword_set(loaded_vars) then begin
        dprint, dlevel=-1, string(probe, drate, format='Failed loading MMS%s %s FEEPS data.')
        error = 1
        return
    endif else undefine, loaded_vars
    has_data = has_data && tdexists(f1d_feeps_name, trange[0], trange[1])
endif
if ~has_data then begin
    dprint, dlevel=-1, 'Data not available!'
    error = 1
    return
endif

; -- Transpose DF shape to (N = time, Nw = energy, Nt = theta, Np = phi)
get_data, f3d_name, data=f3d_data, dlim=dlim
f3d = transpose(f3d_data.y, [0, 3, 2, 1])  ; s3 cm-6
s = size(f3d)
N = s[1]
Nw = s[2]
Nt = s[3]
Np = s[4]

; -- Broadcast DF support to be of the same shape as f3d
dprint, dlevel=-1, 'Broadcasting 3D DF support...'
t_fpi = f3d_data.x
W = dblarr(N, Nw, Nt, Np)       ; eV
theta = dblarr(N, Nw, Nt, Np)   ; rad
phi = dblarr(N, Nw, Nt, Np)     ; rad
for it = 0, Nt - 1 do for ip = 0, Np - 1 do W[*, *, it, ip] = f3d_data.v3
for n = 0, N - 1 do for iw = 0, Nw - 1 do for ip = 0, Np - 1 do theta[n, iw, *, ip] = (!dpi / 180d) * f3d_data.v2
if fpi_drate eq 'brst' then $
    for iw = 0, Nw - 1 do for it = 0, Nt - 1 do phi[*, iw, it, *] = (!dpi / 180d) * f3d_data.v1 else $
    for n = 0, N - 1 do for iw = 0, Nw - 1 do for it = 0, Nt - 1 do phi[n, iw, it, *] = (!dpi / 180d) * f3d_data.v1
V = SQRTEV_to_CMS * sqrt(W)     ; cm/s

; -- Get spacecraft potential and broadcast to W shape
dprint, dlevel=-1, 'Masking spacecraft potential...'
_Vsc = tsample(Vsc_name, trange, time=t_edp)
dt_fpi = mean(t_fpi[1:-1] - t_fpi[0:-2])
dt_edp = mean(t_edp[1:-1] - t_edp[0:-2])
_Vsc = interpol(movavg(_Vsc, dt_fpi / dt_edp), t_edp, t_fpi)
Vsc = dblarr(N, Nw, Nt, Np)
for iw = 0, Nw - 1 do for it = 0, Nt - 1 do for ip = 0, Np - 1 do Vsc[*, iw, it, ip] = _Vsc
idx = where(W le Vsc)
if idx[0] ne -1 then f3d[idx] = !values.d_nan

; -- Unpack FEEPS data
if defined(feeps) then begin
    tinterpol, f1d_feeps_name, f3d_name, suffix='_i', /nan_extrapolate
    f1d_omni_feeps = tsample(f1d_feeps_name + '_i', trange, values=W_feeps)
    W_feeps = 1d+3 * W_feeps[1:*]  ; eV
    f1d_omni_feeps = 1d-3 * f1d_omni_feeps[*, 1:*]  ; cm-2 s-1 sr-1 eV-1
    Nw_feeps = (size(f1d_omni_feeps))[2]
    V_feeps = SQRTEV_TO_CMS * sqrt(W_feeps)
endif

if undefined(keep_raw) then remove_tvar, string(probe, format='mms%s_*raw*')
; ----

; ---- Filter 3D DF based on reduced 1D DF
dprint, dlevel=-1, 'Filtering background radiation...'

; 1D DF from raw f3d
Omega = total(total(sin(theta), 4, /nan, /double), 3, /nan, /double)  ; sr
f1d_omni = total(total(0.5d * V^4 * f3d * sin(theta), 4, /nan, /double), 3, /nan, /double) / Omega  ; cm-2 s-1 sr-1

; Bg subtraction
f1d_bg = dblarr(N, Nw)
for n = 0, N - 1 do f1d_bg[n, *] = mean((f1d_omni[n, sort(f1d_omni[n, *])])[0:5], /nan, /double)
f1d_omni -= f1d_bg

; Filter zeroes
idx = where(f1d_omni le F_thresh)
if idx[0] ne -1 then for it = 0, Nt - 1 do for ip = 0, Np - 1 do begin
    _f3d = f3d[*, *, it, ip]
    _f3d[idx] = !values.d_nan
    f3d[*, *, it, ip] = _f3d
endfor
f1d_omni = total(total(0.5d * V^4 * f3d * sin(theta), 4, /nan, /double), 3, /nan, /double) / Omega  ; cm-2 s-1 sr-1

if defined(feeps) then begin
    dprint, dlevel=-1, 'Combining FEEPS...'

    ; Placeholder
    W_fpi = W
    V_fpi = V
    f3d_fpi = f3d
    f1d_omni_fpi = f1d_omni
    theta_fpi = theta
    phi_fpi = phi

    ; Combined velocity
    V = dblarr(N, Nw + N_ext + Nw_feeps, Nt, Np)
    V[*, 0:Nw - 1, *, *] = V_fpi
    for n = 0, N - 1 do for it = 0, Nt - 1 do for ip = 0, Np - 1 do begin
        V_ext = logspace(alog10(V_fpi[n, -1, it, ip]), alog10(V_feeps[0]), N_ext + 2)
        V[n, Nw:Nw + N_ext - 1, it, ip] = V_ext[1:-2]
        V[n, Nw + N_ext:Nw + N_ext + Nw_feeps - 1, it, ip] = V_feeps
    endfor
    W = SQCMS_TO_EV * V^2

    ; Combined f3d
    f3d = dblarr(N, Nw + N_ext + Nw_feeps, Nt, Np)
    f3d[*, 0:Nw - 1, *, *] = f3d_fpi
    for n = 0, N - 1 do for it = 0, Nt - 1 do for ip = 0, Np - 1 do begin
        ; Apply linear extrapolation
        if defined(extrapolate) then begin
            _f_feeps = f1d_omni_feeps[n, 0] * W_feeps[0]
            _f_fpi = f1d_omni_fpi[n, -1]
            _W_feeps = W_feeps[0]
            _W_fpi = W_fpi[n, -1, it, ip]
            W_ext = logspace(alog10(_W_fpi), alog10(_W_feeps), N_ext + 2)
            V_ext = SQRTEV_TO_CMS * sqrt(W_ext)
            slope = alog10(_f_feeps / _f_fpi) / alog10(_W_feeps / _W_fpi)
            f_ext = 10d^(alog10(_f_fpi) - alog10(_W_fpi) * slope) * W_ext^slope
            f3d[n, Nw:Nw + N_ext - 1, it, ip] = (2d / V_ext^4) * f_ext[1:-2]
        endif
        f3d[n, Nw + N_ext:Nw + N_ext + Nw_feeps - 1, it, ip] = (2d / V_feeps^4) * f1d_omni_feeps[n, *] * W_feeps
    endfor

    ; Recalculate f1d_omni
    theta = dblarr(N, Nw + N_ext + Nw_feeps, Nt, Np)
    phi = dblarr(N, Nw + N_ext + Nw_feeps, Nt, Np)
    for iw = 0, Nw + N_ext + Nw_feeps - 1 do begin
        theta[*, iw, *, *] = theta_fpi[*, 0, *, *]
        phi[*, iw, *, *] = phi_fpi[*, 0, *, *]
    endfor
    Omega = total(total(sin(theta), 4, /nan, /double), 3, /nan, /double)  ; sr
    f1d_omni = total(total(0.5d * V^4 * f3d * sin(theta), 4, /nan, /double), 3, /nan, /double) / Omega  ; cm-2 s-1 sr-1

endif

; Remove zeroes
idx = where(f3d eq 0)
if idx[0] ne -1 then f3d[idx] = !values.d_nan
idx = where(f1d_omni eq 0)
if idx[0] ne -1 then f1d_omni[idx] = !values.d_nan
; ----

; ---- Store data
if defined(extrapolate) then sfx = ' with EXT' else sfx = ''
W_fpi_thresh = mean(W[*, Nw - 1, *, *], /nan)
W_feeps_thresh = mean(W[*, Nw + N_ext - 1, *, *], /nan)
store_data, out_f1d_omni_name + suffix, data={x: t_fpi, y: f1d_omni, v: reform(W[*, *, 0, 0])}
options, out_f1d_omni_name + suffix, $
    spec=1, $
    ylog=1, $
    zlog=1, $
    yrange=[min(W, /nan), max(W, /nan)], $
    ystyle=1, $
    zrange=[min(f1d_omni, /nan), max(f1d_omni, /nan)], $
    zstyle=1, $
    ytitle=string(probe, sfx, format='MMS%s ion OMNI %s'), $
    ysubtitle='(eV)', $
    ztitle='(cm!U-2!N s!U-1!N sr!U-1!N)', $
    databar={yval: [W_fpi_thresh, W_feeps_thresh], color: [1, 1], linestyle: 2, thick: 2}
if fpi_drate eq 'brst' then phi_save = phi[*, 0, 0, *] else phi_save = phi[0, 0, 0, *]
store_data, out_f3d_name + suffix, $
    data={x: t_fpi, y: f3d, v1: reform(W[*, *, 0, 0]), v2: reform(theta[0, 0, *, 0]), v3: reform(phi_save)}, $
    dlim=dlim
; ----

end
