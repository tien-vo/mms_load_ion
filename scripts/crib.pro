tvo_mms_init

probe = '1'
drate = 'brst'
;drate = 'srvy'
;trange = ['2017-07-26T07:20:22', '2017-07-26T07:50:00']
trange = ['2017-07-26T07:28:00', '2017-07-26T07:29:30']

tvo_mms_load_ion_dist, trange=trange, probe=probe, drate=drate, error=error, /feeps
tvo_mms_load_ion_dist, trange=trange, probe=probe, drate=drate, error=error, suffix='_ext', /feeps, /extrapolate, $
    /keep_raw

; Calculate reduced angle spectra for fun (f3d in DBCS, need to regrid to get PAD)
get_data, string(probe, drate, format='mms%s_ion_dist_3d_%s_ext'), data=f3d_data
t = f3d_data.x
f3d = f3d_data.y ; s3 cm-6
theta = f3d_data.v2 * (180d / !dpi)  ; deg
phi = f3d_data.v3 * (180d / !dpi) ; deg
f1d_T = 1d+15 * mean(mean(f3d, dimension=4, /nan), dimension=2, /nan)  ; s3 km-3 cm-3
f1d_P = 1d+15 * mean(mean(f3d, dimension=2, /nan), dimension=2, /nan)  ; s3 km-3 cm-3
store_data, 'theta_distribution', data={x: t, y: f1d_T, v: theta}
options, 'theta_distribution', spec=1, zlog=1, yrange=[0, 180], ystyle=1, ytitle='theta', ysubtitle='(deg)', $
    ztitle='(s!U3!N km!U-3!N cm!U-3!N)', zrange=[1d-17, 1d-8], zstyle=1
store_data, 'phi_distribution', data={x: t, y: f1d_P, v: phi}
options, 'phi_distribution', spec=1, zlog=1, yrange=[0, 360], ystyle=1, ytitle='phi', ysubtitle='(deg)', $
    ztitle='(s!U3!N km!U-3!N cm!U-3!N)', zrange=[1d-17, 1d-8], zstyle=1

; For comparison
mms_load_fgm, trange=trange, probes=probe, data_rate=drate, /time_clip, /no_split_vars, /latest_version
mms_load_fpi, trange=trange, probes=probe, data_rate=drate, datatype='dis-moms', $
    /latest_version, /time_clip, /center_measurement

zlim, string(probe, drate, format='mms%s_dis_energyspectr_omni_%s'), [1e3, 1e6]
zlim, string(probe, drate, format='mms%s_ion_dist_omni_%s'), [1e3, 1e6]
zlim, string(probe, drate, format='mms%s_ion_dist_omni_%s_ext'), [1e3, 1e6]
window, 0, xsize=1200, ysize=1000
tplot, [string(probe, drate, format='mms%s_fgm_b_gse_%s_l2'), $
        string(probe, drate, format='mms%s_epd_feeps_%s_l2_ion_intensity_omni_raw'), $
        string(probe, drate, format='mms%s_dis_energyspectr_omni_%s'), $
        string(probe, drate, format='mms%s_ion_dist_omni_%s'), $
        string(probe, drate, format='mms%s_ion_dist_omni_%s_ext')], window=0
tplot_apply_databar

