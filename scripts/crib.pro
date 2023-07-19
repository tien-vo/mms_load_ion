tvo_mms_init

probe = '1'
drate = 'brst'
;drate = 'srvy'
;trange = ['2017-07-26T07:20:22', '2017-07-26T07:50:00']
trange = ['2017-07-26T07:28:00', '2017-07-26T07:29:30']

tvo_mms_load_ion_dist, trange=trange, probe=probe, drate=drate, error=error, /feeps
tvo_mms_load_ion_dist, trange=trange, probe=probe, drate=drate, error=error, suffix='_ext', /feeps, /extrapolate, $
    /keep_raw

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
