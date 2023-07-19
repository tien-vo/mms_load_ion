pro tvo_mms_init, verbose=verbose, debug=debug, local=local, auth=auth

; ---- Defaults
if undefined(verbose) then verbose = 0
if undefined(auth) then auth = './mms_auth_info.sav'
!quiet = ~verbose

; ---- Spedas initialization
; Debug
if defined(debug) then dprint, setdebug=debug
; Data directory
if defined(local) then if is_string(local) then data_dir = local else data_dir = './cdf/'
; Log into LASP SDC
_ = mms_login_lasp(login_info=auth)
mms_init, local_data_dir=data_dir
mms_set_verbose, verbose

; ---- Plots
; Fix colors
tvlct, red, green, blue, /get
red[6]      = 225
green[4]    = 135
red[5]      = 255
green[5]    = 200
tvlct, red, green, blue
!p.charsize     = 1.5
; Define circular markers
a = findgen(17) * (!pi * 2 / 16.)
usersym, cos(a), sin(a), /fill

end
