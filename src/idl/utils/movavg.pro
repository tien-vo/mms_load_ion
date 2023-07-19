function movavg, x, w, nan=nan
;+
; FUNCTION movavg(x, w)
;
; PURPOSE:
;       Calculate the moving average of a 1D array.
;
; INPUTS:
;       x (required): A 1D array.
;       w: The width over which the average is performed, default to 1
;
; OUTPUT:
;       y: <x>
;-

    if ~defined(w) then w = 1
    if w le 1 then return, x else return, smooth(x, w, /edge_truncate, nan=nan)

end
