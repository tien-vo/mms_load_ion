;+
; FUNCTION logspace(x1, x2, N)
;
; PURPOSE:
;       Create a real array of length N, containing equidistant values
;   between x1 and x2 inclusively in log space, mimicking numpy's logspace
;   function.
;
; INPUTS:
;       X1 (required): The starting point of the array.
;       X2 (required): The starting point of the array.
;       N: The array length (default to 100).
;
; KEYWORDS:
;       BASE: The base of the log scale. Default: 10
;
; OUTPUT:
;       ARR: Equidistant array in log space from x1 to x2. 
;-
function logspace, x1, x2, N, base=base
    if ~keyword_set(N) then N = 100
    if ~keyword_set(base) then base = 10d

    arr = base^(x1 + dindgen(N) * (x2 - x1) / (N - 1))

    return, arr
end
