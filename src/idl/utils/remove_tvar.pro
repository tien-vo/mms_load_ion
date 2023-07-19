;+
; PROCEDURE: remove_tvar, '*'
;
; PURPOSE:
;       Remove tplot variables from memory.
;
; INPUT:
;       varnames: Names to remove (accept wildcards)
;
; KEYWORD:
;       all: Toggle to remove everything.
;−
pro remove_tvar, varnames, all=all

compile_opt idl2

if defined(all) then tplot_names, names=names, /silent else tplot_names, varnames, names=names, /silent
if undefined(names) then return
if (size(names))[0] eq 0 then if names eq '' then return

store_data, delete=names

end
