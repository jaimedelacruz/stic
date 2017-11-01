;;
;; Function to read named variables from a netcdf-4 file
;; Written by J. de la Cruz Rodriguez (ISP-SU 2017)
;;
function readNetcdf, file, varname, exists=exists, printout=printout
  
  if(n_elements(varname) eq 0) then varname = ' '
  
  id = NCDF_OPEN(file)
  finq = NCDF_INQUIRE( id )
  
  defined = 0
  vnames = strarr(finq.nvars)
  for k=0, finq.nvars-1 do begin
     var = NCDF_VARINQ( id, k )
     vnames[k] = var.name
     if(var.name eq varname) then defined = 1
  endfor


  if(keyword_set(exists)) then return, defined
  
  if(defined eq 0) then begin
     if(~keyword_set(printout)) then print, 'readNetcdf: ERROR, variable [', varname,'] not defined. The file contains:'
     for k=0, finq.nvars-1 do print, '[',string(k, format='(I3)'),'] '+vnames[k]
     if(keyword_set(printout)) then return, 0
     
     if(finq.nvars eq 1) then begin
        print, 'readNetcdf: falling back to ['+vnames[0],']'
        varname = vnames[0]
        NCDF_VARGET, id, varname, m
        ncdf_close, id
        return, m
     endif else begin
        sel = 0
        read, sel, prompt='Select var id: '
        varname = vnames[sel]
        NCDF_VARGET, id, varname, m
        ncdf_close, id
        return, m
     endelse
     
     
     return, 0
  endif

  print, 'loadvar: loading [', varname,'] from ', file
  NCDF_VARGET, id, varname, m
  ncdf_close, id

  
  return, temporary(m)
end
