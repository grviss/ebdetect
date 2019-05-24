function lp_get, filename, indx, header=header, printheader=printheader,$
 verbose=verbose
;
; read La Palma data stored as IDL assoc file
; this is a one-image version of lp_read reading only one frame (indx) 
; of the series
;
; First 512 bytes are assumed to contain header info containing:
;    datatype=2 (integer), dims=2, nx=2027, ny=2042
;
; see lp_write for writing right header style
;
 IF n_params() LT 2 THEN BEGIN
    message, /info, 'image = lp_get(filename, indx [, header=header])'
    retall
 ENDIF

 IF n_elements(printheader) EQ 0 THEN printheader=0
 IF keyword_set(verbose) THEN printheader=1

 lp_header, filename, header=header, datatype=datatype, $
            dims=dims, nx=nx, ny=ny, nt=nt, endian=endian_file

if ((byte(1L, 0, 1))[0] eq 1) then endian = 'l' else endian='b'
if(n_elements(endian_file) eq 0) then endian_file=endian
if(datatype gt 1) and (endian ne endian_file) then swap_endian=1 else swap_endian=0
 openr, lur, filename, /get_lun, swap_endian=swap_endian
 if printheader then message, /info, header

 ; read actual data

 case datatype of
     1: begin 
         rec = assoc(lur, bytarr(nx,ny), 512)
         image = rec[indx]
     end
     2: begin
         rec = assoc(lur, intarr(nx,ny), 512)
         image = rec[indx]
     end
     3: begin
         rec = assoc(lur, lonarr(nx,ny), 512)
         image = rec[indx]
     end
     4: begin
         rec = assoc(lur, fltarr(nx,ny), 512)
         image = rec[indx]
     end
     else: begin
         message, /info, 'datatype not supported '
         print, ' datatype = ', datatype
         free_lun, lur
         return, 0
     end
 end
 free_lun, lur
 return, image

end
