;+
; NAME:
;	  EBDETECT_GET_CENTROID
;
; PURPOSE:
;	  Get centroid of detection contour
;
; CATEGORY:
;   Data analysis
;	
; CALLING SEQUENCE:
;   Result = EBDETECT_GET_CENTROID(Det)
;
; INPUTS:
;   SelDet  - EBDETECT detection structure (e.g., *(*sel_detection[0]).det[0])
;
; OPTIONAL INPUTS:
;   Mask    - 2D byte mask with selected pixels set to 1. Must be set if DIMS is
;             not set.
;
; KEYWORD PARAMETERS:
;   DIMS    - XY-dimensions of the full field of view. Defaults to not set and
;             is ignored if Mask is set.
;   FLUX    - Compute the flux-weighted average centroid. Defaults to not set.
;
; OUTPUTS:
;   2-element array containing the x- and y-coordinates of the centroid
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   Inspiration from https://en.wikipedia.org/wiki/Centroid
;   2017 Aug 22 Gregal Vissers: First version
;   2017 Dec 4 GV: Fixed error with optional (keyword) parameters
;-
FUNCTION EBDETECT_GET_CENTROID, SelDet, Mask, DIMS=dims, FLUX=flux

  id_string = '$Id$'
  IF (N_PARAMS() LT 1) THEN BEGIN
    MESSAGE, id_string, /INFO
    MESSAGE, 'Syntax: Result = EBDETECT_GET_CENTROID(Path_xy)', /INFO
    RETURN, -1
  ENDIF
 
  IF (N_ELEMENTS(MASK) LT 2) THEN BEGIN
    IF (N_ELEMENTS(DIMS) LT 2) THEN BEGIN
      MESSAGE, 'ERROR: DIMS keyword must be set if Mask is not provided!', /INFO
      RETURN,-1
    ENDIF ELSE BEGIN
      mask = BYTARR(dims[0], dims[1])
      mask[SelDet.pos] = 1B
    ENDELSE
  ENDIF
  CONTOUR, Mask, PATH_XY=path_xy, /PATH_DATA_COORDS, PATH_INFO=path_info
  sel_path_xy = FLTARR(2,path_info[0].n+1)
  sel_path_xy[*,0:((path_info[0]).n-1)] = path_xy[*,0:((path_info[0]).n-1)]
  sel_path_xy[*,path_info[0].n] = path_xy[*,0]

	a = 0.
	cx = 0.
	cy = 0.
	FOR i=0,(path_info[0]).n-1 DO BEGIN
		a += (path_xy[0,i] * path_xy[1,i+1] - path_xy[0,i+1] * path_xy[1,i])
		cx += (path_xy[0,i] + path_xy[0,i+1]) * (path_xy[0,i] * path_xy[1,i+1] - $
      path_xy[0,i+1] * path_xy[1,i])
		cy += (path_xy[1,i] + path_xy[1,i+1]) * (path_xy[0,i] * path_xy[1,i+1] - $
      path_xy[0,i+1] * path_xy[1,i])
	ENDFOR
	a /= 2.
	cx /= (6.*a)
	cy /= (6.*a)
  xy_geom = [cx, cy]

  IF KEYWORD_SET(FLUX) THEN BEGIN
    npos = N_ELEMENTS(SelDet.pos)
    xys = LONARR(2,npos)
    FOR j=0,npos-1 DO $
      xys[*,j] = ARRAY_INDICES(Mask, SelDet.pos[j])
    cx_flux = TOTAL(xys[0,*] * SelDet.int) / FLOAT(TOTAL(SelDet.int))
    cy_flux = TOTAL(xys[1,*] * SelDet.int) / FLOAT(TOTAL(SelDet.int))
    xy_flux = [cx_flux, cy_flux]
  ENDIF

  ; Construct and return result x- and y-coordinates
  result = {xy_geom:xy_geom, xy_flux:xy_flux}

  RETURN, result

END
