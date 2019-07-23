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
;   NWSUMS  - Number of wing blocks that EBDETECT thresholds over. Defaults to 1.
;
; OUTPUTS:
;   Structure with geometric and flux-weighted centroids, each a 2-element
;   array containing the respective x- and y-coordinates of the centroid
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   Inspiration from https://en.wikipedia.org/wiki/Centroid
;   2017 Aug 22 Gregal Vissers: First version
;   2017 Dec 4 GV:  Fixed error with optional (keyword) parameters
;   2019 Jul 23 GV: Added NWSUMS keyword fixing flux-weighted centroid for 2D
;                   SelDet.int
;-
FUNCTION EBDETECT_GET_CENTROID, SelDet, Mask, DIMS=dims, FLUX=flux, $
  NWSUMS=nwsums

  IF (N_PARAMS() LT 1) THEN BEGIN
    MESSAGE, 'Syntax: Result = EBDETECT_GET_CENTROID(SelDet [, Mask] '+$
      '[, DIMS=dims] [, /FLUX])', /INFO
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
		a += (sel_path_xy[0,i] * sel_path_xy[1,i+1] - sel_path_xy[0,i+1] * sel_path_xy[1,i])
		cx += (sel_path_xy[0,i] + sel_path_xy[0,i+1]) * (sel_path_xy[0,i] * sel_path_xy[1,i+1] - $
      sel_path_xy[0,i+1] * sel_path_xy[1,i])
		cy += (sel_path_xy[1,i] + sel_path_xy[1,i+1]) * (sel_path_xy[0,i] * sel_path_xy[1,i+1] - $
      sel_path_xy[0,i+1] * sel_path_xy[1,i])
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
    IF (nwsums GT 1) THEN BEGIN
      IF (npos GT 1) THEN $
        int_max = MAX(MEAN(SelDet.int, DIM=2, /NAN), wheremax) $
      ELSE $
        int_max = MAX(SelDet.int, /NAN, wheremax)
      int_tmp = SelDet.int[wheremax,*]
    ENDIF ELSE $
      int_tmp = SelDet.int
    cx_flux = TOTAL(xys[0,*] * int_tmp) / FLOAT(TOTAL(int_tmp))
    cy_flux = TOTAL(xys[1,*] * int_tmp) / FLOAT(TOTAL(int_tmp))
    xy_flux = [cx_flux, cy_flux]
  ENDIF ELSE xy_flux = !VALUES.F_NAN

  ; Construct and return result x- and y-coordinates
  result = {xy_geom:xy_geom, xy_flux:xy_flux}

  RETURN, result

END
