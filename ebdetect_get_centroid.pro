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
;   Result = EBDETECT_GET_CENTROID(Path_xy, Npath)
;
; INPUTS:
;   Path_xy - Contour edge of a detection (e.g., return value of PATH_XY keyword
;             to IDL's CONTOUR procedure)
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
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
;-
FUNCTION EBDETECT_GET_CENTROID, Path_xy

  id_string = '$Id$'
  IF (N_PARAMS() LT 1) THEN BEGIN
    MESSAGE, id_string, /INFO
    MESSAGE, 'Syntax: Result = EBDETECT_GET_CENTROID(Path_xy)', /INFO
    RETURN, -1
  ENDIF

	a = 0.
	cx = 0.
	cy = 0.
  npath = N_ELEMENTS(Path_xy[0,*])
	FOR i=0,npath-1 DO BEGIN
		a += (path_xy[0,i] * path_xy[1,i+1] - path_xy[0,i+1] * path_xy[1,i])
		cx += (path_xy[0,i] + path_xy[0,i+1]) * (path_xy[0,i] * path_xy[1,i+1] - $
      path_xy[0,i+1] * path_xy[1,i])
		cy += (path_xy[1,i] + path_xy[1,i+1]) * (path_xy[0,i] * path_xy[1,i+1] - $
      path_xy[0,i+1] * path_xy[1,i])
	ENDFOR
	a /= 2.
	cx /= (6.*a)
	cy /= (6.*a)
  result = [cx,cy]

  ; Return x- and y-coordinates
  RETURN, result

END
