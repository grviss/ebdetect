;+
; NAME:
;	  EBDETECT_ARRAY_APPEND
;
; PURPOSE:
;	  Append values to an array (creating one if it doesn't yet exist)
;
; CATEGORY:
;   Array manipulation
;	
; CALLING SEQUENCE:
;   result = EBDETECT_ARRAY_APPEND(Inputarray, Appendix)
;
; INPUTS:
;   Inputarray  - Variable name of the array to be appended to
;   Appendix    - Value (or array) to be appended
;
; KEYWORD PARAMETERS:
;   DIMS  - Dimensionality of the array to be appended.
;
; OUTPUTS:
;   Array with value/array "Appendix" appended to the existing input array. If
;   input array does not exist, it will be initialised with value/array
;   "Appendix".
;
; OPTIONAL OUTPUTS:
;   None
;
; COMMON BLOCKS:
;   None
;
; SIDE EFFECTS:
;   None
;
; RESTRICTIONS:
;   None
;
; EXAMPLE:
;   IDL> FOR i=0,9 DO values = EBDETECT_ARRAY_APPEND(values, i)
;   IDL> HELP, values
;     VALUES    INT     = ARRAY[10]
;   IDL> PRINT, values
;     0   1   2   3   4   5   6   7   8   9
;
; MODIFICATION HISTORY:
;   2017 Dec 05 Gregal Vissers: Tailored and simplified version from own
;                               ARRAY_APPEND procedure 
;-
;
FUNCTION EBDETECT_ARRAY_APPEND, Inputarray, Appendix, DIMS=dims

	IF (N_PARAMS() LT 2) THEN BEGIN
		MESSAGE, 'Syntax: result = EBDETECT_ARRAY_APPEND(Inputarray, Appendix '+$
      '[, DIMS=dims])', /INFO
		RETURN, -1
	ENDIF

  ; If Inputarray is not undefined, append the Appendix, otherwise initialize
  ; Inputarry with Appendix as first value(s)
  IF NOT (SIZE(Inputarray, /TYPE) EQ 0) THEN BEGIN
    IF (N_ELEMENTS(DIMS) NE 1) THEN $
      result = [Inputarray, Appendix] $
    ELSE BEGIN
      CASE dims OF
        0:  result = [Inputarray, Appendix] 
        1:  result = [Inputarray, Appendix] 
        2:  result = [[Inputarray], [Appendix]] 
        3:  result = [[[Inputarray]], [[Appendix]]] 
        ELSE: BEGIN
                MESSAGE, 'ERROR: The DIMS keyword may only be supplied with '+$
                  'values between 0-3', /INFO
                STOP
              END
      ENDCASE
    ENDELSE
  ENDIF ELSE $
    result = Appendix

  RETURN, result

END
