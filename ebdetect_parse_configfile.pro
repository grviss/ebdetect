;+
; NAME:
;	  EBDETECT_PARSE_CONFIGFILE
;
; PURPOSE:
;	  Handles parsing the lines of the configuration file for EBDETECT
;
; CATEGORY:
;   Data analysis
;	
; CALLING SEQUENCE:
;	  result = EBDETECT_PARSE_CONFIGFILE(Inputline)
;
; INPUTS:
;	  Inputline  - input line from a file containing keyword information (the
;                equal sign ('=') is assumed to be used as delimiter between
;                keyword name and its value).
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   Structure with the keyword 'field' name and its 'value' extracted from the
;   line. 
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
;   result = EBDETECT_PARSE_CONFIGFILE('field_txt = value_txt')
;
;   result is then a structure with result.field = 'field_txt' and 
;   result.value = 'value_txt'
;
; MODIFICATION HISTORY:
;   2016 Nov 30 Gregal Vissers: First version
;-
;

FUNCTION EBDETECT_PARSE_CONFIGFILE, Inputline
  
  id_string = '$Id$'
  IF (N_PARAMS() LT 1) THEN BEGIN
    MESSAGE, id_string, /INFO
    MESSAGE, 'Syntax: result = EBDETECT_PARSE_CONFIGFILE(Inputline)', $
      /INFO
    RETURN, -1
  ENDIF 

  ; Set defaults in case of empty line
  field = ''
  value = 0B

  IF (N_ELEMENTS(Inputline) EQ 1) THEN BEGIN
    ; Check that the line is not commented out or a white line
    first_char = STRMID(Inputline, 0, 1)
    IF ((first_char NE '#') AND (STRCOMPRESS(first_char) NE '')) THEN BEGIN
      splitline = STRCOMPRESS(STRSPLIT(Inputline, '=', /EXTRACT), /REMOVE_ALL)
      field = splitline[0]
      tmp_value = splitline[1]
      ; Check whether value is actually an array of values, if so reconstruct
      IF (STRMID(tmp_value,0,1) EQ '[') THEN BEGIN
        tmp_value = STRMID(tmp_value,1,STRLEN(tmp_value)-2)
        value = STRCOMPRESS(STRSPLIT(tmp_value, ',', /EXTRACT), /REMOVE_ALL)
      ENDIF ELSE $
        value = tmp_value
    ENDIF 
  ENDIF ELSE BEGIN
    MESSAGE,'Input line is empty. Returning.',/INFO
  ENDELSE

  ; Create result structure and return
  result = {field:field, value:value}

  RETURN, result
END
