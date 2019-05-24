;+
; NAME:
;	  EBDETECT_PARSE_CONFIGFILE
;
; PURPOSE:
;	  Parse line from configuration file for EBDETECT
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
;   2017 Oct 16 GV: Added check and parsing of 2D arrays as value_txt
;-
;

FUNCTION EBDETECT_PARSE_CONFIGFILE, Inputline
  
  IF (N_PARAMS() LT 1) THEN BEGIN
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
        ; Check whether it's a 2D array, if so reconstruct
        IF (STRMID(tmp_value,0,1) EQ '[') THEN BEGIN
          frstidx = 0L
          lastidx = STRPOS(tmp_value,']',0)
          prev = 0
          WHILE (frstidx[N_ELEMENTS(frstidx)-1] GE 0) DO BEGIN
            frstidx = [frstidx,STRPOS(tmp_value,'[',frstidx[prev]+1)]
            lastidx = [lastidx,STRPOS(tmp_value,']',lastidx[prev]+1)]
            prev += 1
          ENDWHILE
          frstidx = frstidx[0:N_ELEMENTS(frstidx)-2]
          lastidx = lastidx[0:N_ELEMENTS(lastidx)-2]
          value = STRCOMPRESS(STRSPLIT(STRMID(tmp_value,frstidx[0]+1,$
            lastidx[0]-frstidx[0]-1), ',', /EXTRACT), /REMOVE_ALL)
          FOR ii = 1,N_ELEMENTS(frstidx)-1 DO BEGIN
            loc_value =  $
              STRCOMPRESS(STRSPLIT(STRMID(tmp_value,frstidx[ii]+1,$
              lastidx[ii]-frstidx[ii]-1), ',', /EXTRACT), /REMOVE_ALL)
            value = [[value], [loc_value]]
          ENDFOR
        ENDIF ELSE $
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
