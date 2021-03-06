;+
; NAME:
;	  EBDETECT_INITIALIZE
;
; PURPOSE:
;	  Initialise variables from configuration file for EBDETECT
;
; CATEGORY:
;   Data analysis
;	
; CALLING SEQUENCE:
;	  result = EBDETECT_INITIALIZE(ConfigFile)
;
; INPUTS:
;	  ConfigFile  - input text file setting all parameters
;
; OUTPUTS:
;   Structure with parameters and switches used to run EBDETECT
;
; RESTRICTIONS:
;   Requires the following procedures and functions:
;   Procedures: N/A
;   Functions:  EBDETECT_PARSE_CONFIGFILE(), 
;
; EXAMPLE:
;   result = EBDETECT_INITIALIZE('ebdetect_config.txt')
;
; MODIFICATION HISTORY:
;   2016 Nov 30 Gregal Vissers: First version
;-
;
FUNCTION EBDETECT_INITIALIZE, ConfigFile, VERBOSE=verbose

  IF (N_PARAMS() LT 1) THEN BEGIN
    MESSAGE, 'Syntax: result = EBDETECT_INITIALIZE(ConfigFile [, /VERBOSE])', $
      /INFO
    RETURN, -1
  ENDIF 

  IF KEYWORD_SET(VERBOSE) THEN t_init = SYSTIME(/SECONDS)
  ; Create structure with default values
  result = { $
    ; File and directory names
    inputfile:'', lcsum_cube:'', comparison_mask:'', $
    inputdir:'./', outputdir:'./', $
    ; Detection parameters
    nw:1L, nx:1L, ny:1L, nt:1L, wsum_pos:0L, lcsum_pos:0L, $
    asecpix:[1.,1.], intensity_constraint:!VALUES.F_NAN, $
    sdev_mult_constraint:0., mean_mult_constraint:1., lc_sigma:0., lifetime_constraint:[0L,1L], $
    size_constraint:[1L,1L], overlap_constraint:1L, t_skip_constraint:0L, $
    limit_group_search:0L, $
    kernel_size:[1L,1L], running_mean:0L, $
    region_threshold:[0L,0L,0L,0L], remove_detections:-1, $
    override_merge:[-1,-1,-1,-1], $
    ; Switches
    sum_cube:0B, get_kernels:0B, get_centroids:1B, $
    lc_constraint:0B, merge_check:1B, split_check:1B, $
    read_detect_init:0B, write_detect_init:1B, $
    read_detect_overlap:0B, write_detect_overlap:1B,$
    write_detect_final:1B, write_mask:1B, write_inplace:1B, $
    exit_status:0B }
  tag_names_orig = TAG_NAMES(result)
  ntag_names_orig = N_ELEMENTS(tag_names_orig)
  dtypes = BYTARR(ntag_names_orig)
  result_orig = result
  result = CREATE_STRUCT(result, 'read_from_file', BYTARR(ntag_names_orig))
  FOR i=0,N_ELEMENTS(dtypes)-1 DO dtypes[i] = SIZE(result.(i), /TYPE)
 
  ; Checking existence of ConfigFile and if it does, process
  IF (N_ELEMENTS(ConfigFile) NE 1) THEN BEGIN
    ConfigFile = 'ebdetect_config.txt'
    EBDETECT_FEEDBACK, /WARNING, $
      'No configuration file has been supplied. '+$
      'Checking for default file in current working directory...'
  ENDIF
  result.exit_status = FILE_TEST(ConfigFile)
  IF (result.exit_status EQ 0) THEN BEGIN
    EBDETECT_FEEDBACK, /ERROR, /TERMINATE, $
      'No configuration file '+ConfigFile+$
      ' exists in the current working directory' 
    RETURN, result
  ENDIF ELSE BEGIN
    IF KEYWORD_SET(VERBOSE) THEN EBDETECT_FEEDBACK, 'Parsing '+ConfigFile+'...'
    ; Open ConfigFile for reading
    OPENR, lun, ConfigFile, /GET_LUN 
    nl = 0L 
    WHILE ~EOF(lun) DO BEGIN line = ''
      READF, lun, line    
      parsed_line = EBDETECT_PARSE_CONFIGFILE(line)
      IF (parsed_line.field NE '') THEN BEGIN
        IF KEYWORD_SET(VERBOSE) THEN EBDETECT_FEEDBACK, '   '+line
        ; Check which tag the current parsed_line.field corresponds to
        wheretag = WHERE(STRLOWCASE(tag_names_orig) EQ $
          STRLOWCASE(parsed_line.field), count)
        IF (count EQ 1) THEN BEGIN
          result.read_from_file[wheretag] = 1
          CASE dtypes[wheretag] OF
            1:  value = BYTE(FIX(parsed_line.value))
            2:  value = FIX(parsed_line.value)
            3:  value = LONG(parsed_line.value)
            4:  value = FLOAT(parsed_line.value)
            5:  value = DOUBLE(parsed_line.value)
            7:  value = parsed_line.value
            ELSE: value = FLOAT(parsed_line.value)    ; Failsafe
          ENDCASE
          ; Failsafe for REGION_THRESHOLD, which can be filename or 4-element
          ; array
          IF (STRLOWCASE(parsed_line.field) EQ 'region_threshold') THEN BEGIN
            IF (N_ELEMENTS(parsed_line.value) NE 1) THEN $
              value = LONG(parsed_line.value) $
            ELSE $
              value = parsed_line.value
          ENDIF

          ; If the size (i.e., number of elements) of the parsed value doesn't
          ; correspond with the default, delete the tag and reconstruct the
          ; result array with the read-in value
          IF (N_ELEMENTS(parsed_line.value) NE $
              N_ELEMENTS(result_orig.(wheretag))) THEN BEGIN
            result = EBDETECT_TAG_DELETE(result, parsed_line.field)
            result = CREATE_STRUCT(result, parsed_line.field, value)
          ENDIF ELSE BEGIN
            wheretag_new = WHERE(STRLOWCASE(TAG_NAMES(result)) EQ $
              STRLOWCASE(parsed_line.field))
            result.(wheretag_new) = value
          ENDELSE
        ENDIF
      ENDIF
    ENDWHILE
    
    ; Output default settings for variables not read from file
    IF KEYWORD_SET(VERBOSE) THEN BEGIN
      wherezero = WHERE(result.read_from_file EQ 0, count)
      IF (count GE 1) THEN BEGIN
        EBDETECT_FEEDBACK, ' '
        msg1 = '   The following keywords and settings were not found in '
        msg2 = '   the configuration file and were set to their defaults:'
        EBDETECT_FEEDBACK, STRJOIN(REPLICATE('=', STRLEN(msg1)))
        EBDETECT_FEEDBACK, msg1
        EBDETECT_FEEDBACK, msg2
        EBDETECT_FEEDBACK, STRJOIN(REPLICATE('=', STRLEN(msg1)))
        ; Skip over last item (as that is exit_status)
        FOR i=0,count-2 DO BEGIN
          value = result_orig.(wherezero[i])
          IF (SIZE(value, /TYPE) EQ 1) THEN value = FIX(value)
          IF (N_ELEMENTS(value) GT 1) THEN $
            value = '['+STRJOIN(STRTRIM(value,2),',')+']'
          EBDETECT_FEEDBACK, '   '+tag_names_orig[wherezero[i]]+$
            ' = '+STRTRIM(value,2)
        ENDFOR
      ENDIF
    ENDIF
  ENDELSE

  ; Clean up and return
  FREE_LUN, lun
  IF KEYWORD_SET(VERBOSE) THEN $
    EBDETECT_FEEDBACK, /DONE, T_INIT=t_init

  RETURN, result
                
END
