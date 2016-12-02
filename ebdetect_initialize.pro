;+
; NAME:
;	  EBDETECT_INITIALIZE
;
; PURPOSE:
;	  Handles initialising the configuration file for EBDETECT
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
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   Structure with parameters and switches used to run EBDETECT
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;   Requires the following procedures and functions:
;   Procedures: N/A
;   Functions:  EBDETECT_PARSE_CONFIGFILE(), 
;
; PROCEDURE:
;
; EXAMPLE:
;   result = EBDETECT_INITIALIZE('ebdetect_config.txt')
;
; MODIFICATION HISTORY:
;   2016 Nov 30 Gregal Vissers: First version
;-
;
FUNCTION EBDETECT_INITIALIZE, ConfigFile, VERBOSE=verbose

  id_string = '$Id$'
  IF (N_PARAMS() LT 1) THEN BEGIN
    MESSAGE, id_string, /INFO
    MESSAGE, 'Syntax: result = EBDETECT_INITIALIZE(ConfigFile [, /VERBOSE])', $
      /INFO
    RETURN, -1
  ENDIF 

  ; Create structure with default values
  result = { $
    ; File and directory names
    inputfile:'', lcsum_cube:'', comparison_mask:'', $
    detect_init_file:'', inputdir:'./', outputdir:'./', $
    ; Detection parameters
    nlp:1L, nx:1L, ny:1L, nt:1L, wsum_pos:0L, lcsum_pos:0L, $
    sigma_constraint:0., lc_sigma:0., lifetime_constraint:[0L,1L], $
    size_constraint:[1L,1L], overlap_constraint:1L, t_skip_constraint:0L, $
    limit_group_search:0L, $
    loose_hysteresis:0B, kernel_size:[1L,1L], running_mean:0L, $
    region_threshold:[0L,0L,0L,0L], remove_detections:-1, $
    override_merge:[-1,-1,-1,-1], $
    ; Switches
    sum_cube:0B, factor_sigma:1B, double_set:1B, get_kernels:0B, $
    lc_constraint:0B, pad:1B, merge_check:1B, split_check:1B, $
    write_detect_init:1B, write_detect_overlap:1B, write_detect_final:1B, $
    write_mask:1B, write_inplace:1B, $
    exit_status:0B }
  dtypes = BYTARR(N_ELEMENTS(TAG_NAMES(result)))
  result_orig = result
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
        wheretag = WHERE(STRLOWCASE(TAG_NAMES(result_orig)) EQ $
          STRLOWCASE(parsed_line.field), count)
        IF (parsed_line.field EQ 'SIGMA_CONSTRAINT') THEN STOP
        IF (count EQ 1) THEN BEGIN
          CASE dtypes[wheretag] OF
            1:  value = BYTE(parsed_line.value)
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
              N_ELEMENTS(result.(wheretag))) THEN BEGIN
            result = EBDETECT_TAG_DELETE(result, parsed_line.field)
            result = CREATE_STRUCT(result, parsed_line.field, value)
          ENDIF ELSE $
            result.(wheretag) = value
        ENDIF
      ENDIF
    ENDWHILE
  ENDELSE

  ; Clean up and return
  FREE_LUN, lun
  IF KEYWORD_SET(VERBOSE) THEN $
    EBDETECT_FEEDBACK, /DONE

  RETURN, result
                
END
