;+
; NAME:
;	  EBDETECT_TAG_DELETE
;
; PURPOSE:
;	  Delete structure tags and corresponding values for EBDETECT
;
; CATEGORY:
;   Data analysis
;	
; CALLING SEQUENCE:
;	  result = EBDETECT_TAG_DELETE(Structure, Tags)
;
; INPUTS:
;	  Structure  - input structure
;   Tags       - single tag or array of tags to be removed from Structure
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;   VERBOSE    - output information on whether a certain tag is found or not
;
; OUTPUTS:
;   Structure where Tags have been removed.
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
;   result = EBDETECT_TAG_DELETE(Structure, 'TAG_NAME1')
;
; MODIFICATION HISTORY:
;   2016 Nov 30 Gregal Vissers: First version
;-
;
FUNCTION EBDETECT_TAG_DELETE, Structure, Tags, VERBOSE=verbose

  IF (N_PARAMS() LT 2) THEN BEGIN
    MESSAGE, 'Syntax: result = EBDETECT_TAG_DELETE(Structure, Tags, [, /VERBOSE])', $
      /INFO
    RETURN, -1
  ENDIF 
  
  ; Deletes Tags from Structure
  ntags = N_ELEMENTS(Tags) 
  IF (ntags GE 1) THEN BEGIN
    tagnames = TAG_NAMES(Structure)
    seltags = REPLICATE(1,N_ELEMENTS(tagnames))
    ; Check where the tags are located and set select array to 0 correspondingly
    FOR i=0,ntags-1 DO BEGIN
      wheretag = WHERE(STRLOWCASE(tagnames) EQ STRLOWCASE(tags[i]), count)
      IF (count NE 0) THEN BEGIN
        seltags[wheretag] = 0
        IF KEYWORD_SET(VERBOSE) THEN $
          MESSAGE, tags[i]+' found. Excluding from structure.',/INFO
      ENDIF ELSE IF KEYWORD_SET(VERBOSE) THEN $
          MESSAGE, tags[i]+' not found in structure. Skipping tag.', /INFO
    ENDFOR
    seltags = WHERE(seltags EQ 1)
    nseltags = N_ELEMENTS(seltags)
    ; Loop over remainder of tags and reconstruct the structure
    FOR i=0,nseltags-1 DO BEGIN
      IF (i NE 0) THEN $
        newstructure = CREATE_STRUCT(newstructure, $
          tagnames[seltags[i]], structure.(seltags[i])) $
      ELSE $
        newstructure = CREATE_STRUCT(tagnames[seltags[i]], $
                        structure.(seltags[i]))
    ENDFOR
  ENDIF
  RETURN, newstructure
END
