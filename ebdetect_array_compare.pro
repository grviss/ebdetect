;+
; NAME:
;	  EBDETECT_ARRAY_COMPARE
;
; PURPOSE:
;	  Compares two arrays and returns common elements as well as their indices.
;
; CATEGORY:
;   Data analysis
;	
; CALLING SEQUENCE:
;   result = EBDETECT_ARRAY_COMPARE(Array1, Array2)
;
; INPUTS:
;   Array1  - input array of nArray1 elements
;   Array2  - input array of nArray2 elements
;
; OPTIONAL INPUTS:
;	  
;
; KEYWORD PARAMETERS:
;   VERBOSE - Prints out information on arrays. Defaults to 0.
;
; OUTPUTS:
;   Structure with 4 elements: 
;     - result.common_array:  the array of elements that are common to both
;                             Array1 and Array2
;     - result.ncommon_array: number of elements of result.common_array
;     - result.index_array1:  array with indices to subscript Array1 with to get
;                             result.common_array
;     - result.index_array2:  array with indices to subscript Array2 with to get
;                             result.common_array
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;   Requires the following function:
;     EBDETECT_ARRAY_APPEND()
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   2016 Dec 01 Gregal Vissers: Tailored and simplified version from own
;                               ARRAY_COMPARE procedure (with init 2011 Oct 04)
;-
;
FUNCTION EBDETECT_ARRAY_COMPARE, Array1, Array2, VERBOSE=verbose

	IF (N_PARAMS() LT 2) THEN BEGIN
    MESSAGE, '$Id$', /INFO
		MESSAGE, 'Syntax: result = EBDETECT_ARRAY_COMPARE(Array1, Array2 '+$
      '[, /VERBOSE])', /INFO
		RETURN, -1
	ENDIF

	narray1 = N_ELEMENTS(array1)
	narray2 = N_ELEMENTS(array2)
	IF KEYWORD_SET(VERBOSE) THEN BEGIN
		MESSAGE,'Array1 has '+STRTRIM(narray1,2)+' elements.', /INFO
		MESSAGE,'Array2 has '+STRTRIM(narray2,2)+' elements.', /INFO
	ENDIF

	i = 0L
	j = 0L
	common_array = !NULL
  index_array1 = !NULL
  index_array2 = !NULL
  ncommon_array = 0
	WHILE ((i LT narray1) AND (j LT narray2)) DO BEGIN
		IF (array1[i] EQ array2[j]) THEN BEGIN
      common_array = EBDETECT_ARRAY_APPEND(common_array, array1[i])
      index_array1 = EBDETECT_ARRAY_APPEND(index_array1, i)
      index_array2 = EBDETECT_ARRAY_APPEND(index_array2, j)
			i += 1L
			j += 1L
		ENDIF ELSE IF (array1[i] LT array2[j]) THEN $
      i += 1L $
    ELSE $
      j += 1L
	ENDWHILE

  IF (common_array NE !NULL) THEN BEGIN
    ncommon_array = N_ELEMENTS(COMMON_ARRAY) 
    result = {common_array:common_array, ncommon_array:ncommon_array, $
              index_array1:index_array1, index_array2:index_array2}
  ENDIF ELSE $
    result = {ncommon_array:ncommon_array}

  IF KEYWORD_SET(VERBOSE) THEN PRINT, common_array
	
  RETURN, result

END
