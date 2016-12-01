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
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   2016 Dec 01 Gregal Vissers: Taylored and simplified version from own
;                               ARRAY_COMPARE procedure (with init 2011 Oct 04)
;-
;
FUNCTION EBDETECT_ARRAY_COMPARE, Array1, Array2, VERBOSE=verbose

  id_string = '$Id'
	IF (N_PARAMS() LT 2) THEN BEGIN
    MESSAGE, '$Id', /INFO
		MESSAGE, 'Syntax: result = EBDETECT_ARRAY_COMPARE(Array1, Array2 '+$
      '[, /VERBOSE]', /INFO
		RETURN, -1
	ENDIF

	narray1 = N_ELEMENTS(array1)
	narray2 = N_ELEMENTS(array2)
	IF KEYWORD_SET(VERBOSE) THEN BEGIN
		PRINT,'Array1 has '+STRTRIM(narray1,2)+' elements.'
		PRINT,'Array2 has '+STRTRIM(narray2,2)+' elements.'
	ENDIF

	sorted_array1 = SORT(array1)
	sorted_array2 = SORT(array2)
	i = 0L
	j = 0L
	common_array = -1
	WHILE ((i LT narray1) AND (j LT narray2)) DO BEGIN
		IF (array1[i] EQ array2[j]) THEN BEGIN
			IF (common_array[0] NE -1) THEN BEGIN
        common_array = [common_array, array1[i]] 
        index_array1 = [index_array1, i]
        index_array2 = [index_array2, j]
      ENDIF ELSE BEGIN
        common_array = array1[i]
        index_array1 = i
        index_array2 = j
      ENDELSE
			i += 1L
			j += 1L
		ENDIF ELSE IF (array1[i] LT array2[j]) THEN $
      i += 1L $
    ELSE $
      j += 1L
	ENDWHILE
	
  ncommon_array = N_ELEMENTS(COMMON_ARRAY)

  IF KEYWORD_SET(VERBOSE) THEN PRINT, common_array
  result = {common_array:common_array, ncommon_array:ncommon_array, $
            index_array1:index_array1, index_array2:index_array2}
	
  RETURN, result

END
