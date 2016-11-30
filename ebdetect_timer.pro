;+
; NAME:
;	  EBDETECT_TIMER
;
; PURPOSE:
;   This routine outputs the expected time until completion for a loop
;
; CATEGORY:
;   Utilities, timer
;
; CALLING SEQUENCE:
;   PROCESS_TIMER, pass, npass, t0
;
; INPUTS:
;   pass  - counter of passes through the loop
;   npass - total number of passes of the loop
;   t0    - reference time in seconds, i.e., time before starting the loop
;
; KEYWORD PARAMETERS:
;   EXTRA_OUTPUT  - formatted string to be appended to the PROCESS_TIMER output
;   CALLBY        - scalar string specifying the procedure/function name calling PROCESS_TIMER.
;
; OUTPUTS:
;   String containing the current status of a process plus estimated time of completion, e.g.:
;     
;     Pass 3/300, or 1.000%. Running 00:00:02/00:03:20.
;
; RESTRICTIONS:
;   Requires the following procedures and functions:
;   Functions:  TIME2STRING()
;
; EXAMPLE:
;   EBDETECT_TIMER, t, nt, t0
;
; MODIFICATION HISTORY:
; 	2016 Nov 30 Gregal Vissers: Taylored version from PROCESS_TIMER
;-

PRO EBDETECT_TIMER, pass, npass, t0, EXTRA_OUTPUT=extra_output, $
  CALLBY=callby

  id_string = '$Id$'
	IF (N_PARAMS() LT 3) THEN BEGIN
    MESSAGE, id_string, /INFO
		MESSAGE, 'Syntax: EBDETECT_TIMER, pass, npass, t0 '+$
      '[, EXTRA_OUTPUT=extra_output] [, CALLBY=callby]', /INFO
		RETURN
	ENDIF
	
	IF (N_ELEMENTS(EXTRA_OUTPUT) NE 1) THEN extra_output = ''
  IF (N_ELEMENTS(CALLBY) NE 1) THEN $
    callby = '' $
  ELSE $
    callby = STRUPCASE(callby)+': '

	timer_t = SYSTIME(/SECONDS)                             ; Determine current time in seconds
	accumsectime = (timer_t-t0)                             ; Determine accumulated time in seconds
	totalsectime = (timer_t-t0)/FLOAT(pass)*FLOAT(npass)    ; Estimate ETA in seconds
	ndig = FLOOR(ALOG10(npass))+1                           ; Determine number of digits in completion
  leftsectime = totalsectime - accumsectime
;	accumsectime = STRMID(TIME2STRING(accumsectime),0,8)    ; Convert accumulated time to time string
	totalsectime = STRMID(TIME2STRING(totalsectime),0,8)    ; Convert estimated time to time string
	leftsectime = STRMID(TIME2STRING(leftsectime),0,8)    ; Convert estimated time to time string

  ; Write out results
	WRITEU, -1, STRING(FORMAT='(%"\r'+callby+'Progress ",f7.3,"% (",i'+STRTRIM(ndig,2)+',"/",i'+STRTRIM(ndig,2)+$
              ',"). Estimated time left: ",a8,". ",a'+STRTRIM(STRLEN(extra_output),2)+')',$
              ; Values
              STRTRIM(100.*pass/FLOAT(npass),2),STRTRIM(LONG(pass),2),STRTRIM(LONG(npass),2),$
              STRTRIM((leftsectime),2), extra_output)
  IF (pass EQ npass) THEN BEGIN
    PRINT,' Total time: '+STRTRIM((totalsectime),2)
  ENDIF

END
