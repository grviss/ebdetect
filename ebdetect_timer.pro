;+
; NAME:
;	  EBDETECT_TIMER
;
; PURPOSE:
;   Output expected time until completion for a loop
;
; CATEGORY:
;   Utilities, timer
;
; CALLING SEQUENCE:
;   EBDETECT_TIMER, Pass, nPass, t0
;
; INPUTS:
;   Pass  - counter of passes through the loop
;   nPass - total number of passes of the loop
;   t0    - reference time in seconds, i.e., time before starting the loop
;
; KEYWORD PARAMETERS:
;   EXTRA_OUTPUT  - formatted string to be appended to the EBDETECT_TIMER output
;   CALLBY        - scalar string specifying the procedure/function name calling
;                   EBDETECT_TIMER, or set as flag (in which case EBDETECT_TIMER 
;                   will find out by itself)
;   DONE          - add "done!" at end of extra output once done
;   TOTAL_TIME    - print out total time spend once done
;
; OUTPUTS:
;   String containing the current status of a process plus estimated time of
;   completion, e.g.:
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
;   2016 Dec 01 GV: Added DONE and TOTAL_TIME keywords
;-

PRO EBDETECT_TIMER, Pass, nPass, t0, EXTRA_OUTPUT=extra_output, $
  CALLBY=callby, DONE=done, TOTAL_TIME=total_time

	IF (N_PARAMS() LT 3) THEN BEGIN
		MESSAGE, 'Syntax: EBDETECT_TIMER, Pass, nPass, t0 '+$
      '[, EXTRA_OUTPUT=extra_output] [, CALLBY=callby]'+$
      '[, /DONE] [, /TOTAL_TIME]', /INFO
		RETURN
	ENDIF
	
	IF (N_ELEMENTS(EXTRA_OUTPUT) NE 1) THEN extra_output = ''
  IF (N_ELEMENTS(CALLBY) NE 1) THEN $
    callby = '> ' $
  ELSE BEGIN
    IF (SIZE(CALLBY, /TYPE) NE 7) THEN BEGIN
      traceback = SCOPE_TRACEBACK(/STRUCTURE)
      callby = '> '+traceback[N_ELEMENTS(traceback)-2].ROUTINE 
    ENDIF ELSE $
      callby = '> '+STRUPCASE(callby)
    callby += ': '
  ENDELSE

	timer_t = SYSTIME(/SECONDS)                             ; Determine current time in seconds
	accumsectime = (timer_t-t0)                             ; Determine accumulated time in seconds
	totalsectime = (timer_t-t0)/FLOAT(Pass)*FLOAT(nPass)    ; Estimate ETA in seconds
	ndig = FLOOR(ALOG10(nPass))+1                           ; Determine number of digits in completion
  leftsectime = totalsectime - accumsectime
	totalsectime = STRMID(TIME2STRING(totalsectime),0,8)    ; Convert estimated time to time string
	leftsectime = STRMID(TIME2STRING(leftsectime),0,8)    ; Convert estimated time to time string

  ; Write out results
	WRITEU, -1, STRING(FORMAT='(%"\r'+callby+'Progress ",f7.3,"% (",i'+STRTRIM(ndig,2)+',"/",i'+STRTRIM(ndig,2)+$
              ',"). Estimated time left: ",a8,". ",a'+STRTRIM(STRLEN(extra_output),2)+')',$
              ; Values
              STRTRIM(100.*Pass/FLOAT(nPass),2),STRTRIM(LONG(Pass),2),STRTRIM(LONG(nPass),2),$
              STRTRIM((leftsectime),2), extra_output)
  IF (Pass EQ nPass) THEN BEGIN
    IF KEYWORD_SET(DONE) THEN BEGIN
      IF (STRCOMPRESS(extra_output) EQ '') THEN $
        PRINT, ' Done!' $
      ELSE $
        PRINT, ' done!'
    ENDIF
    IF KEYWORD_SET(TOTAL_TIME) THEN $
      PRINT,' Total time: '+STRTRIM((totalsectime),2)
    IF ~KEYWORD_SET(DONE) AND ~KEYWORD_SET(TOTAL_TIME) THEN $
      PRINT,' ' ; Enter in case no final message is printed
  ENDIF

END
