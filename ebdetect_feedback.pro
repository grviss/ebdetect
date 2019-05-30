;+
; NAME:
;	  EBDETECT_FEEDBACK
;
; PURPOSE:
;	  Return feedback message during EBDETECT running.
;
; CATEGORY:
;   Data analysis
;	
; CALLING SEQUENCE:
;   EBDETECT_FEEDBACK, MessageString
;
; INPUTS:
;   None
;
; OPTIONAL INPUTS:
;	  MessageString - message to be displayed. Default: '' 
;
; KEYWORD PARAMETERS:
;   DONE      - Appends "Finished" to the displayed message.
;   STATUS    - Prepends "STATUS" to the displayed message.
;   ERROR     - Prepends "ERROR" to the displayed message.
;   WARNING   - Prepends "WARNING" to the displayed message.
;   TERMINATE - Appends "Terminating procedure" to the displayed message.
;   T_INIT    - Initial time used to compute elapsed time.
;
; OUTPUTS:
;   Message printed to the command line.
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
; PROCEDURE:
;   None
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   2016 Nov 30 Gregal Vissers: First version
;-
PRO EBDETECT_FEEDBACK, MessageString, DONE=done, STATUS=status, ERROR=error, $
  WARNING=warning, T_INIT=t_init, TERMINATE=terminate

  IF ((N_PARAMS() LT 1) AND ~KEYWORD_SET(DONE)) THEN BEGIN
    MESSAGE, 'Syntax: EBDETECT_FEEDBACK [, MessageString] [, /DONE] '+$
      '[, /STATUS] [, /ERROR] [, /WARNING] [, /TERMINATE] [, T_INIT=t_init]', $
      /INFO
    RETURN
  ENDIF 
  
  IF (N_ELEMENTS(MessageString) NE 1) THEN messagestring = ''
  traceback = SCOPE_TRACEBACK(/STRUCTURE)
  called_by = traceback[N_ELEMENTS(traceback)-2].ROUTINE + ': '
  type = ''
  IF KEYWORD_SET(STATUS) THEN type = 'STATUS : '
  IF KEYWORD_SET(ERROR) THEN type = 'ERROR : '
  IF KEYWORD_SET(WARNING) THEN type = 'WARNING : '
  messagestring_final = called_by + type + messagestring

  IF KEYWORD_SET(DONE) THEN BEGIN
    messagestring_final += ' Finished' 
    IF (N_ELEMENTS(T_INIT) EQ 1) THEN $
      messagestring_final += ' in '+ STRTRIM(SYSTIME(/SECONDS)-t_init,2) + $
      ' sec' $
    ELSE $
      messagestring_final += '!'
  ENDIF 

  IF KEYWORD_SET(TERMINATE) THEN $
    messagestring_final += ' Terminating procedure.'
    
  MESSAGE, messagestring_final, /INFO, /NONAME
  IF KEYWORD_SET(DONE) THEN PRINT, ' '

END
