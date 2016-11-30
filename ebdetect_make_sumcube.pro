;+
; NAME:
;	  EBDETECT_MAKE_SUMCUBE
;
; PURPOSE:
;	  This procedure produces a summed data cube from an input data cube by summing over given 
;   wavelength indices.
;
; CATEGORY:
;   Data analysis
;	
; CALLING SEQUENCE:
;   EBDETECT_MAKE_SUMCUBE, Inputfile, Sum_positions
;	
;
; INPUTS:
;	  Inputfile - Inputfile in CRISPEX-ready format
;
; OPTIONAL INPUTS:
;   Sum_positions - Wavelength positions over which to sum. Defaults to all wavelengths (i.e.,
;                   INDGEN(NLP).
;
; KEYWORD PARAMETERS:
;   NLP           - Number of wavelength positions in inputfile. Defaults to 2.
;   NS            - Number of Stokes parameters in the input cube. Defaults to 1.
;   SET_NS        - Selected Stokes parameter for the summing of the cube. Defaults to 0.
;   FITS          - Flag identifying Inputfile as FITS file. Defaults to 0
;                   (i.e., file is legacy CRISPEX cube)
;   WRITE_INPLACE - Write frame by frame to disk. Defaults to 0 (i.e., cube is
;                   kept in memory before writing to disk)
;   OUTPUTFILENAME- Filename of output file. Defaults to
;                   'sum_'+FILE_BASENAME(Inputfile)
;   OUTDIR        - Output directory. Defaults to './'
;
; OUTPUTS:
;
; RESTRICTIONS:
;   Requires the following procedures and functions:
;   Procedures: LP_WRITE, EBDETECT_TIMER                [general]
;               LP_HEADER                               [in case of legacy cube]
;   Functions:  CRISPEX_FITSPOINTER(), FITS2IDL_TYPE()  [in case of FITS cube]
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   2016 Nov 30 Gregal Vissers: Taylored version from MK_SUMMED_CUBE
;-


PRO EBDETECT_MAKE_SUMCUBE, inputfile, sum_positions, NLP=nlp, NS=ns, SET_NS=set_ns, $
  NT=nt, FITS=fits, WRITE_INPLACE=write_inplace, OUTPUTFILENAME=outputfilename,$
  OUTDIR=outdir

  id_string = '$Id$'
  IF (N_PARAMS() LT 2) THEN BEGIN
    MESSAGE, id_string, /INFO
    MESSAGE,'Syntax: MK_SUMMED_CUBE, inputfile, sum_positions, NLP=nlp, '+$
      'NS=ns, SET_NS=set_ns, FITS=fits, WRITE_INPLACE=write_inplace, '+$
      'OUTPUTFILENAME=outputfilename, OUTDIR=outdir',/INFO
    RETURN
  ENDIF
	IF (N_ELEMENTS(NLP) NE 1) THEN nlp = 2
	IF (N_ELEMENTS(NS) NE 1) THEN ns = 1
	IF (N_ELEMENTS(SET_NS) NE 1) THEN set_ns = 0
	IF (N_ELEMENTS(SUM_POSITIONS) LT 1) THEN sum_positions = INDGEN(nlp)
  IF KEYWORD_SET(FITS) THEN BEGIN
    offset = CRISPEX_FITSPOINTER(inputfile, EXTEN_NO=0, header, /SILENT)
    nx = SXPAR(header, 'NAXIS1')
    ny = SXPAR(header, 'NAXIS2')
    nlp = SXPAR(header, 'NAXIS3')
    nt = SXPAR(header, 'NAXIS4')
    datatype = FITS2IDL_TYPE(header, /HEADER)
  ENDIF ELSE BEGIN
	  LP_HEADER, inputfile, NX=nx, NY=ny, NT=imnt, DATATYPE=datatype
    IF (N_ELEMENTS(NT) NE 1) THEN $;BEGIN
  	  nt = imnt/nlp/ns
    offset = 512
  ENDELSE

  IF (N_ELEMENTS(OUTDIR) NE 1) THEN $
    outdir = './'
  IF (N_ELEMENTS(OUTPUTFILENAME) NE 1) THEN $
    outputfilename = 'sum_'+FILE_BASENAME(inputfile)

  OPENR, lun, inputfile, /GET_LUN, SWAP_ENDIAN=KEYWORD_SET(FITS)
  readfile = ASSOC(lun,MAKE_ARRAY(nx,ny, TYPE=datatype),offset)
  IF ~KEYWORD_SET(WRITE_INPLACE) THEN $
    summed_cube = MAKE_ARRAY(nx,ny,nt, TYPE=datatype)
;  summed_cube = FLTARR(nx,ny,nt);MAKE_ARRAY(nx,ny,nt, TYPE=datatype)
	npass = nt*N_ELEMENTS(sum_positions)
	pass = 0
	t0 = SYSTIME(/SECONDS)
	FOR t=0,nt-1 DO BEGIN
    tmp_im = FLTARR(nx,ny) ;MAKE_ARRAY(nx,ny,TYPE=datatype)
		FOR lp=0,N_ELEMENTS(sum_positions)-1 DO BEGIN
      tmp_im += FLOAT(readfile[t*nlp*ns+sum_positions[lp]*ns+set_ns])
			pass += 1
			EBDETECT_TIMER, pass, npass, t0, EXTRA_OUTPUT='(t,lp,s)=('+STRTRIM(t,2)+','+$
                     STRTRIM(sum_positions[lp],2)+','+STRTRIM(set_ns,2)+').', $
                     /CALLBY
		ENDFOR
	  tmp_im /= FLOAT(N_ELEMENTS(sum_positions))
    IF KEYWORD_SET(WRITE_INPLACE) THEN $
      LP_PUT, tmp_im, outdir+outputfilename, t, nt=nt, $
        KEEP_OPEN=(t NE nt-1) $
    ELSE $
      summed_cube[*,*,t] = tmp_im
	ENDFOR
;	summed_cube /= FLOAT(N_ELEMENTS(sum_positions))
  IF ~KEYWORD_SET(WRITE_INPLACE) THEN $
  	LP_WRITE, summed_cube, outdir+outputfilename
	EBDETECT_FEEDBACK, /STATUS, 'Written: '+outputfilename
  
  FREE_LUN, lun

END
