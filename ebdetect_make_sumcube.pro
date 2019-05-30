;+
; NAME:
;	  EBDETECT_MAKE_SUMCUBE
;
; PURPOSE:
;	  Create summed data cube from an input data cube by summing over given ;
;   wavelength indices.
;
; CATEGORY:
;   Data analysis
;	
; CALLING SEQUENCE:
;   EBDETECT_MAKE_SUMCUBE, Inputfile, SumWavelengths
;	
;
; INPUTS:
;	  Inputfile       - Inputfile in CRISPEX-ready format
;
; OPTIONAL INPUTS:
;   SumWavelengths  - Wavelength positions over which to sum. Defaults to all
;                     wavelengths, i.e., INDGEN(NW).
;
; KEYWORD PARAMETERS:
;   NW              - Number of wavelength positions in inputfile. Defaults to 2.
;   NS              - Number of Stokes parameters in the input cube. Defaults to 1.
;   SET_NS          - Selected Stokes parameter for the summing of the cube. Defaults to 0.
;   FITS            - Flag identifying Inputfile as FITS file. Defaults to 0
;                     (i.e., file is legacy CRISPEX cube)
;   WRITE_INPLACE   - Write frame by frame to disk. Defaults to 0 (i.e., cube is
;                     kept in memory before writing to disk)
;   OUTPUTFILENAME  - Filename of output file. Defaults to
;                     'sum_'+FILE_BASENAME(Inputfile)
;   OUTDIR          - Output directory. Defaults to './'
;   VERBOSE         - User feedback verbosity. Defaults to not set.
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
;   2016 Nov 30 Gregal Vissers: Tailored version from MK_SUMMED_CUBE
;   2017 Oct 16 GV: Added handling of 2D array for SumWavelengths
;-


PRO EBDETECT_MAKE_SUMCUBE, Inputfile, SumWavelengths, NW=nw, NS=ns, SET_NS=set_ns, $
  NT=nt, FITS=fits, WRITE_INPLACE=write_inplace, OUTPUTFILENAME=outputfilename,$
  OUTDIR=outdir, VERBOSE=verbose

  IF (N_PARAMS() LT 2) THEN BEGIN
    MESSAGE,'Syntax: EBDETECT_MAKE_SUMCUBE, Inputfile [, SumWavelengths] '+$
      '[, NW=nw] [, NS=ns] [, SET_NS=set_ns] [, /FITS] [, /WRITE_INPLACE] '+$
      '[, OUTPUTFILENAME=outputfilename] [, OUTDIR=outdir] '+$
      '[, /VERBOSE]',/INFO
    RETURN
  ENDIF

	IF (N_ELEMENTS(NW) NE 1) THEN nw = 2
	IF (N_ELEMENTS(NS) NE 1) THEN ns = 1
	IF (N_ELEMENTS(SET_NS) NE 1) THEN set_ns = 0
	IF (N_ELEMENTS(SumWavelengths) LT 1) THEN SumWavelengths = INDGEN(nw)
  ; Failsafe against out of range SumWavelengths
  SumWavelengths = SumWavelengths > 0 < (nw-1)
  n_dims = SIZE(SumWavelengths, /N_DIMENSIONS)
  IF (n_dims LE 2) THEN BEGIN
    IF (n_dims EQ 2) THEN $
      nsums = (SIZE(SumWavelengths, /DIMENSIONS))[n_dims-1] $
    ELSE $
      nsums = 1

    ; Handle header information
    IF KEYWORD_SET(FITS) THEN BEGIN
      offset = CRISPEX_FITSPOINTER(Inputfile, EXTEN_NO=0, header, /SILENT)
      nx = SXPAR(header, 'NAXIS1')
      ny = SXPAR(header, 'NAXIS2')
      nw = SXPAR(header, 'NAXIS3')
      nt = SXPAR(header, 'NAXIS4')
      datatype = FITS2IDL_TYPE(header, /HEADER)
    ENDIF ELSE BEGIN
	    LP_HEADER, Inputfile, NX=nx, NY=ny, NT=imnt, DATATYPE=datatype
      IF (N_ELEMENTS(NT) NE 1) THEN nt = imnt/nw/ns
      offset = 512
    ENDELSE

    ; Set output location and filename
    IF (N_ELEMENTS(OUTDIR) NE 1) THEN $
      outdir = './'
    IF (N_ELEMENTS(OUTPUTFILENAME) NE 1) THEN $
      outputfilename = 'sum_'+FILE_BASENAME(Inputfile)

    ; Output file information 
    IF KEYWORD_SET(VERBOSE) THEN BEGIN
      EBDETECT_FEEDBACK, '  Input file dimensions: [nx,ny,nt,nw,ns]=['+$
        STRTRIM(nx,2)+','+STRTRIM(ny,2)+','+STRTRIM(nt,2)+','+$
        STRTRIM(nw,2)+','+STRTRIM(ns,2)+']'
      EBDETECT_FEEDBACK, '  Summing over positions: ['+$
        STRCOMPRESS(STRJOIN(SumWavelengths,','),/REMOVE_ALL)+']'
      EBDETECT_FEEDBACK, '  Writing result file "'+outputfilename+'" to "'+$
        outdir+'"'
      EBDETECT_FEEDBACK, /DONE 
    ENDIF

    ; Start reading input file
    OPENR, lun, Inputfile, /GET_LUN, SWAP_ENDIAN=KEYWORD_SET(FITS)
    readfile = ASSOC(lun,MAKE_ARRAY(nx,ny, TYPE=datatype),offset)
    IF ~KEYWORD_SET(WRITE_INPLACE) THEN $
      summed_cube = MAKE_ARRAY(nx,ny,nt*nsums, TYPE=datatype)
	  npass = nt*N_ELEMENTS(SumWavelengths)
	  pass = 0
	  t0 = SYSTIME(/SECONDS)
	  FOR t=0,nt-1 DO BEGIN
      tmp_im = FLTARR(nx,ny) 
      FOR ss=0,nsums-1 DO BEGIN
	  	  FOR w=0,N_ELEMENTS(SumWavelengths[*,ss])-1 DO BEGIN
          tmp_im += FLOAT(readfile[t*nw*ns+SumWavelengths[w,ss]*ns+set_ns])
	  	  	pass += 1
	  	  	EBDETECT_TIMER, pass, npass, t0, EXTRA_OUTPUT='(t,w,s)=('+STRTRIM(t,2)+','+$
                         STRTRIM(SumWavelengths[w],2)+','+STRTRIM(set_ns,2)+').', $
                         /CALLBY
	  	  ENDFOR
	      tmp_im /= FLOAT(N_ELEMENTS(SumWavelengths[*,ss]))
        IF KEYWORD_SET(WRITE_INPLACE) THEN $
          LP_PUT, tmp_im, outdir+outputfilename, t*nsums+ss, nt=nt*nsums, $
            KEEP_OPEN=(t NE nt-1) $
        ELSE $
          summed_cube[*,*,t*nsums+ss] = tmp_im
      ENDFOR
	  ENDFOR
    IF ~KEYWORD_SET(WRITE_INPLACE) THEN $
    	LP_WRITE, summed_cube, outdir+outputfilename
	  EBDETECT_FEEDBACK, /STATUS, 'Written: '+outputfilename
    
    FREE_LUN, lun
  ENDIF ELSE BEGIN
    EBDETECT_FEEDBACK, /ERROR, $
      'Can only parse 1D or 2D arrays as input to SumWavelengths variable!'
    RETURN
  ENDELSE

END
