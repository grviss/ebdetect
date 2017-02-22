;+
; NAME:
;	  EBDETECT
;
; PURPOSE:
;	  This routine is meant for Ellerman bomb detection, but detects anything
;   above a certain intensity threshold and abiding by size and lifetime
;   constraints.
;
; CATEGORY:
;   Data analysis
;	
; CALLING SEQUENCE:
;   EBDETECT, ConfigFile	
;
; INPUTS:
;
; OPTIONAL INPUTS:
;	  ConfigFile  - input configuration file that contains information on all 
;                 data files, detection parameters and switches. If not provided, 
;                 EBDETECT will look for ebdetect_config.txt in the working
;                 directory 
;
; KEYWORD PARAMETERS:
;   OVERRIDE_PARAMS - Structure with tag names as defined in ConfigFile that
;                     should be overridden by other values. Should only be used
;                     when running a detection grid where up to a handful of
;                     parameters change for each run while the majority remains
;                     the same. 
;   VERBOSE         - Set verbosity level:
;                       0 = no feedback
;                       1 = initial parameters and progress timers
;                       2 = as 1, plus interim status reports, feedback movies
;                           and detection statistics
;                       3 = as 2, plus stopping in between major steps for
;                           debugging. Defaults to 0.
;   NO_PLOT         - Set if plots should be suppressed. Has no effect unless 
;                     VERBOSE is set to 2 or higher. Defaults to 0.
;
; OUTPUTS:
;   Depending on the switches set in the configuration file, the code will
;   output detection files and/or mask cubes with initial, intermediate and
;   final detections. 
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;   Requires the following procedures and functions:
;   Procedures: LP_HEADER, LP_WRITE, EBDETECT_MAKE_SUMCUBE, EBDETECT_TIMER,
;               EBDETECT_FEEDBACK
;   Functions:  EBDETECT_INITIALIZE(), EBDETECT_ARRAY_COMPARE(), LP_GET()
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   2012 March Gregal Vissers: First version
;   2012 Dec 12 GV: Updated version of old routine
;   2016 Sep 16 GV: Imported current version into CVS
;   2016 Nov 30 GV: Cleaned away all input keywords that are now populated from
;                   the configuration file and updated handling of input keywords
;   2016 Dec 02 GV: Version 1.0 
;   2017 Feb 21 GV: Added OVERRIDE_PARAMS keyword
;-
;
PRO EBDETECT, ConfigFile, OVERRIDE_PARAMS=override_params, VERBOSE=verbose, $
  NO_PLOT=no_plot

;============================================================================== 

  id_string = '$Id$'
  IF (N_PARAMS() LT 1) THEN BEGIN
    MESSAGE, id_string, /INFO
    MESSAGE, 'Syntax: EBDETECT, ConfigFile [, VERBOSE=verbose]', /INFO
    RETURN
  ENDIF

  t_init_overall = SYSTIME(/SECONDS)
  t_init = t_init_overall
  IF (N_ELEMENTS(VERBOSE) NE 1) THEN verbose = 0
 
  ; Get parameters and exit on error
  params = EBDETECT_INITIALIZE(ConfigFile, VERBOSE=(verbose GE 2))
  IF (params.exit_status EQ 0) THEN RETURN

  ; Check for overriding parameters from configuration file
  ; Can be useful when running grids of detections
  IF (N_ELEMENTS(OVERRIDE_PARAMS) GE 1) THEN BEGIN
    IF (SIZE(OVERRIDE_PARAMS, /TYPE) EQ 8) THEN BEGIN
      override_tags = TAG_NAMES(override_params)
      init_tags = TAG_NAMES(params)
      pass = 0
      FOR i=0,N_ELEMENTS(override_tags)-1 DO BEGIN
        wheretag = WHERE(STRLOWCASE(init_tags) EQ $
          STRLOWCASE(override_tags[i]), count)
        IF (count EQ 1) THEN BEGIN
          IF (pass EQ 0) THEN $
            EBDETECT_FEEDBACK, 'Overriding results from EBDETECT_INITIALIZE'
          IF (N_ELEMENTS(params.(wheretag)) EQ 1) THEN $
            oldvalue = STRTRIM(params.(wheretag),2) $
          ELSE $
            oldvalue = '['+STRCOMPRESS(STRJOIN(params.(wheretag),','),/REMOVE_ALL)+']'
          IF (N_ELEMENTS(override_params.(i)) EQ 1) THEN $
            newvalue = STRTRIM(override_params.(i),2) $
          ELSE $
            newvalue = '['+STRCOMPRESS(STRJOIN(override_params.(i),','),/REMOVE_ALL)+']'
          EBDETECT_FEEDBACK, STRJOIN(REPLICATE(' ',14))+$
            STRUPCASE(override_tags[i])+' = '+oldvalue+' --> '+newvalue
          IF (N_ELEMENTS(override_params.(i)) NE $
            N_ELEMENTS(params.(wheretag))) THEN BEGIN
            params = EBDETECT_TAG_DELETE(params, override_tags[i])
            params = CREATE_STRUCT(params, override_tags[i], override_params.(i))
          ENDIF ELSE BEGIN
            wheretag_new = WHERE(STRLOWCASE(TAG_NAMES(params)) EQ $
              STRLOWCASE(override_tags[i]))
            params.(wheretag_new) = override_params.(i)
          ENDELSE
          pass += 1
        ENDIF ELSE BEGIN
          EBDETECT_FEEDBACK, 'Tag '+override_tags[i]+' could not be found in '+$
            'EBDETECT parameter structure. Please check that the tag name '+$
            'supplied to OVERRIDE_PARAMS is correct.', /WARNING
        ENDELSE
      ENDFOR
    ENDIF ELSE BEGIN 
      EBDETECT_FEEDBACK,'OVERRIDE_PARAMS is not a structure. Parameters from '+$
        'configuration file will not be overridden.', /WARNING
    ENDELSE
  ENDIF

  ; Read in variables
  feedback_txt = 'processing variables and switches.'
  EBDETECT_FEEDBACK, feedback_txt+'..', /STATUS
  IF (STRMID(params.inputdir,0,1,/REVERSE) NE PATH_SEP()) THEN $
    params.inputdir += PATH_SEP()
  IF (STRMID(params.outputdir,0,1,/REVERSE) NE PATH_SEP()) THEN $
    params.outputdir += PATH_SEP()
  ; Check that directory exists, and if not create it
  IF (FILE_TEST(params.outputdir) EQ 0) THEN $
    FILE_MKDIR, params.outputdir
  ; Set summed cube file names
  IF KEYWORD_SET(params.sum_cube) THEN BEGIN
    wsum_cube_filename = params.inputfile 
    full_wsum_cube_filename = params.inputdir + wsum_cube_filename 
  ENDIF ELSE BEGIN
    wsum_cube_filename = 'wsum_'+FILE_BASENAME(params.inputfile)
    full_wsum_cube_filename = params.outputdir + wsum_cube_filename 
  ENDELSE
  lcsum_cube_filename = 'lcsum_'+FILE_BASENAME(params.inputfile)
  ; Check existence of input files
  inputfile_exists = 0
  wsum_cube_exists = 0
  lcsum_cube_exists= 0
  detect_init_file_exists = 0
  comparison_mask_exists = 0
  wsum_cube_exists = FILE_TEST(full_wsum_cube_filename)
  IF (params.inputfile NE '') THEN BEGIN
    IF (params.sum_cube EQ 0) THEN BEGIN
      inputfile_exists = FILE_TEST(params.inputdir+params.inputfile)
      IF inputfile_exists THEN BEGIN
        ; Get input file dimensions
      	LP_HEADER,params.inputdir+params.inputfile, NX=nx, NY=ny, NT=imnt 
    	  nt = imnt/params.nlp
      ENDIF ELSE BEGIN
        EBDETECT_FEEDBACK, /ERROR, /TERMINATE, $
          'No inputfile '+params.inputfile+' exists in directory '+$
          params.inputdir
        RETURN 
      ENDELSE
    ENDIF ELSE BEGIN
      IF wsum_cube_exists THEN $
        LP_HEADER, full_wsum_cube_filename, NX=nx, NY=ny, NT=nt  $
      ELSE BEGIN
        EBDETECT_FEEDBACK, /ERROR, /TERMINATE, $
          'No inputfile '+FILE_BASENAME(full_wsum_cube_filename)+' exists in directory '+$
          FILE_DIRNAME(full_wsum_cube_filename)
        RETURN 
      ENDELSE
    ENDELSE
  ENDIF ELSE BEGIN
    EBDETECT_FEEDBACK, /ERROR, /TERMINATE, $
      'Either inputfile or sum_cube should be supplied with an existing file.'+$
      ' Please check your input in '+ConfigFile+'.'
    RETURN
  ENDELSE
  params.nx = nx
  params.ny = ny
  params.nt = nt
  IF (params.lcsum_cube NE '') THEN $
    lcsum_cube_exists = FILE_TEST(params.inputdir + params.lcsum_cube) $
  ELSE $
    lcsum_cube_exists = FILE_TEST(params.outputdir + lcsum_cube_filename)
	IF (params.detect_init_file NE '') THEN detect_init_file_exists = $
    FILE_TEST(params.outputdir + params.detect_init_file)
  IF (params.comparison_mask NE '') THEN $
    comparison_mask_exists = FILE_TEST(params.inputdir + params.comparison_mask)

  ; Read-in of detection thresholds/constraints
  ; Intensity thresholds and switches
  intensity_constraint_set = 0
  intensity_constraint_idx = WHERE(STRLOWCASE(TAG_NAMES(params)) EQ $
    'intensity_constraint', count)
  IF (count NE 0) THEN $
    intensity_constraint_set = FINITE(params.intensity_constraint)
  IF (N_ELEMENTS(params.sigma_constraint) GE 1) THEN $
    sigma_constraint = params.sigma_constraint[SORT(params.sigma_constraint)] 
  IF intensity_constraint_set THEN $
    nlevels = N_ELEMENTS(params.intensity_constraint) $
  ELSE $
    nlevels = N_ELEMENTS(sigma_constraint)
  IF (N_ELEMENTS(params.lc_sigma) GT 1) THEN $
    lc_sigma = params.lc_sigma[SORT(params.lc_sigma)]
  IF (N_ELEMENTS(params.region_threshold) EQ 4) THEN $
    region_threshold_set = (TOTAL(params.region_threshold) NE 0) $
  ELSE $
    region_threshold_set = (STRCOMPRESS(params.region_threshold) NE '')
  ; Spatial thresholds
  min_size = params.size_constraint[0]
	IF (N_ELEMENTS(params.size_constraint) GT 1) THEN $
    max_size = params.size_constraint[1] $
  ELSE $
    max_size = LONG(params.nx)*LONG(params.ny)
  dataratio = params.nx/FLOAT(params.ny)
  ; Lifetime thresholds
  min_lifetime = params.lifetime_constraint[0]
  max_lifetime = params.nt
	IF (N_ELEMENTS(params.lifetime_constraint) EQ 2) THEN $
    max_lifetime = params.lifetime_constraint[1] $
  ELSE $
    max_lifetime = params.nt

  ; Set output filenames
  suffix = '.'+(STRSPLIT(FILE_BASENAME(params.inputfile), '.', /EXTRACT))[-1] 
  outfilename_base = 'ebdetect'
  sigma_label = '_stdev'
  IF (N_ELEMENTS(params.sigma_constraint) GT 1) THEN $
    sigma_label +=STRJOIN(STRTRIM(params.sigma_constraint,2),'-') $
  ELSE $
    sigma_label +=STRTRIM(params.sigma_constraint,2)
  IF params.lc_constraint THEN $
    lc_label = '_lcsdv'+STRTRIM(params.lc_sigma,2) $
  ELSE lc_label = ''
  time_label = '_lifet' 
  IF (N_ELEMENTS(params.lifetime_constraint) GT 1) THEN $
    time_label +=STRJOIN(STRTRIM(params.lifetime_constraint,2),'-') $
  ELSE $
    time_label +=STRTRIM(params.lifetime_constraint,2)

  
  outfilename_base += sigma_label+lc_label+time_label+$
    '_'+FILE_BASENAME(params.inputfile, suffix)
  running_mean_idlsave = outfilename_base + '_running_mean_sdev.idlsave'
  detect_init_idlsave = outfilename_base + '_init.idlsave'
  detect_overlap_idlsave = outfilename_base + '_overlap.idlsave'
  detect_final_idlsave = outfilename_base +'_final.idlsave'
  overlap_mask_filename = outfilename_base + '_overlap.maskcube'
  final_mask_filename = outfilename_base + '_final.maskcube'
  kernel_mask_filename = outfilename_base + '_final_kernels.maskcube'

  IF (verbose GE 2) THEN BEGIN
    IF params.write_detect_init THEN $
      EBDETECT_FEEDBACK, 'Initial detections will be written to:     '+$
        detect_init_idlsave
    IF params.write_detect_overlap THEN BEGIN
      EBDETECT_FEEDBACK, 'Overlapping detections will be written to: '+$
        detect_overlap_idlsave
      IF params.write_mask THEN $
        EBDETECT_FEEDBACK, 'Overlap detection mask will be written to: '+$
          overlap_mask_filename
    ENDIF
    IF params.write_detect_final THEN BEGIN
      EBDETECT_FEEDBACK, 'Final detections will be written to:       '+$
        detect_final_idlsave
      IF params.write_mask THEN $
        EBDETECT_FEEDBACK, 'Final detection mask will be written to:   '+$
          final_mask_filename
    ENDIF
  ENDIF

  EBDETECT_FEEDBACK, feedback_txt, /STATUS, /DONE, T_INIT=t_init
  IF (verbose EQ 3) THEN STOP

;===============================================================================
;======================== Create/set summed cube files =========================
;===============================================================================
  ; Create summed wing cube if multiple sum_positions are given
  IF (verbose GE 2) THEN BEGIN
    t_init = SYSTIME(/SECONDS)
    feedback_txt = 'Creating/setting summed cube files.'
    EBDETECT_FEEDBACK, feedback_txt+'..', /STATUS
  ENDIF
	IF ((wsum_cube_exists NE 1) AND (N_ELEMENTS(params.wsum_pos) GE 1)) THEN BEGIN
    IF (verbose GE 2) THEN EBDETECT_FEEDBACK, /STATUS, $
      '> Creating summed wing cube...'
    EBDETECT_MAKE_SUMCUBE, params.inputdir+params.inputfile, $
      params.wsum_pos, NLP=params.nlp, $
      OUTPUTFILENAME=params.outputdir+wsum_cube_filename, $
      WRITE_INPLACE=params.write_inplace, OUTDIR=params.outputdir
    sum_cube = params.outputdir+wsum_cube_filename
    IF (verbose GE 2) THEN EBDETECT_FEEDBACK, /STATUS, /DONE
  ENDIF ELSE $
    sum_cube = full_wsum_cube_filename
     
  ; Create summed "line-center" cube if multiple positions around line center
  ; are given
  IF (params.lc_constraint AND (lcsum_cube_exists NE 1) AND $
    (N_ELEMENTS(params.lcsum_pos) GE 1)) THEN BEGIN
    IF (verbose GE 2) THEN EBDETECT_FEEDBACK, /STATUS, $
      '> Creating summed line center cube...'
    EBDETECT_MAKE_SUMCUBE, params.inputdir+params.inputfile, $
      params.lcsum_pos, NLP=params.nlp, $
      OUTPUTFILENAME=params.outputdir+'lcsum_'+FILE_BASENAME(params.inputfile), $
      WRITE_INPLACE=params.write_inplace, OUTDIR=params.outputdir
    IF (verbose GE 2) THEN EBDETECT_FEEDBACK, /STATUS, /DONE
  ENDIF ELSE $
    lcsum_cube = params.inputdir+params.lcsum_cube

  IF (verbose GE 2) THEN BEGIN
    EBDETECT_FEEDBACK, feedback_txt, /STATUS, /DONE, T_INIT=t_init
  	IF (verbose EQ 3) THEN STOP
  ENDIF

;===============================================================================
;======================== Determine average intensities ========================
;===============================================================================
; Run first detection based on the intensity and size thresholds
; Supply SUM_CUBE with filename if not continuing from before
; Set keyword WRITE_FIRST_DETECT to write detections to file
	IF (detect_init_file_exists NE 1) THEN BEGIN
    IF (verbose GE 2) THEN BEGIN
      t_init = SYSTIME(/SECONDS)
      feedback_txt = 'Determining average intensities.'
      EBDETECT_FEEDBACK, feedback_txt+'..', /STATUS
    ENDIF
    IF (params.running_mean EQ 1) THEN BEGIN
      IF (SIZE(params.running_mean,/TYPE) NE 7) THEN BEGIN
        running_mean_summed_cube = FLTARR(params.nt)
        IF ~KEYWORD_SET(params.factor_sigma) THEN $
          running_sdev = FLTARR(params.nt) $
        ELSE $
          running_sdev = 'N/A'
        t0 = SYSTIME(/SECONDS)
        FOR t=0L,params.nt-1 DO BEGIN
          tlow = (t-params.running_mean) > 0
          tupp = (tlow + params.running_mean) < (params.nt-1)
          ; determine average 
          FOR tt=tlow,tupp DO BEGIN
            selpix = WHERE(LP_GET(params.region_threshold,tt) EQ 1)
            IF (t NE 0) THEN $
              tmp_mean_summed_cube = [running_mean_summed_cube, $
                                      (LP_GET(sum_cube,tt))[selpix]] $
            ELSE $
              tmp_mean_summed_cube = [(LP_GET(sum_cube,tt))[selpix]]
          ENDFOR
          running_mean_summed_cube[t] = MEAN(tmp_mean_summed_cube, /DOUBLE ,/NAN)
          ; Determine the standard deviation in the cube
          IF ~KEYWORD_SET(params.factor_sigma) THEN $
      		  running_sdev[t] = STDDEV(DOUBLE(tmp_mean_summed_cube),/NAN) 
          IF (verbose GE 1) THEN $
            EBDETECT_TIMER,t+1,params.nt,t0, /DONE, TOTAL_TIME=(verbose GE 2), $
              EXTRA_OUTPUT='Determining running mean...'
        ENDFOR
  			SAVE, running_mean_summed_cube, running_sdev, $
          FILENAME=params.outputdir+running_mean_idlsave
        EBDETECT_FEEDBACK, /SATUS, 'Written: '+$
          params.outputdir+running_mean_idlsave
      ENDIF ELSE RESTORE, params.running_mean, VERBOSE=(verbose GT 1)
    ENDIF ELSE BEGIN
      ; Read in summed wing cube
      summed_cube = FLTARR(params.nx,params.ny,params.nt)
      t0 = SYSTIME(/SECONDS)
      FOR t=0,params.nt-1 DO BEGIN
        summed_cube[*,*,t] = LP_GET(sum_cube,t)
        IF (verbose GE 1) THEN EBDETECT_TIMER,t+1,params.nt,t0, /DONE, $
          EXTRA_OUTPUT='Getting summed wing cube in memory...', $
          TOTAL_TIME=(verbose GE 2)
      ENDFOR
      ; Select only pixels as defined by REGION_THRESHOLD
      IF KEYWORD_SET(region_threshold_set) THEN BEGIN
        IF (N_ELEMENTS(params.region_threshold) EQ 4) THEN $
          sel_summed_cube = $
            summed_cube[params.region_threshold[0]:params.region_threshold[2],$
                        params.region_threshold[1]:params.region_threshold[3],*] $
        ELSE BEGIN
          t0 = SYSTIME(/SECONDS)
          FOR t=0L,params.nt-1 DO BEGIN
            selpix = WHERE(LP_GET(params.region_threshold,t) EQ 1, count)
            IF (count NE 0) THEN BEGIN
              IF (t EQ 0) THEN $
                sel_summed_cube = [(LP_GET(sum_cube,t))[selpix]] $
              ELSE $
                sel_summed_cube = [sel_summed_cube, (LP_GET(sum_cube,t))[selpix]]
            ENDIF
            IF (verbose GE 1) THEN $
              EBDETECT_TIMER,t+1,params.nt,t0, /DONE, TOTAL_TIME=(verbose GE 2), $
                EXTRA_OUTPUT='Determining selected summed wing cube pixels...'
          ENDFOR            
        ENDELSE
      ENDIF ELSE $
        sel_summed_cube = summed_cube
      IF ~KEYWORD_SET(params.factor_sigma) THEN $
        ; Determine the standard deviation in the cube
  		  sdev = STDDEV(DOUBLE(sel_summed_cube),/NAN) 
      ; Determine the average of the cube
  		mean_summed_cube = MEAN(sel_summed_cube, /DOUBLE,/NAN)             
    ENDELSE

    ; Determine line center constraints if any given
    IF KEYWORD_SET(params.lc_constraint) THEN BEGIN
      lc_summed_cube = FLTARR(params.nx,params.ny,params.nt)
      IF lcsum_cube_exists THEN $
        FOR t=0,params.nt-1 DO lc_summed_cube[0,0,t] = LP_GET(lcsum_cube,t) $
      ELSE BEGIN
        FOR t=0,params.nt-1 DO $
          lc_summed_cube[0,0,t] = LP_GET(inputfile,t*nlp+params.lcsum_pos)
      ENDELSE
      IF (N_ELEMENTS(params.region_threshold) GE 1) THEN BEGIN
        IF (N_ELEMENTS(params.region_threshold) EQ 4) THEN $
          sel_lc_summed_cube = $
            lc_summed_cube[params.region_threshold[0]:params.region_threshold[2],$
                           params.region_threshold[1]:params.region_threshold[3],*] $
        ELSE BEGIN
          t0 = SYSTIME(/SECONDS)
          FOR t=0L,params.nt-1 DO BEGIN
            lc_selpix = WHERE(LP_GET(params.region_threshold,t) EQ 1, count)
            IF (count NE 0) THEN BEGIN
              IF (t EQ 0) THEN $
                sel_lc_summed_cube = [(LP_GET(lcsum_cube,t))[lc_selpix]] $
              ELSE $
                sel_lc_summed_cube = [sel_lc_summed_cube, (LP_GET(lcsum_cube,t))[lc_selpix]]
            ENDIF
            IF (verbose GE 1) THEN $
              EBDETECT_TIMER,t+1,params.nt,t0, /DONE, TOTAL_TIME=(verbose GE 2), $
              EXTRA_OUTPUT='Determining selected summed line center cube pixels...'
          ENDFOR
        ENDELSE
      ENDIF ELSE $
        sel_lc_summed_cube = lc_summed_cube
      sdev_lc_cube = STDDEV(sel_lc_summed_cube, /DOUBLE, /NAN)
      mean_lc_cube = MEAN(sel_lc_summed_cube, /DOUBLE, /NAN)
    ENDIF
    IF (verbose GE 2) THEN $
      EBDETECT_FEEDBACK, feedback_txt, /STATUS, /DONE, T_INIT=t_init

;===============================================================================
;======================== Apply intensity thresholding =========================
;===============================================================================
    ; Define empty mask cube
		mask_cube = BYTARR(params.nx,params.ny,params.nt)                     
		totalpixels = LONG(params.nx)*LONG(params.ny)
		results = PTRARR(params.nt,/ALLOCATE_HEAP)
		pass = 0L
		totnstructs = 0L
		totnlabels = 0L
		IF (verbose GE 2) THEN BEGIN
      t_init = SYSTIME(/SECONDS)
      feedback_txt = 'Applying intensity thresholding.'
      EBDETECT_FEEDBACK, feedback_txt+'..', /STATUS
      IF ~KEYWORD_SET(NO_PLOT) THEN $
        WINDOW, XSIZE=750*dataratio, YSIZE=750, $
        TITLE='EBDETECT: Intensity thresholding'
    ENDIF
		t0 = SYSTIME(/SECONDS)
		FOR t=0L,params.nt-1 DO BEGIN
			mask = BYTARR(params.nx,params.ny)                           
			select_summed_cube = LP_GET(sum_cube,t)
      ; Select the pixels where cube intensity > mean intensity + sigma * stdev
      IF KEYWORD_SET(RUNNING_MEAN) THEN BEGIN
        mean_summed_cube = running_mean_summed_cube[t]
        sdev = running_sdev[t]
      ENDIF
      FOR s=0,nlevels-1 DO BEGIN                      
        ; Allow for hysteresis constraints
        IF intensity_constraint_set THEN $
          threshold = params.intensity_constraint[s] $
        ELSE BEGIN
          IF KEYWORD_SET(params.factor_sigma) THEN  $
            threshold = mean_summed_cube*params.sigma_constraint[s] $
          ELSE $
            threshold = mean_summed_cube+params.sigma_constraint[s]*sdev
        ENDELSE
  			wheregt = WHERE(select_summed_cube GT threshold, count)			
        ; Increase mask pixels gt constraint with 1
        IF (count NE 0) THEN $
    			mask[wheregt] += 1B                              
      ENDFOR

      ; Process line center condition, i.e., I < lc_threshold
      IF KEYWORD_SET(params.lc_constraint) THEN BEGIN
        IF lcsum_cube_exists THEN $
          select_lc_cube = LP_GET(lcsum_cube,t) $
        ELSE $
          select_lc_cube = LP_GET(inputfile,t*nlp+params.lcsum_pos)
        lc_threshold = mean_lc_cube + params.lc_sigma[0]*sdev_lc_cube
        wheregtlc = WHERE(select_lc_cube GT lc_threshold, lc_count)
        IF (lc_count NE 0) THEN mask[wheregtlc] = 0B
      ENDIF

      ; Pad mask
      pad_mask = BYTARR(params.nx+2,params.ny+2)
      pad_mask[1:params.nx,1:params.ny] = mask
      ; Where pixels gt lower threshold
      wheregt0 = WHERE(pad_mask GT 0, nwheregt0)          
      IF (nlevels GT 1) THEN $
        ; Where pixels gt upper threshold
        wheregt1 = WHERE(pad_mask GT 1, nwheregt1) $          
      ELSE BEGIN
        wheregt1  = wheregt0
        nwheregt1 = nwheregt0
      ENDELSE
      struct_mask = BYTARR(params.nx+2,params.ny+2) 
			IF (wheregt1[0] NE -1) THEN BEGIN
				nstructs = 0L
        ; Pixel coordinates to be discarded, init val
				discard_pix = -1                            
				FOR i=0L,nwheregt1-1 DO BEGIN									; Loop over all selected pixels
          ; If the considered pixel is not a discarded pixel, begin growing region
					IF (TOTAL(discard_pix EQ wheregt1[i]) LE 0) THEN BEGIN					
            IF (nwheregt1 GT 1) THEN BEGIN
            ; Grow the region of selected pixels touching the selected pixel
              IF KEYWORD_SET(params.loose_hysteresis) THEN $
    						structpix = REGION_GROW(pad_mask,wheregt1[i],/ALL)	$
              ELSE BEGIN
  						  structpix = REGION_GROW(pad_mask,wheregt1[i],/ALL, THRESH=[1,2])
              ENDELSE
            ENDIF ELSE structpix = wheregt1[i]
						nstructpix = N_ELEMENTS(structpix)
						IF (N_ELEMENTS(SIZE_CONSTRAINT) GE 1) THEN BEGIN
							IF ((nstructpix GE min_size) AND (nstructpix LE max_size)) THEN BEGIN
    						IF (nstructpix NE 1) THEN $     ; Added check for limb-to-limb
                  discard_pix = structpix[1:nstructpix-1] $
                ELSE $
                  discard_pix = -1
    						struct_mask[structpix] = 1B					; Add the region to the mask
    						nstructs += 1L
    						totnstructs += 1L
              ENDIF
						ENDIF ELSE BEGIN
    					IF (nstructpix NE 1) THEN $       ; Added check for limb-to-limb
                discard_pix = structpix[1:nstructpix-1] $
              ELSE $
                discard_pix = -1
							struct_mask[structpix] = 1B						; Add the region to the mask
							nstructs += 1L
							totnstructs += 1L
						ENDELSE
					ENDIF
				ENDFOR
        ; If there is a valid detection, initiate labelling and writing to
        ; results structure
        IF (TOTAL(WHERE(struct_mask GT 0)) NE -1) THEN BEGIN
          ; Label the pixels of all regions
					labels = LABEL_REGION(struct_mask,/ALL_NEIGHBORS)								
          labels = labels[1:params.nx,1:params.ny]                  
					nlabels = MAX(labels,/NAN)
					nlabels_pix = N_ELEMENTS(WHERE(labels GT 0))							
					nstruct_pix = N_ELEMENTS(WHERE(struct_mask GT 0))
					IF (nstruct_pix EQ nlabels_pix) THEN BEGIN							
            ; If nlabel_pix = nstructs_pix, then start writing results
						label_vals = LINDGEN(nlabels)+1
						structs = PTRARR(nlabels,/ALLOCATE_HEAP)
            ; Loop over all labels
						FOR j=0,nlabels-1 DO BEGIN								
              ; Select the pixels corresponding to current label
							positions = WHERE(labels EQ label_vals[j])					
              ; Write results to pointer
							*structs[j] = CREATE_STRUCT('label',label_vals[j]+totnlabels,'pos',positions)	
						ENDFOR
					ENDIF ELSE BEGIN										
            ; If number of structure pixels != the number of label pixels, stop with error
            EBDETECT_FEEDBACK, /ERROR, $
						  'Something is very wrong here... '+$
              STRTRIM(nstruct_pix,2)+' NE '+STRTRIM(nlabels_pix,2)
          ENDELSE
					IF ((verbose GE 2) AND ~KEYWORD_SET(NO_PLOT)) THEN $
            TVSCL,CONGRID(labels,512*dataratio,512)
				ENDIF ELSE BEGIN
					nlabels = 0
					structs = 0
				ENDELSE
				totnlabels += nlabels
			ENDIF ELSE BEGIN
				struct_mask = mask
				nlabels = 0
				structs = 0
			ENDELSE
			mask_cube[*,*,t] = mask
      ; Write time marker, number of detections and structures pointer to results pointer
			*results[t] = CREATE_STRUCT('t',t,'ndetect',nlabels,'structs',structs) 
			; Update timer and output extra information
      pass += 1L
			IF (verbose GE 1) THEN $
        EBDETECT_TIMER, pass, params.nt, t0, /DONE, TOTAL_TIME=(verbose GE 2), $
        EXTRA_OUTPUT=' Pixels detected: '+STRTRIM(nwheregt1,2)+'+'+$
        STRTRIM(nwheregt0,2)+'/'+STRTRIM(totalpixels,2)+$
        '. Detected structures: '+STRTRIM(nlabels,2)+'/'+STRTRIM(totnlabels,2)+'.'
		ENDFOR
    IF (verbose GE 2) THEN $
      EBDETECT_FEEDBACK, feedback_txt, /STATUS, /DONE, T_INIT=t_init
   ; Write thresholding detections to file
		IF KEYWORD_SET(params.write_detect_init) THEN BEGIN									
			ndetections = totnlabels
			SAVE, results, ndetections, FILENAME=params.outputdir+detect_init_idlsave
			EBDETECT_FEEDBACK,'Written: '+params.outputdir+detect_init_idlsave, /STATUS
		ENDIF
	  IF (verbose EQ 3) THEN STOP
	ENDIF


;================================================================================ 
;============================ Apply overlapping check ===========================
;================================================================================ 
; Overlap filter: only propagate cases for which certain overlap criteria are obeyed
; All detections at t=0 are "true"
; If a first detection file is supplied, restore it now
  IF (verbose GE 2) THEN BEGIN
    t_init = SYSTIME(/SECONDS)
    feedback_txt = 'Applying overlap check.'
    EBDETECT_FEEDBACK, feedback_txt+'..', /STATUS
  ENDIF
	IF detect_init_file_exists THEN BEGIN
    RESTORE, params.outputdir+params.detect_init_file							
    EBDETECT_FEEDBACK, '> Restored detections from file: '+$
      params.outputdir+params.detect_init_file
  ENDIF
	pass = 0L
;	totpasses = 0L
	tt = 0
	first_detect = 0
  ; Necessary while loop as the first frame(s) might not contain any detection
	WHILE (first_detect EQ 0) DO BEGIN
		IF ((*results[tt]).ndetect GT 0) THEN BEGIN
			detect_counter = LONG((*(*results[tt]).structs[(*results[tt]).ndetect-1]).label)
			first_detect = 1
		ENDIF
		tt += 1
	ENDWHILE
;	FOR t=0L,params.nt-2 DO totpasses += LONG((*results[t]).ndetect)
	t0 = SYSTIME(/SECONDS)
  ; Loop over all but the last time step
	FOR t=0L,params.nt-1 DO BEGIN													
    ; Loop over all detections at the current time step
		FOR j=0,(*results[t]).ndetect-1 DO BEGIN									
			pass += 1L
      ; If the label of the current detection is bigger than the detection counter
			IF ((*(*results[t]).structs[j]).label GT detect_counter) THEN BEGIN					
        ; Increase the detection counter by 1
				detect_counter += 1L										
        ; Relabel the current detection with the updated detection counter
				(*(*results[t]).structs[j]).label = detect_counter						
			ENDIF
			orig_detection = (*(*results[t]).structs[j]).pos

			;;; Check for splitting events ;;;
			overlapped = 0
			ncor = 0
			k_array = -1
			ncomarr = -1
			t_usel = t
			t_ubound = (t + params.t_skip_constraint) < (params.nt-1)
			WHILE ((overlapped EQ 0) AND (t_usel LT t_ubound)) DO BEGIN
				IF (t_usel LT t_ubound) THEN t_usel += 1
        ; Loop over all detections at the next time step
				FOR k=0,(*results[t_usel]).ndetect-1 DO BEGIN								
					comp_detection = (*(*results[t_usel]).structs[k]).pos
          ; Find the common elements between the considered detections
          array_compare = EBDETECT_ARRAY_COMPARE(orig_detection, comp_detection)
;					position_label = ' (t,det_orig,t_comp,det_comp)=('+STRTRIM(t,2)+','+$
;            STRTRIM(j,2)+','+STRTRIM(t_usel,2)+','+STRTRIM(k,2)+$
;            '). Single detections: '+STRTRIM(detect_counter,2)+'.'
          ; If the number of common elements >= overlap constraint
					IF ((N_ELEMENTS(array_compare.common_array) GE params.overlap_constraint) AND $
              (TOTAL(array_compare.common_array) NE -1)) THEN BEGIN		
						IF (TOTAL(k_array) NE -1) THEN $
              k_array = [k_array,k] $
            ELSE $
              k_array = k 
						IF (TOTAL(ncomarr) NE -1) THEN $
              ncomarr = [ncomarr,array_compare.ncommon_array] $
            ELSE $
              ncomarr = array_compare.ncommon_array 
						overlapped = 1
						ncor += 1
					ENDIF
				ENDFOR
			ENDWHILE

      ; If overlap occurs, assign the labels
			IF overlapped THEN BEGIN											
        ; If there is only one that overlaps
        IF (ncor EQ 1) THEN BEGIN
					oldlabel = (*(*results[t_usel]).structs[k_array[0]]).label 
          ; Assign the next detection the current detection label
					(*(*results[t_usel]).structs[k_array[0]]).label = $
            (*(*results[t]).structs[j]).label 		
          extraout = 'Overlap: '+STRTRIM(oldlabel,2)+' > '+$
            STRTRIM((*(*results[t_usel]).structs[k_array[0]]).label,2)+','+$
            STRTRIM(ncomarr[0],2)
        ENDIF ELSE BEGIN
          ; If there are multiple that overlap, determine which one has the
          ; biggest overlap
					wheremaxoverlap = WHERE(ncomarr EQ MAX(ncomarr,/NAN),$
            COMPLEMENT=wherenotmaxoverlap,NCOMPLEMENT=nwherenotmaxoverlap)
					oldlabel = (*(*results[t_usel]).structs[k_array[wheremaxoverlap[0]]]).label 
          ; Assign the next detection the current detection label
					(*(*results[t_usel]).structs[k_array[wheremaxoverlap[0]]]).label = $
            (*(*results[t]).structs[j]).label 		
          extraout = 'Overlap: '+STRTRIM(oldlabel,2)+' > '+$
            STRTRIM((*(*results[t_usel]).structs[k_array[wheremaxoverlap[0]]]).label,2)+','+$
            STRTRIM(ncomarr[0],2)
					FOR kk=0,nwherenotmaxoverlap-1 DO BEGIN
						oldlabel = (*(*results[t_usel]).structs[k_array[wherenotmaxoverlap[kk]]]).label 
						detect_counter += 1L										; - Increase the detection counter by 1
            ; Assign the next detection the current detection label
						(*(*results[t_usel]).structs[$
              k_array[wherenotmaxoverlap[kk]]]).label = detect_counter 		
					ENDFOR
				ENDELSE
			ENDIF ELSE $
        extraout = 'No overlap: '+STRTRIM((*(*results[t]).structs[j]).label,2)
      ; Output timer and extra information
      IF (verbose GE 1) THEN EBDETECT_TIMER, t+1, params.nt, t0, $
        EXTRA=extraout, TOTAL_TIME=(verbose GE 2)
		ENDFOR
	ENDFOR
	last_detect_counter = detect_counter
  IF (verbose GE 2) THEN BEGIN
    PRINT, ''
    EBDETECT_FEEDBACK, feedback_txt, /STATUS, /DONE, T_INIT=t_init
  	IF (verbose EQ 3) THEN STOP
  ENDIF

;================================================================================
;=========================== Check for merging events ===========================
;================================================================================
	IF KEYWORD_SET(params.merge_check) THEN BEGIN
    IF (verbose GE 1) THEN BEGIN
      IF (verbose GE 2) THEN BEGIN
        t_init = SYSTIME(/SECONDS)
        feedback_txt = 'Performing merge check.'
        EBDETECT_FEEDBACK, feedback_txt+'..', /STATUS
      ENDIF ELSE $
        PRINT,''
    ENDIF
		pass = 0L
		totpasses = 0L
		FOR t=0,params.nt-1 DO totpasses += LONG((*results[t]).ndetect)
		t0 = SYSTIME(/SECONDS)
    ; Loop over all but the last time step
		FOR t_dum=0,params.nt-1 DO BEGIN													
			t = params.nt-t_dum-1
      ; Loop over all detections at the current time step
			FOR j=0,(*results[t]).ndetect-1 DO BEGIN									
				pass += 1L
				orig_detection = (*(*results[t]).structs[j]).pos
				overlapped = 0
				ncor = 0
				k_array = -1
				ncomarr = -1
				t_lsel = t
				t_lbound = (t - params.t_skip_constraint) > 0
        ; Check labels and overlap
				WHILE ((overlapped EQ 0) AND (t_lsel GT t_lbound)) DO BEGIN						
					IF (t_lsel GT t_lbound) THEN t_lsel -= 1
					FOR k=0,(*results[t_lsel]).ndetect-1 DO BEGIN								
            ; Loop over all detections at the next time step
						comp_detection = (*(*results[t_lsel]).structs[k]).pos
            ; Find the common elements between the considered detections
            array_compare = EBDETECT_ARRAY_COMPARE(orig_detection, $
              comp_detection)
            ; If the number of common elements >= overlap constraint
						IF ((N_ELEMENTS(array_compare.common_array) GE $
                params.overlap_constraint) AND $
                (TOTAL(array_compare.common_array) NE -1)) THEN BEGIN		
						  IF (TOTAL(k_array) NE -1) THEN $
                k_array = [k_array,k] $
              ELSE $
                k_array = k 
						  IF (TOTAL(ncomarr) NE -1) THEN $
                ncomarr = [ncomarr,array_compare.ncommon_array] $
              ELSE $
                ncomarr = array_compare.ncommon_array 
							overlapped = 1
							ncor += 1
						ENDIF
					ENDFOR
				ENDWHILE			
        
        ; If overlap occured, assign labels
				IF overlapped THEN BEGIN										
					IF (ncor GT 1) THEN BEGIN
						wheremaxoverlap = WHERE(ncomarr EQ MAX(ncomarr,/NAN))
						oldlabel = (*(*results[t]).structs[j]).label
						newlabel = (*(*results[t_lsel]).structs[k_array[wheremaxoverlap[0]]]).label 
            extraout = 'Overlap: '+STRTRIM(oldlabel,2)+' > '+$
              STRTRIM(newlabel,2)+','+STRTRIM(ncomarr[wheremaxoverlap[0]],2)
            ;Assign the current detection the previous detection label
						(*(*results[t]).structs[j]).label = newlabel						
						FOR tt = t+1,nt-1 DO BEGIN
							kk = 0
							newlabel_set = 0
							WHILE ((newlabel_set EQ 0) AND (kk LT (*results[tt]).ndetect-1)) DO BEGIN
								kk += 1	
								IF ((*(*results[tt]).structs[kk]).label EQ oldlabel) THEN BEGIN
									(*(*results[tt]).structs[kk]).label = newlabel
									newlabel_set = 1
								ENDIF
							ENDWHILE
						ENDFOR
					ENDIF
				ENDIF
        IF (verbose GE 1) THEN EBDETECT_TIMER,t_dum+1,nt,t0,EXTRA=extraout, $
          TOTAL_TIME=(verbose GE 2)
			ENDFOR
		ENDFOR
	ENDIF

;================================================================================
;=========================== Override merging events ============================
;================================================================================
	IF (N_ELEMENTS(params.override_merge) EQ 4) THEN BEGIN
    IF (TOTAL(params.override_merge) NE -4) THEN BEGIN
  		IF (verbose GE 2) THEN $
        EBDETECT_FEEDBACK, /STATUS, 'Overriding merge detection:'
  		FOR t=params.override_merge[0],params.override_merge[1] DO BEGIN
  			FOR k=0,(*results[t]).ndetect-1 DO BEGIN
  				oldlabel = (*(*results[t]).structs[k]).label 
  				IF (oldlabel EQ params.override_merge[2]) THEN BEGIN
  					(*(*results[t]).structs[k]).label = params.override_merge[3]
            IF (verbose GE 2) THEN $
              EBDETECT_FEEDBACK, '(t,detection,oldlabel,newlabel) => ('+$
  					    STRTRIM(t,2)+','+STRTRIM(k,2)+','+STRTRIM(oldlabel,2)+$
                STRTRIM(params.override_merge[3],2)+')'
  				ENDIF ELSE $
            EBDETECT_FEEDBACK, 'No detection labelled '+STRTRIM(oldlabel,2)+$
            ' found at '+STRTRIM(t,2)
  			ENDFOR
  		ENDFOR
    ENDIF
	ENDIF

  IF (verbose GE 2) THEN BEGIN
    EBDETECT_FEEDBACK, feedback_txt, /STATUS, /DONE, T_INIT=t_init
    EBDETECT_FEEDBACK, /STATUS, $
      'Final number of single detections: '+STRTRIM(detect_counter,2)
	  IF (verbose EQ 3) THEN STOP
  ENDIF

;================================================================================
;=========================== Display and write results ==========================
;================================================================================
	IF ((verbose GE 2) OR KEYWORD_SET(params.write_detect_overlap)) THEN BEGIN
    ; Prep before looping over time
		IF (verbose GE 2) THEN BEGIN
      t_init = SYSTIME(/SECONDS)
      feedback_txt = 'Writing interim detection results.'
      EBDETECT_FEEDBACK, feedback_txt+'..', /STATUS
			IF ~KEYWORD_SET(NO_PLOT) THEN $
        WINDOW,XSIZE=750*dataratio,YSIZE=750, $
        TITLE='EBDETECT: Continuity constraints'
		ENDIF
    ; Define cube to be saved if not writing in place
		IF (KEYWORD_SET(params.write_detect_overlap) AND $
        KEYWORD_SET(params.write_mask)) THEN BEGIN
      IF ~KEYWORD_SET(params.write_inplace) THEN $
        overlap_mask_cube = BYTARR(params.nx,params.ny,params.nt)
    ENDIF
    ; Loop over time and display results and/or write results to disk
		FOR t=0,params.nt-1 DO BEGIN
			mask = BYTARR(params.nx,params.ny)
			FOR j=0,(*results[t]).ndetect-1 DO $
        mask[(*(*results[t]).structs[j]).pos] = 1B
			IF ((verbose GE 2) AND ~KEYWORD_SET(NO_PLOT)) THEN BEGIN
				TV,CONGRID(BYTSCL(LP_GET(sum_cube,t),/NAN),750*dataratio,750)
				LOADCT,13,/SILENT
				CONTOUR,CONGRID(mask,750*dataratio,750),COLOR=255, LEVELS = 1, /ISOTROPIC, $
          XS=13,YS=13,POSITION=[0,0,1,1],/NORMAL,/NOERASE
        IF comparison_mask_exists THEN $
          CONTOUR,REFORM(LP_GET(params.inputdir+params.comparison_mask,t)), $
            COLOR=200, LEVELS=1, /ISO, XS=13, YS=13, $
                  POS=[0,0,1,1], /NORMAL, /NOERASE
				LOADCT,0,/SILENT
				FOR j=0,(*results[t]).ndetect-1 DO BEGIN
					xyout_pos = ARRAY_INDICES(mask,((*(*results[t]).structs[j]).pos)[0])+[-5,-20]
					XYOUTS,xyout_pos[0]/FLOAT(nx),xyout_pos[1]/FLOAT(ny),$
            STRTRIM((*(*results[t]).structs[j]).label,2), $
            COLOR=255, /NORMAL, CHARSIZE=2
				ENDFOR
				XYOUTS,10.,10.,'All detections at t='+STRTRIM(t,2),/DATA,COLOR=255,CHARSIZE=2
				WAIT,0.05
			ENDIF
			IF (KEYWORD_SET(params.write_detect_overlap) AND $
          KEYWORD_SET(params.write_mask)) THEN BEGIN
        IF KEYWORD_SET(params.write_inplace) THEN $
          LP_PUT, mask, params.outputdir+overlap_mask_filename, t, nt=params.nt, $
            KEEP_OPEN=(t NE params.nt-1) $
        ELSE $
          overlap_mask_cube[*,*,t] = mask
      ENDIF
		ENDFOR
    ; Save results
		IF KEYWORD_SET(params.write_detect_overlap) THEN BEGIN
      ; Write detection save file
			ndetections = detect_counter
			SAVE, results, ndetections, FILENAME=params.outputdir+detect_overlap_idlsave
			EBDETECT_FEEDBACK,'> Written: '+params.outputdir+detect_overlap_idlsave, /STATUS
      ; Write cube if needed
      IF (KEYWORD_SET(params.write_mask) AND $
          ~KEYWORD_SET(params.write_inplace)) THEN BEGIN
  			LP_WRITE, overlap_mask_cube, params.outputdir+overlap_mask_filename
        EBDETECT_FEEDBACK, /STATUS, $
          '> Written: '+params.outputdir+overlap_mask_filename
      ENDIF
		ENDIF
	ENDIF
  IF (verbose GE 2) THEN BEGIN
    EBDETECT_FEEDBACK, feedback_txt, /STATUS, /DONE, T_INIT=t_init
  	IF (verbose EQ 3) THEN STOP
  ENDIF

;================================================================================
;============== Group detections and apply lifetime constraints =================
;================================================================================
	; Group detections by label AKA determine lifetimes
  IF (verbose GE 2) THEN BEGIN
    t_init = SYSTIME(/SECONDS)
    feedback_txt = 'Grouping detections and applying lifetime constraints.'
    EBDETECT_FEEDBACK, feedback_txt+'..', /STATUS
  ENDIF
	detections = PTRARR(detect_counter,/ALLOCATE_HEAP)
	sel_detect_idx = -1
	t0 = SYSTIME(/SECONDS)
  lifetime_max = 0L
	FOR d=0L,detect_counter-1 DO BEGIN											; Loop over all single detections
		t_arr = -1
		j_arr = -1
    label_check = d+1L
    IF (params.limit_group_search NE 0) THEN BEGIN
      ; Find first occurrence of current detection counter
      t_first = -1L
      t=0L
      WHILE (t_first EQ -1) DO BEGIN
        ; Loop over all detections at each time step
			  FOR j=0,(*results[t]).ndetect-1 DO BEGIN								
          ; If the detection label equals the current detection counter
				  IF ((*(*results[t]).structs[j]).label EQ label_check) THEN t_first = t 
        ENDFOR
        t += 1L
      ENDWHILE
      t_first_last = t_first
      t_low = t_first - LONG(params.limit_group_search/2.)
      IF (t_low LT 0) THEN BEGIN
        t_low = 0L
        t_upp = LONG(params.limit_group_search)
      ENDIF ELSE BEGIN
        t_upp = t_first + LONG(params.limit_group_search/2.)
        IF (t_upp GT (nt-1)) THEN BEGIN
          t_upp = LONG(nt)-1L
          t_low = t_upp - LONG(params.limit_group_search)
        ENDIF
      ENDELSE
    ENDIF ELSE BEGIN
      t_low = 0L
      t_upp = LONG(nt)-1L
    ENDELSE
    ; Loop over all time steps
		FOR t=t_low,t_upp DO BEGIN												
      ; Loop over all detections at each time step
			FOR j=0,(*results[t]).ndetect-1 DO BEGIN								
        ; If the detection label equals the current detection counter
				IF ((*(*results[t]).structs[j]).label EQ label_check) THEN BEGIN					
          ; Add the time step to the time step array
					IF (TOTAL(t_arr) NE -1) THEN $
            t_arr = [t_arr,t] $
          ELSE $
            t_arr = t
          ; - Add the detection number for that time step to an array
					IF (TOTAL(j_arr) NE -1) THEN $
            j_arr = [j_arr,j] $	
          ELSE $
            j_arr = j 
				ENDIF
			ENDFOR
		ENDFOR
		nt_arr = N_ELEMENTS(t_arr)
		lifetime = t_arr[nt_arr-1] - t_arr[0] + 1									; Determine lifetime
    IF (lifetime GT lifetime_max) THEN lifetime_max = lifetime
		; Checking lifetime constraint
    ; If the detections lifetime >= lifetime constraint
		IF ((lifetime GE min_lifetime) AND (lifetime LE max_lifetime)) THEN BEGIN		
      ; Add the detection label to the array of selected detections
			IF (TOTAL(sel_detect_idx) NE -1) THEN $
        sel_detect_idx = [sel_detect_idx,d]	$
      ELSE $
        sel_detect_idx = d 
      extra = 'Selected:    '
		ENDIF ELSE extra = 'Not selected:'
		det = PTRARR(nt_arr,/ALLOCATE_HEAP)
    ; Loop over all time steps where detection is present
		FOR tt=0L,nt_arr-1 DO $
			*det[tt] = CREATE_STRUCT('pos',(*(*results[t_arr[tt]]).structs[j_arr[tt]]).pos)
    ; Write results grouped by detection with lifetime information
		*detections[d] = CREATE_STRUCT('label',label_check,'t',t_arr,$
      'lifetime',lifetime,'det',det)				
		IF (verbose GE 1) THEN $
      EBDETECT_TIMER, label_check, detect_counter, t0, $
        EXTRA=extra+' d='+STRTRIM(d,2)+', nt='+$
        STRTRIM(nt_arr,2)+', t_upp='+STRTRIM(t_arr[nt_arr-1],2)+$
        ', t_low='+STRTRIM(t_arr[0],2)+$
        ', t='+STRTRIM(lifetime,2)+'. So far t_max='+STRTRIM(lifetime_max,2), $
        TOTAL_TIME=(verbose GE 2)
	ENDFOR
  IF (verbose GE 2) THEN BEGIN
    PRINT,''
    EBDETECT_FEEDBACK, feedback_txt, /STATUS, /DONE, T_INIT=t_init
    EBDETECT_FEEDBACK, /STATUS, $
	    'Final number of detections after lifetime constraint: '+$
      STRTRIM(N_ELEMENTS(sel_detect_idx),2)
    EBDETECT_FEEDBACK, /STATUS, $
      'Detection indices: ['+STRJOIN(STRTRIM(sel_detect_idx,2),',')+']'
	  IF (verbose EQ 3) THEN STOP
  ENDIF
	
;================================================================================
;===================== Construct final detection output =========================
;================================================================================
; Includes removing detections if so specified as well as getting kernels
	nsel_detections_orig = N_ELEMENTS(sel_detect_idx)
  IF (params.remove_detections[0] NE -1) THEN $
    nremove_detections = N_ELEMENTS(params.remove_detections) $
  ELSE $
    nremove_detections = 0
  nsel_detections = nsel_detections_orig - nremove_detections
	sel_detections = PTRARR(nsel_detections,/ALLOCATE_HEAP)
	sel_detect_mask = BYTARR(params.nx,params.ny,params.nt)
	IF (verbose GE 2) THEN BEGIN
    t_init = SYSTIME(/SECONDS)
    feedback_txt1 = 'Constructing final detection output.'
    EBDETECT_FEEDBACK, feedback_txt1+'..', /STATUS
		IF ~KEYWORD_SET(NO_PLOT) THEN BEGIN
      WINDOW,XSIZE=750*dataratio,YSIZE=750, $
        TITLE='EBDETECT: Lifetime constraints'
		  PLOT,INDGEN(params.nx),INDGEN(params.ny),POS=[0,0,1,1],XRANGE=[0,nx-1],$
        YRANGE=[0,ny-1],/XS,/YS,/NODATA
    ENDIF
	ENDIF
  ; Create final selection of detections based on lifetime constraints
  detpass = 0L
  IF KEYWORD_SET(params.get_kernels) THEN BEGIN
    totnkernellabels = 0L
    totnkernels = 0L
    sel_kernel_mask = BYTARR(params.nx,params.ny,params.nt)
    last_kernel_detect_counter = 0L
    sum_kernel_detect_counter = 0L
    IF (verbose GE 2) THEN BEGIN
      t_init2 = SYSTIME(/SECONDS)
      feedback_txt2 = 'Getting high-intensity kernels.'
      EBDETECT_FEEDBACK, feedback_txt2+'..', /STATUS
    ENDIF
  ENDIF
	FOR dd=0L,nsel_detections_orig-1 DO BEGIN											
    detlabel = (*detections[sel_detect_idx[dd]]).label
;    print,dd,detpass,detlabel
    ; Compare the detection label with those of the detections to be removed
    IF (nremove_detections GE 1) THEN $
      whereremove = WHERE(params.remove_detections EQ detlabel) $
    ELSE $
      whereremove = -1
    IF (whereremove[0] EQ -1) THEN BEGIN
      ; If detection is not to be removed, continue selecting
  		*sel_detections[detpass] = *detections[sel_detect_idx[dd]]				
  		nt_loc = N_ELEMENTS((*sel_detections[detpass]).t)
      IF KEYWORD_SET(params.get_kernels) THEN kernelresults = PTRARR(nt_loc,/ALLOCATE_HEAP)
  		FOR tt=0,nt_loc-1 DO BEGIN 
        t_real = ((*sel_detections[detpass]).t)[tt]
	  		mask = BYTARR(nx,ny)
	  		mask[(*(*sel_detections[detpass]).det[tt]).pos] = 1B
	  		sel_detect_mask[*,*,((*sel_detections[detpass]).t)[tt]] += mask
        ; Check for kernel pixels within the detection
        IF KEYWORD_SET(params.get_kernels) THEN BEGIN
          kernel_mask = BYTARR(params.nx+2,params.ny+2) 
          loc_detmask = BYTARR(params.nx+2,params.ny+2)
          loc_detmask[1:params.nx,1:params.ny] = mask
          loc_detpos = WHERE(mask EQ 1)
  			  nkernels = 0L
          ; Pixel coordinates to be discarded, init val
				  discard_kernelpix = -1                            
          ; mask_cubes contains mask with 1s & 2s (=kernels)
          tmp_mask = mask_cube[*,*,t_real]     
          tmp_pad_mask = BYTARR(params.nx+2,params.ny+2)
          tmp_pad_mask[1:params.nx,1:params.ny] = tmp_mask
          loc_detpos_pad = WHERE(loc_detmask EQ 1)
          ; Kernels are only those pixels that pass the high intensity threshold
          wherekernelpix = loc_detpos_pad[$
            WHERE(tmp_pad_mask[loc_detpos_pad] EQ 2,nwherekernel)]
          FOR kk=0,nwherekernel-1 DO BEGIN
            ; If the considered pixel is not a discarded pixel, begin growing region
  					IF (TOTAL(discard_kernelpix EQ wherekernelpix[kk]) LE 0) THEN BEGIN					
              IF (nwherekernel GT 1) THEN BEGIN
                ; Grow the region of selected pixels touching the selected pixel
    		        kernelpix = REGION_GROW(tmp_pad_mask,wherekernelpix[kk],/ALL,THRESH=2)	
              ENDIF ELSE kernelpix = wherekernelpix[kk]
  						nkernelpix = N_ELEMENTS(kernelpix)
      					IF (nkernelpix NE 1) THEN $       ; Added check for limb-to-limb
                  discard_kernelpix = kernelpix[1:nkernelpix-1] $
                ELSE $
                  discard_kernelpix = -1
  							kernel_mask[kernelpix] = 1B						; Add the region to the mask
                sel_kernel_mask[*,*,t_real] += kernel_mask[1:nx,1:ny]
  							nkernels += 1L
  							totnkernels += 1L
  					ENDIF
          ENDFOR
  				IF (TOTAL(WHERE(kernel_mask GT 0)) NE -1) THEN BEGIN
  					kernellabels = LABEL_REGION(kernel_mask,/ALL_NEIGHBORS)								
            kernellabels = kernellabels[1:params.nx,1:params.ny]                  
  					nkernellabels = MAX(kernellabels,/NAN)
  					nkernellabels_pix = N_ELEMENTS(WHERE(kernellabels GT 0))							
  					nkernel_pix = N_ELEMENTS(WHERE(kernel_mask GT 0))
  					IF (nkernel_pix NE nkernellabels_pix) THEN BEGIN							
              ; If number of structure pixels != the number of label pixels, stop with error
              EBDETECT_FEEDBACK, /ERROR, $
  						  'Something is very wrong here... '+STRTRIM(nkernel_pix,2)+$
                ' NE '+STRTRIM(nkernellabels_pix,2)
  						STOP
  					ENDIF ELSE BEGIN										
              ; If nlabel_pix = nkernels_pix, then start writing results
  						kernellabel_vals = LINDGEN(nkernellabels)+1
  						kernels = PTRARR(nkernellabels,/ALLOCATE_HEAP)
  						kernels = PTRARR(nkernels,/ALLOCATE_HEAP)
              ; Loop over all labels
  						FOR j=0,nkernellabels-1 DO BEGIN								
                ; Select the pixels corresponding to current label
  							kernelpositions = WHERE(kernellabels EQ kernellabel_vals[j])					
                ; Write results to pointer
  							*kernels[j] = CREATE_STRUCT('label',kernellabel_vals[j]+totnkernellabels,$
                  'pos',kernelpositions)
  						ENDFOR
  					ENDELSE
  				ENDIF ELSE BEGIN
  					nkernellabels = 0
  					kernels = 0
  				ENDELSE
  				totnkernellabels += nkernellabels
          *kernelresults[tt] = CREATE_STRUCT('t',t_real,'nkernels',$
            nkernellabels,'kernels',kernels)					
        ENDIF
	  	ENDFOR
      IF KEYWORD_SET(params.get_kernels) THEN BEGIN
      	pass = 0L
      	totpasses = 0L
      	tt = 0
      	kernel_first_detect = 0
        IF (dd EQ 0) THEN BEGIN
        	WHILE (kernel_first_detect EQ 0) DO BEGIN
        		IF ((*kernelresults[tt]).nkernels GT 0) THEN BEGIN
        			kernel_detect_counter = LONG((*(*kernelresults[tt]).kernels[$
                                          (*kernelresults[tt]).nkernels-1]).label)
        			kernel_first_detect = 1
        		ENDIF
        		tt += 1
        	ENDWHILE
        ENDIF
        t0 = SYSTIME(/SECONDS)
      	FOR t=0L,nt_loc-1 DO BEGIN													
          ; Loop over all detections at the current time step
      		FOR j=0,(*kernelresults[t]).nkernels-1 DO BEGIN									
      			pass += 1L
            ; If the label of the current detection is bigger than the detection counter
      			IF ((*(*kernelresults[t]).kernels[j]).label GT $
                kernel_detect_counter) THEN BEGIN					
              ; Increase the detection counter by 1
      				kernel_detect_counter += 1L										
              ; Relabel the current detection with the updated detection counter
      				(*(*kernelresults[t]).kernels[j]).label = kernel_detect_counter						
      			ENDIF
      			orig_detection = (*(*kernelresults[t]).kernels[j]).pos
      			;;; Check for splitting events ;;;
      			kernel_overlapped = 0
      			ncor = 0          ; number of detections with which there is overlap
      			k_array = -1      ; kernel index array
      			ncomarr = -1      ; array with number of common elements
      			t_usel = t
      			t_ubound = (t + params.t_skip_constraint) < (nt_loc-1)
      			WHILE ((kernel_overlapped EQ 0) AND (t_usel LT t_ubound)) DO BEGIN
      				IF (t_usel LT t_ubound) THEN t_usel += 1
              ; Loop over all detections at the next time step
      				FOR k=0,(*kernelresults[t_usel]).nkernels-1 DO BEGIN								
      					comp_detection = (*(*kernelresults[t_usel]).kernels[k]).pos
                ; Find the common elements between the considered detections
                array_compare = EBDETECT_ARRAY_COMPARE(orig_detection, $
                  comp_detection)
                ; If the number of common elements >= overlap constraint
      					IF ((N_ELEMENTS(array_compare.common_array) GE $
                    params.overlap_constraint) AND $
                    (TOTAL(array_compare.common_array) NE -1)) THEN BEGIN		
						      IF (TOTAL(k_array) NE -1) THEN $
                    k_array = [k_array,k] $
                  ELSE $
                    k_array = k 
						      IF (TOTAL(ncomarr) NE -1) THEN $
                    ncomarr = [ncomarr,array_compare.ncommon_array] $
                  ELSE $
                    ncomarr = array_compare.ncommon_array 
      						kernel_overlapped = 1
      						ncor += 1
      					ENDIF
      				ENDFOR
      			ENDWHILE

            ; If overlap occurs, assign the labels
      			IF kernel_overlapped THEN BEGIN											
              ; overlap
              IF (ncor EQ 1) THEN BEGIN
                ; If there is only one detection overlapping, assign that
                ; detection the current label 
                oldlabel = (*(*kernelresults[t_usel]).kernels[k_array[0]]).label 
                ; Assign the next detection the current detection label
      					(*(*kernelresults[t_usel]).kernels[k_array[0]]).label = $
                  (*(*kernelresults[t]).kernels[j]).label 		
                extraout = 'Overlap: '+STRTRIM(oldlabel,2)+' > '+$
                  STRTRIM((*(*kernelresults[t_usel]).kernels[k_array[0]]).label,2)+','+$
                  STRTRIM(ncomarr[0],2)
              ENDIF ELSE BEGIN
              ; If there is more than one detection overlapping, find out which
              ; one has the max
      					wheremaxoverlap = WHERE(ncomarr EQ MAX(ncomarr,/NAN),$ 
                  COMPLEMENT=wherenotmaxoverlap,NCOMPLEMENT=nwherenotmaxoverlap)
      					oldlabel = (*(*kernelresults[t_usel]).kernels[$
                  k_array[wheremaxoverlap[0]]]).label 
                ; Assign the next detection the current detection label
      					(*(*kernelresults[t_usel]).kernels[k_array[wheremaxoverlap[0]]]).label = $
                  (*(*kernelresults[t]).kernels[j]).label 		
                extraout = 'Overlap: '+STRTRIM(oldlabel,2)+' > '+$
                  STRTRIM((*(*kernelresults[t_usel]).kernels[$
                  k_array[wheremaxoverlap[0]]]).label,2)+','+STRTRIM(ncomarr[0],2)
      					FOR kk=0,nwherenotmaxoverlap-1 DO BEGIN
      						oldlabel = (*(*kernelresults[t_usel]).kernels[$
                    k_array[wherenotmaxoverlap[kk]]]).label 
                  ; - Increase the detection counter by 1
      						kernel_detect_counter += 1L										
                  ; Assign the next detection the current detection label
      						(*(*kernelresults[t_usel]).kernels[$
                    k_array[wherenotmaxoverlap[kk]]]).label = kernel_detect_counter 		
      					ENDFOR
      				ENDELSE
      			ENDIF ELSE $
              extraout = 'No overlap: '+STRTRIM((*(*kernelresults[t]).kernels[j]).label,2)
            IF (verbose GE 1) THEN EBDETECT_TIMER, t+1, nt_loc, t0, $
              EXTRA=extraout, TOTAL_TIME=(verbose GE 2)
      		ENDFOR
      	ENDFOR
;      	last_kernel_detect_counter = kernel_detect_counter
;      	IF (verbose EQ 2) THEN STOP

      	;;; Check for merging events ;;;
      	IF KEYWORD_SET(params.merge_check) THEN BEGIN
          IF (verbose GE 1) THEN BEGIN
            PRINT,''
            IF (verbose GE 2) THEN $
              EBDETECT_FEEDBACK, /STATUS, 'Performing merge check for kernels...'
          ENDIF
      		pass = 0L
      		totpasses = 0L
      		FOR t=0,params.nt-1 DO totpasses += LONG((*kernelresults[t]).nkernels)
      		t0 = SYSTIME(/SECONDS)
          ; Loop over all but the last time step
      		FOR t_dum=0,nt_loc-1 DO BEGIN													
      			t = nt_loc-t_dum-1
            ; Loop over all detections at the current time step
      			FOR j=0,(*kernelresults[t]).nkernels-1 DO BEGIN									
      				pass += 1L
      				orig_detection = (*(*kernelresults[t]).kernels[j]).pos
      				kernel_overlapped = 0
      				ncor = 0
      				k_array = -1
      				ncomarr = -1
      				t_lsel = t
      				t_lbound = (t - params.t_skip_constraint) > 0
              ; Check labels and overlap
      				WHILE ((kernel_overlapped EQ 0) AND (t_lsel GT t_lbound)) DO BEGIN						
      					IF (t_lsel GT t_lbound) THEN t_lsel -= 1
                ; Loop over all detections at the next time step
      					FOR k=0,(*kernelresults[t_lsel]).nkernels-1 DO BEGIN								
      						comp_detection = (*(*kernelresults[t_lsel]).kernels[k]).pos
                  array_compare = EBDETECT_ARRAY_COMPARE(orig_detection, $
                    comp_detection)
                  ; If the number of common elements >= overlap constraint
						      IF ((N_ELEMENTS(array_compare.common_array) GE $
                      params.overlap_constraint) AND $
                      (TOTAL(array_compare.common_array) NE -1)) THEN BEGIN		
						        IF (TOTAL(k_array) NE -1) THEN $
                      k_array = [k_array,k] $
                    ELSE $
                      k_array = k 
						        IF (TOTAL(ncomarr) NE -1) THEN $
                      ncomarr = [ncomarr,array_compare.ncommon_array] $
                    ELSE $
                      ncomarr = array_compare.ncommon_array 
      							kernel_overlapped = 1
      							ncor += 1
      						ENDIF
      					ENDFOR
      				ENDWHILE			

              ; If overlap occured, assign labels
      				IF kernel_overlapped THEN BEGIN										
      					IF (ncor GT 1) THEN BEGIN
      						wheremaxoverlap = WHERE(ncomarr EQ MAX(ncomarr,/NAN))
      						oldlabel = (*(*kernelresults[t]).kernels[j]).label
      						newlabel = (*(*kernelresults[t_lsel]).kernels[$
                    k_array[wheremaxoverlap[0]]]).label 
                  extraout = 'Overlap: '+STRTRIM(oldlabel,2)+' > '+$
                    STRTRIM(newlabel,2)+','+STRTRIM(ncomarr[wheremaxoverlap[0]],2)
                  ;Assign the current detection the previous detection label
      						(*(*kernelresults[t]).kernels[j]).label = newlabel						
      						FOR tt = t+1,nt_loc-1 DO BEGIN
      							kk = 0
      							newlabel_set = 0
      							WHILE ((newlabel_set EQ 0) AND $
                            (kk LT (*kernelresults[tt]).nkernels-1)) DO BEGIN
      								kk += 1	
      								IF ((*(*kernelresults[tt]).kernels[kk]).label EQ oldlabel) THEN BEGIN
      									(*(*kernelresults[tt]).kernels[kk]).label = newlabel
      									newlabel_set = 1
      								ENDIF
      							ENDWHILE
      						ENDFOR
      					ENDIF
      				ENDIF
              IF (verbose GE 1) THEN EBDETECT_TIMER,t_dum+1,nt,t0,EXTRA=extraout, $
                TOTAL_TIME=(verbose GE 2)
      			ENDFOR
      		ENDFOR
      	ENDIF
      ENDIF
      detpass += 1
    ENDIF ELSE BEGIN
      IF (verbose GE 2) THEN $
        EBDETECT_FEEDBACK, /STATUS, 'Removed detection '+STRTRIM(dd,2)+$
          ' with label '+STRTRIM(detlabel,2)
    ENDELSE

    IF KEYWORD_SET(GET_KERNELS) THEN BEGIN
    	; Group kernel detections by label
      new_kernel_counter = kernel_detect_counter - last_kernel_detect_counter
    	kernel_detections = PTRARR(new_kernel_counter,/ALLOCATE_HEAP)
    	sel_kernel_detect_idx = -1
    	t0 = SYSTIME(/SECONDS)
      lifetime_max = 0L
      ; Loop over all single detections
    	FOR d=0L,new_kernel_counter-1 DO BEGIN											
    		t_arr = -1
        t_real_arr = -1
    		j_arr = -1
        label_check = d+1L+last_kernel_detect_counter
        ; Loop over all time steps
    		FOR t=0,nt_loc-1 DO BEGIN												
          ; Loop over all detections at each time step
    			FOR j=0,(*kernelresults[t]).nkernels-1 DO BEGIN								
            ; If the detection label equals the current detection counter
    				IF ((*(*kernelresults[t]).kernels[j]).label EQ label_check) THEN BEGIN					
    					IF (TOTAL(t_arr) NE -1) THEN BEGIN
                t_real_arr = [t_real_arr,(*kernelresults[t]).t]				
                t_arr = [t_arr,t]
              ENDIF ELSE BEGIN
                t_real_arr = (*kernelresults[t]).t 
                t_arr = t
              ENDELSE
    					IF (TOTAL(j_arr) NE -1) THEN $
                j_arr = [j_arr,j] $
              ELSE $
                j_arr = j 
    				ENDIF
    			ENDFOR
    		ENDFOR
    		nt_arr = N_ELEMENTS(t_real_arr)
        ; Determine lifetime
    		kernel_lifetime = t_real_arr[nt_arr-1] - t_real_arr[0] + 1									
        IF (kernel_lifetime GT kernel_lifetime_max) THEN $
          kernel_lifetime_max = kernel_lifetime
    		kernel_det = PTRARR(nt_arr,/ALLOCATE_HEAP)
        ; Loop over all time steps where detection is present
    		FOR tt=0L,nt_arr-1 DO BEGIN											
    			*kernel_det[tt] = $
            CREATE_STRUCT('pos',(*(*kernelresults[t_arr[tt]]).kernels[j_arr[tt]]).pos)
    		ENDFOR
        ; Write results grouped by detection with lifetime information
    		*kernel_detections[d] = CREATE_STRUCT('label',label_check,$
          't',t_real_arr,'lifetime',kernel_lifetime, 'det',kernel_det)				
        IF (verbose GE 1) THEN $
    		EBDETECT_TIMER, label_check, detect_counter, t0, $
          EXTRA='d='+STRTRIM(d,2)+', nt='+$
          STRTRIM(nt_arr,2)+', t_upp='+STRTRIM(t_arr[nt_arr-1],2)+$
          ', t_low='+STRTRIM(t_arr[0],2)+$
          ', t='+STRTRIM(kernel_lifetime,2)+'. So far t_max='+$
          STRTRIM(kernel_lifetime_max,2), TOTAL_TIME=(verbose GE 2)
    	ENDFOR
      ; Add kernel detections to selected detections
  	  IF KEYWORD_SET(params.get_kernels) THEN BEGIN
        *sel_detections[dd] = CREATE_STRUCT(*sel_detections[dd], 'kernel_detections', $
                kernel_detections)
      ENDIF
      last_kernel_detect_counter = kernel_detect_counter
      sum_kernel_detect_counter += new_kernel_counter
    ENDIF
	ENDFOR

	IF ((verbose GE 2) AND KEYWORD_SET(params.get_kernels)) THEN BEGIN
    EBDETECT_FEEDBACK, feedback_txt2, /STATUS, /DONE, T_INIT=t_init2
    EBDETECT_FEEDBACK, /STATUS, $
      'Final number of single kernel detections: '+STRTRIM(kernel_detect_counter,2)
  ENDIF
  replay = 0
  replay_point:
  off = [-5,-20]  ; Offset for label overlays
  ; Create mask from final selection of detections
	FOR t=0,params.nt-1 DO BEGIN													
		IF ((verbose GE 2) AND ~KEYWORD_SET(NO_PLOT)) THEN BEGIN
	    TV,CONGRID(BYTSCL(LP_GET(sum_cube,t),/NAN),750*dataratio,750)
			LOADCT,13,/SILENT
			CONTOUR,CONGRID(sel_detect_mask[*,*,t],750*dataratio,750),COLOR=255, $
        LEVELS=1, /ISOTROPIC, XS=13,YS=13,POSITION=[0,0,1,1], /NORMAL, /NOERASE
      IF comparison_mask_exists THEN $
        CONTOUR,REFORM(LP_GET(paramms.comparison_mask,t)), COLOR=200, LEVELS=1, /ISO, $
          XS=13, YS=13, POS=[0,0,1,1], /NORMAL, /NOERASE
			LOADCT,0,/SILENT
			FOR dd=0,nsel_detections-1 DO BEGIN
				where_detect = WHERE((*sel_detections[dd]).t EQ t, nwhere_detect)
				IF ((TOTAL(where_detect) NE -1) AND (nwhere_detect EQ 1)) THEN BEGIN
					xyout_pos = ARRAY_INDICES(mask,$
            ((*(*sel_detections[dd]).det[where_detect[0]]).pos)[0])+off
					XYOUTS,xyout_pos[0]/FLOAT(nx),xyout_pos[1]/FLOAT(ny),$
            STRTRIM((*sel_detections[dd]).label,2),$
            COLOR=255, /NORMAL, CHARSIZE=2
				ENDIF
			ENDFOR
			XYOUTS,10.,10.,'Selected detections at t='+STRTRIM(t,2),/DATA,COLOR=255,CHARSIZE=2
		  IF (verbose EQ 2) THEN WAIT,0.5
		ENDIF
	ENDFOR
  IF (verbose GE 2) THEN BEGIN
    EBDETECT_FEEDBACK, feedback_txt1, /STATUS, /DONE, T_INIT=t_init
	  IF (verbose EQ 3) THEN BEGIN
      ans = ''
      READ, ans, PROMPT='Do you wish to replay the movie? [Y/N]'
      IF (STRUPCASE(ans) EQ 'Y') THEN GOTO,replay_point
      STOP
    ENDIF
  ENDIF


;================================================================================
;============================== Write final results =============================
;================================================================================
	IF KEYWORD_SET(params.write_detect_final) THEN BEGIN
    IF (verbose GE 2) THEN BEGIN
      t_init = SYSTIME(/SECONDS)
      feedback_txt = 'Writing final results.'
      EBDETECT_FEEDBACK, feedback_txt+'..', /STATUS
    ENDIF
    ; Write IDL save file
		SAVE,sel_detections,nsel_detections,filename=params.outputdir+detect_final_idlsave
		EBDETECT_FEEDBACK, '> Written: '+params.outputdir+detect_final_idlsave
    IF KEYWORD_SET(params.write_mask) THEN BEGIN
  		LP_WRITE,sel_detect_mask, params.outputdir+final_mask_filename
  		EBDETECT_FEEDBACK,'> Written: '+params.outputdir+final_mask_filename
      IF KEYWORD_SET(params.get_kernels) THEN BEGIN
  	  	LP_WRITE,sel_kernel_mask, params.outputdir+kernel_mask_filename
  	  	EBDETECT_FEEDBACK, '> Written: '+params.outputdir+kernel_mask_filename
      ENDIF
    ENDIF
    IF (verbose GE 2) THEN $
      EBDETECT_FEEDBACK, feedback_txt, /STATUS, /DONE, T_INIT=t_init
	ENDIF

;================================================================================
;=========================== Output final statistics ============================
;================================================================================
  EBDETECT_FEEDBACK, 'Detection statistics:'
  EBDETECT_FEEDBACK, '# after intensity & size thresholds: '+$
    STRTRIM(totnlabels,2)
	EBDETECT_FEEDBACK, '# after continuity constraints:      '+$
    STRTRIM(detect_counter,2)
	EBDETECT_FEEDBACK, '# after lifetime constraint:         '+$
    STRTRIM(N_ELEMENTS(sel_detect_idx),2)
	IF KEYWORD_SET(GET_KERNELS) THEN $
    EBDETECT_FEEDBACK, '# of kernels:                        '+$
      STRTRIM(kernel_detect_counter,2)

  EBDETECT_FEEDBACK, /DONE, T_INIT=t_init_overall
  IF ((verbose GE 2) AND ~KEYWORD_SET(NO_PLOT)) THEN WDELETE, 0
	IF (verbose EQ 3) THEN STOP

END
