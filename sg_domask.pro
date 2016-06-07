:;+
; NAME:
;    sg_domask
;
; PURPOSE:
;    Wrapped to reduce DEIMOS masks
;
; CALLING SEQUENCE:
;    sg_domask,'mymask.plan'
; 
; INPUTS:
;
; OPTIONAL INPUTS:
;       
; KEYWORDS:
;
; OUTPUTS:

;
; COMMENTS:
;
; REVISION HISTORY:
;   Jun16,  MG modify domask for dSph project
;
;----------------------------------------------------------------------
pro sg_domask, planfile, chiplist=chiplist,$
                wave_offset=wave_offset, tweakmethod=tweakmethod, $
                slit_match_problem=slit_match_problem


  if NOT keyword_set(planfile) then message, 'You must specify a plan file!'

  ; READ DESIGN INFO FROM HEADER, CREATE BINTAB ORG FILE
    make_bintab_file, planfile


; write calibSlit files
  mg_deimos_mask_calibrate, planfile, /noplot, chiplist=chiplist, $
    wave_offset=wave_offset, slit_match_problem=slit_match_problem


; write spSlit files
  jds_deimos_2dreduce, planfile


; combine multiple exposures
  list=findfile('spSlit*.fits')
  jds_spslit_combine,list, nlsky = nlsky

  read_planfile, planfile, maskname, rawdatadir, outdatadir, $
             flatnames

  head_flat = headfits(flatnames[0])
  deimos_grating, head_flat, grating, grangle, lambda_c

  minlambda=lambda_c - 1300*(1200./grating) 
  bluelim=5500.-1500.*(minlambda lt 5500)

  slitfiles = findfile('slit*.fits', count=nfiles)
  slitfiles = slitfiles[sort(slitfiles)]
  isred = (strpos(slitfiles, 'R.fits') ne -1)
  trans = max(where(isred and lindgen(nfiles) le nfiles/2))
  if trans eq -1 then trans = floor(nfiles/2)

; Make simple 2d images, 2 per mask
  epos = strpos(slitfiles[0], '.fits')
  masknumber = strmid(slitfiles[0], 4, epos-4-4) ;get mask name '.xxxx.'
;  image = slit_merge_lambda(slitfiles[0:trans], hdr,blue=bluelim)
;  writefits, 'Allslits0'+masknumber+'fits', image, hdr
;  if trans gt 1 then begin
;    image = slit_merge_lambda(slitfiles[trans+1:nfiles-1], $
;                              hdr, blue = bluelim)
;    writefits, 'Allslits1'+masknumber+'fits', image, hdr
;  endif


; make non-local allslits.
  if keyword_set(nlsky) then begin
    img = slit_merge_lambda(slitfiles[0:trans], hdr, /nonlocal,blue=bluelim)
    writefits, 'Allslits0' + masknumber + 'nonlocal.fits', img, hdr
    if trans gt 1 then begin
      img = slit_merge_lambda(slitfiles[trans+1:nfiles-1], $
                              hdr, /nonlocal, blue = bluelim)
      writefits, 'Allslits1' + masknumber + 'nonlocal.fits', img, hdr
    endif
  endif

; do 1d extraction.
  deimos_isdeep, isdeep, maskname
  slitfiles = findfile('slit*.fits', count=nfiles)

  if isdeep then $
    do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.1 $
  else begin
      if strpos(strupcase(maskname), 'KTRS') ge 0 then $
        do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.5 $
      else do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.1
  endelse

; do 1d extraction w/ non-local sky
  if keyword_set(nlsky) then $
    do_extract, files=slitfiles, /nonlocal

  print
  print, systime()+'  Done with slit processing!!!'
  print
  print, 'Total time elapsed: ', (systime(1)-t1)/3600., ' hours.'

; compress the files.
  spawn, 'gzip -1 -vf slit*.*.fits ' ;operations suspend
  if not keyword_set(keep_calibs) then spawn,'rm calib* spSlit*'

  print
  print, systime()+'  Done with gzipping!!!'
  print
  print, 'Total time elapsed: ', (systime(1)-t1)/3600., ' hours.'


; do quality assurance
  if isdeep then qa_check, /doplot

; make a done-processing file to signify the completetion of the
; spec2d pipeline.
  openw, 2, 'doneprocessing.txt'
  printf, 2, spec2d_version()
  close, 2

  exit
end







