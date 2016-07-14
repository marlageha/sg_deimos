;+
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
  sg_make_bintab_file, planfile


; write calibSlit files
  sg_deimos_mask_calibrate, planfile, chiplist=chiplist, $
    wave_offset=wave_offset, slit_match_problem=slit_match_problem


; write spSlit files
  sg_deimos_2dreduce, planfile,tweakmethod=tweakmethod,chiplist=chiplist



; combine multiple exposures
  list=findfile('spSlit*.fits')

  jds_spslit_combine,list

  sg_read_planfile, planfile, maskname, rawdatadir, outdatadir, $
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

; Make combined slit image
  epos = strpos(slitfiles[0], '.fits')
  masknumber = strmid(slitfiles[0], 4, epos-4-4) ;get mask name '.xxxx.'
  image = slit_merge_lambda(slitfiles[0:trans], hdr,blue=bluelim)
  writefits, 'Allslits'+masknumber+'fits', image, hdr
 

; do 1d extraction.  Why these nsigmas?
  isdeep =0
  slitfiles = findfile('slit*.fits', count=nfiles)
  do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.1



  print
  print, systime()+'  Done with slit processing!!!'

; compress the files.
  spawn, 'gzip -1 -vf *.fits ' 
;  if not keyword_set(keep_calibs) then spawn,'rm calib* spSlit*'

 ; Run Spec1d
  sg_spec1d,maskname
  print,'Now run sg_zspec,'+maskname

  exit
end







