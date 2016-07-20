;+
;
; NAME
;      mg_spec1d
;
; PURPOSE
;      A wrapper routine to organize and run the 1d analysis
;      pipeline tuned for stellar data.
;
; SYNTAX
;      mg_spec1d, maskname
;
; INPUTS
;      maskname -- name of mask directory 'N147-1'
;
; KEYWORDS
;
; OUTPUTS
;     result -- output structure giving result of last object in call
;               (or single object)
;   
; PROCEDURES CALLED 
;      splog
;      zfin 
;      zrefind
;      vdispfit
;      fill_gap
;      cpbackup
;      dfpsplot
;      concat_dir
;      
;
; EXAMPLES
;
;
; COMMENTS
;    mg 7/06 = need to fix spec1d/spec2d names
;    mg 2/07 - changed telluric corrections
;    mg 7/09 - ADDED BINTAB file to determine alignment stars
;
; HISTORY
;      Created August 7, 2002 by mcc.
;      various revisions 30sep02 MD
;      modified from DEEP reduce1d for stellar spectra MG 3/06
;-
pro sg_spec1d, maskname, allchips=allchips,plot=plot,nozfile=nozfile,nogalaxy=nogalaxy


   ; SET THIS PATH TO DEIMOS DATA WORKING DIRECTORY
     loadct,12
     data_path = getenv('DEIMOS_DATA')


   ; FIND DIRECTORIES FOR ALL CCD CHIPs (assumes /maskname_ccdnumber)  
     files = findfile(data_path + maskname +'/spec1d.*.fits*')
     if keyword_set(allchips) then files = findfile(data_path + maskname +$
                     '_*/spec1d.*.fits*')
     nfiles = n_elements(files)


   ; CHECK THAT SPEC1D FILES ARE FOUND
     if nfiles eq 0 then message, 'No 1d files specified by user ' + $
     'or found in current directory!'

   ; READ BINTABS FILE TO SET ALIGNMENT STARS
     bintabfile = strcompress(data_path+maskname+'/'+$
                              maskname + '.bintabs.fits.gz',/rem)
     bintab = mrdfits(bintabfile,1)
     

   ; SET OUTPUT DIRECTOY 
     outdir = data_path + '/zresult/'
     zdir = findfile(outdir)
     if (zdir[0] eq '') then spawn,'mkdir '+outdir
     outfile = data_path + '/zresult/mgzresult.'+maskname+'.fits'


   ; GET DIRECTORIES AND NAME OF EIGEN TEMPLATES
     sgidl_path = getenv('SGIDL')
     eigendir = sgidl_path + '/templates/'

;     eigenfile_star = 'deimos-021507.fits'
     eigenfile_star = 'deimos-052715.fits'
     eigenfile_gal =  'deimos-galaxy.fits'
     
     shdr = headfits(djs_filepath(eigenfile_star, root_dir=eigendir))
     nstar = sxpar(shdr, 'NAXIS2') > 1
     subclass = strarr(nstar)   ;types of stars
     for istar=0, nstar-1 do subclass[istar] = $
               strtrim( sxpar(shdr, 'NAME'+strtrim(string(istar),2)), 2)
  

    ;; Create the wavelength header    
       npix = strcompress(sxpar(shdr, 'NAXIS1'),/rem)
       wav0 = strcompress(sxpar(shdr, 'COEFF0'), /rem)
       dlam = strcompress(sxpar(shdr, 'COEFF1'), /rem)
       mkhdr, head, 1, npix
       sxaddpar, head, 'COEFF0', wav0
       sxaddpar, head, 'COEFF1', dlam
       wv_template = 10^(wav0 + findgen(npix)*dlam)


     ; CREATE STRUCTURE FOR REBINNED 1D SPECTRA
        ztemp = {zstar, spec:fltarr(npix),$
                       lambda:fltarr(npix),$
                       ivar:fltarr(npix),$
                       objpos:0.0,$
                       fwhm:0.0,$
                       nsigma:0.0,$
                       r1:0.0,$
                       r2:0.0}

        result = replicate({sg_spec1dstrct}, nfiles)
        result[0:*].eigendir  =eigendir



;--------------------------------------------------
  ; CALCULATE HELIOCENTRIC CORRECTION
    print,'Reading ',strcompress(files[0],/rem)
    spec = mrdfits(strcompress(files[0],/rem),1,hdr)
    mjd = sxpar(hdr, 'MJD-OBS') 
    jd = double(mjd) + 2400000.5
    ra_obj  = strcompress(sxpar(hdr, 'RA'), /rem)
    dec_obj = strcompress(sxpar(hdr, 'DEC'), /rem)
    sradec,ra_obj,dec_obj,r,d

       if (mjd eq 0) then begin
          print,'ERROR IN HEADER!  No MJD tage. Computing from date/time'
          
           date = sxpar(hdr, 'DATE-OBS') 
           dt = strsplit(date,'-',/extract)
           yr = float(dt[0])
           mn =float(dt[1])
           dy =float(dt[2])

           time = sxpar(hdr, 'UTC')
           tm = strsplit(time,':',/extract)
           hr = float(tm[0])
           mnt = float(tm[1])
           mjd = date2mjd(yr,mn,dy)
           jd = double(mjd) + 2400000.5 + hr/24. + mnt/(24.*60)
           print,jd

       endif
   vh =  helio_deimos(r,d, jd=jd)
   print,ra_obj,dec_obj,r,d,vh
   result[0:*].vhelio=vh


   
;--------------------------------------------------
; loop over all of the spec1d files.

   for ii=0,nfiles-1 do begin

        print,files[ii]
        hdr = headfits(files[ii], ext=1, /silent)
        maskname = strcompress(sxpar(hdr, 'SLMSKNAM'), /rem)
        slitname = strcompress(sxpar(hdr, 'SLITNO'), /rem)
        objname = strcompress(sxpar(hdr, 'OBJNO'), /rem)
        ra_obj  = strcompress(sxpar(hdr, 'RA_OBJ'), /rem)
        dec_obj = strcompress(sxpar(hdr, 'DEC_OBJ'), /rem)
        date = strcompress(sxpar(hdr, 'DATE-OBS'), /rem)
        mjd = double(sxpar(hdr, 'MJD-OBS'))
        airmass = sxpar(hdr, 'AIRMASS')

      ; get the observation date from the header.
        date = strcompress(sxpar(hdr, 'DATE-OBS'), /rem)
        mjd = double(sxpar(hdr, 'MJD-OBS'))


      
      ; fill the gap between the blue and red portions of the 1-d spectra.
        boxsprof = 0
        ss1d = sg_fill_gap(files[ii], boxsprof=boxsprof)



      ; ADD FILE NAMES TO STRCT 
        s=strpos(files[ii],maskname)
        result[ii].spec1d_file = strcompress('/'+strmid(files[ii],s),/rem)


        tmp1 = str_sep(result[ii].spec1d_file,'/')
        tmp2 = str_sep(result[ii].spec1d_file,'.')


        result[ii].spec2d_file = strcompress('/'+tmp1[1]+'/slit.'+tmp2[1]+$
                                             '.'+tmp2[2]+'R.fits.gz',/rem)

        result[ii].zspec1d_file = ' '
        if NOT (keyword_set(nozfile)) then begin
           result[ii].zspec1d_file= strcompress('/'+tmp1[1]+'/z'+tmp1[2],/rem)
            zstarfile= strcompress('zspec1d.'+tmp1[1]+'.'+$
                                             tmp2[2]+'.'+tmp2[3]+'.fits',/rem)
        endif



          result[ii].objname = objname
          result[ii].slitname = slitname
          result[ii].maskname = maskname
          result[ii].ra  = ra_obj
          result[ii].dec = dec_obj
          result[ii].airmass = airmass
          result[ii].date = date
          result[ii].mjd = mjd
          result[ii].comment = ' '

          result[ii].sn = 0

        ; INITIALIZE ZQUALITY, SET ALIGNMENT STARS = -1
          result[ii].zquality = -2  
          q=where(strcompress(objname,/rem) eq strcompress(bintab.object,/rem))
          if strcompress(bintab[q].objclass,/rem) eq 'Alignment_Star' then $
             result[ii].zquality = -1


      ; IF THIS SLIT IS MISSING SKIP TO END
       if size(ss1d, /tname) ne 'STRUCT' then begin
          ss1d=mg_fill_gap(files[ii])
       endif
       if size(ss1d, /tname) ne 'STRUCT' then begin
          print, 'Skipping slit ' + slitname + '!'
          if (strcompress(result[ii].objname,/rem) ne 'serendip1' and $
              strcompress(result[ii].objname,/rem) ne 'serendip2') then stop
          goto, skipslit
      endif 

       ; SKIP SLIT IF NO DATA
          q = where( NOT Float( Finite(ss1d.spec*sqrt(ss1d.ivar))),nq)
          if (nq ge 50) then begin
              print,'too many NaNs'
              stop
              goto,skipslit
           endif
          if (q[0] ne -1) then ss1d.ivar[q] = 0


      ; FIX MISMATCH BETWEEN SLITS
        if (n_elements(ss1d.spec) ge 4200) then begin   
          m1=median(ss1d.spec[3950:4050])
          m2=median(ss1d.spec[4150:4250])
          m=m1/m2
          if (m le 0) then m=1.
          if (m le 0.7 or m ge 1.3) then m=1.
          ss1d.spec[4100:*] = ss1d.spec[4100:*] *m
          ss1d.ivar[4100:*] = ss1d.ivar[4100:*] /m^2
          print,m
        endif



          
      ; REBIN DATA ONTO THE SAME WAVELENGTH ARRY
          lambda = float(ss1d.lambda) 
          mgx_specrebin, lambda, ss1d.spec, wv_template, flux
          mgx_specrebin, lambda, ss1d.ivar, wv_template, ivar
          q=where(NOT float(finite(ivar)))
          if (q[0] ne -1) then ivar[q] = 0
          if (q[0] ne -1) then flux[q] = 0


       ; A-Band correction
         sivar = ivar 
         zans_aband = sg_aband(wv_template,flux,ivar,head,plot=plot)
         result[ii].aband = zans_aband.z

       ; REMOVE TELLURIC LINES FROM CROSS CORRELATION
         sky_min1 = [6862, 7588] * (1. + zans_aband.z)  
         sky_max1 = [6950, 7700] * (1. + zans_aband.z)  

         sky_min2 = [7167, 8210, 8900] * (1. + zans_aband.z)  
         sky_max2 = [7315, 8325, 9200] * (1. + zans_aband.z)  


        q1=where(wv_template ge sky_min1[0] and wv_template le sky_max1[0])
        q3=where(wv_template ge sky_min1[1] and wv_template le sky_max1[1])

        q=[q1,q3]
        if (q[0] ne -1) then ivar[q] =ivar[q]*0.0001   ;  THE STRONG A- B-BANDS

        q2=where(wv_template ge sky_min2[0] and wv_template le sky_max2[0])
        q4=where(wv_template ge sky_min2[1] and wv_template le sky_max2[1])
        q5=where(wv_template ge sky_min2[2] and wv_template le sky_max2[2])

        q=[q2,q4,q5]
        if (q[0] ne -1) then ivar[q] =ivar[q]*0.1   ; THERE ARE THE WEAKER ABS LINES
        
;qtest = where(lambda ge 7700)
;ivar[qtest] = ivar[qtest]*0.0001

       ; ESTIMATE S/N NEAR Ca II LINES
          s2n=0
          qsn = where(lambda ge 8000 and lambda le 8700)
          if (qsn[0] ne -1) then s2n = mean(ss1d.spec[qsn]*sqrt(ss1d.ivar[qsn]))
          if (qsn[0] eq -1) then s2n = 0
          splog, ' '
          splog, 'mean S2N for frame: ', files[ii], ' : ', s2n
          result[ii].SN = s2n



       ; WRITE ZSTAR FILE-- REBINNED 1D SPECTRUM
         zstar = replicate(ztemp,1)
         zstar.spec   = flux
         zstar.ivar   = sivar
         zstar.lambda = wv_template
         zstar.objpos = ss1d.objpos
         zstar.fwhm   = ss1d.fwhm
         zstar.nsigma = ss1d.nsigma
         zstar.r1     =ss1d.r1
         zstar.r2     =ss1d.r2

        if NOT (keyword_set(nozfile)) then begin
           print
           print,'Writing rebinned file', result[ii].zspec1d_file
           mwrfits,zstar,zstarfile,/create
           spawn,'gzip -f '+zstarfile
        endif
; _____________________________________________________________
; MEASURE REDSHIFTS
; _____________________________________________________________

     ; DO STELLAR REDSHIFTS 
        npoly=4
        ; SET LOWER POLY ORDER IF MISSING EITHER RED OR BLUE SPECTRUM
        if (min(ss1d.lambda) ge 7000) then npoly =2   ; missing blue
        if (max(ss1d.lambda) le 8100) then npoly =2   ; missing red
        print,minmax(ss1d.lambda)
        print,'Npoly = ',npoly

        zmin = -0.0033          ; -1000 km/sec
        zmax = 0.0033           ; +1000 km/sec
        pspace = 0.25
        nfind = 1


        zans_star=0
        for istar=0, nstar-1 do begin
            subclass = strtrim( sxpar(shdr, 'NAME'+strtrim(string(istar),2)), 2)       
             

            tmp_zans_star = mg_zfind(flux,ivar, hdr=head, eigendir=eigendir,$
                                  eigenfile=eigenfile_star, columns=istar, $
                                  npoly=npoly, zmin=zmin, zmax=zmax, pspace=1, $       
                                   nfind=nfind, width=5*pspace, $
                                   objflux=objflux, objivar=objivar,loglam=loglam, $
                                  debug=debug, /silent)

            print,subclass, ' ',tmp_zans_star.rchi2,tmp_zans_star.z*3e5
        
            tmp_zans_star.class = 'STAR'
            tmp_zans_star.subclass = subclass

             
            if(n_tags(zans_star) gt 0) then $
              zans_star=[zans_star, tmp_zans_star] $
            else $
              zans_star=tmp_zans_star
 
;            if (istar eq 0) then plot,wv_template,flux,xrange=[8300,8700]
;            synth = synthspec(tmp_zans_star, loglam=alog10(wv_template),$
;             eigendir = eigendir)
;            oplot,wv_template,synth,color=180


        endfor
        isort     = sort(zans_star.rchi2 + (zans_star.dof EQ 0)* $
                     max(zans_star.rchi2))

         result[ii].tmpl_rchi2 = zans_star[isort].rchi2
         result[ii].tmpl_z     = zans_star[isort].z
         result[ii].tmpl_tfile = zans_star[isort].subclass
         for jj=istar,34 do result[ii].tmpl_tfile[jj] = ' '


        result[ii].class='STAR'
        result[ii].subclass=zans_star[isort[0]].subclass
        result[ii].z       =zans_star[isort[0]].z
        result[ii].z_err   =zans_star[isort[0]].z_err
        result[ii].vdisp   =zans_star[isort[0]].vdisp
        result[ii].vdisp_err=zans_star[isort[0]].vdisp_err
        result[ii].rchi2   =zans_star[isort[0]].rchi2
        result[ii].dof     =zans_star[isort[0]].dof
        result[ii].tfile   =zans_star[isort[0]].tfile
        result[ii].tcolumn =zans_star[isort[0]].tcolumn
        result[ii].npoly   =zans_star[isort[0]].npoly
        result[ii].theta   =zans_star[isort[0]].theta



if keyword_set(plot) then begin
;      plot,wv_template,flux,xrange=[8100,8700]
;      synth = synthspec(result[ii], loglam=alog10(wv_template),$
;                        eigendir = eigendir)
;      oplot,wv_template,synth,color=180
endif

      print,'v = ',3e5 * result[ii].z,'  ',result[ii].subclass
      print,'SN = ',s2n


      skipslit:a=1

  endfor 

  ; WRITE OUTPUT

          mwrfits, result, outfile, /create
          mwrfits, result, outfile
          print,'Writing output to ',outfile

          
       ; RUN GALAXIES

        ;if NOT keyword_set(nogalaxy) then 
        sg_spec1d_gal,maskname,plot=plot

end








