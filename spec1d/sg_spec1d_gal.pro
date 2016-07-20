;+
;
; NAME
;      mg_galspec
;
; PURPOSE
;      A stand alone code to determine redshifts of galaxies/QSOs
;      from stellar DEIMOS masks-- only runs on Q=1 targets
;
; SYNTAX
;      mg_galspec, maskname
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
;    mg 7/07 - modified from mg_spec1d for galaxies only
;    mg 7/06 - need to fix spec1d/spec2d names
;    mg 2/07 - changed telluric corrections
;
; HISTORY
;      Created August 7, 2002 by mcc.
;      various revisions 30sep02 MD
;      modified from DEEP reduce1d for stellar spectra MG 3/06
;-
pro sg_spec1d_gal, maskname, plot=plot


   ; SET THIS PATH TO DEIMOS DATA WORKING DIRECTORY
     loadct,12
     data_path = getenv('DEIMOS_DATA')
                
     
   ; READ ZSPEC FILE
     zfile = data_path + '/zresult/mgzresult.'+maskname+'.fits'
     zdata = xmrdfits(zfile,1)


   ; GET DIRECTORIES AND NAME OF EIGEN TEMPLATES
     mgidl_path = getenv('SGIDL')
     eigendir = mgidl_path + '/templates/'

     eigenfile_gal =  'deimos-galaxy.fits'
     
     shdr = headfits(djs_filepath(eigenfile_gal, root_dir=eigendir))
     nstar = sxpar(shdr, 'NAXIS2') > 1

    ;; Create the wavelength header    
       npix = strcompress(sxpar(shdr, 'NAXIS1'),/rem)
       wav0 = strcompress(sxpar(shdr, 'COEFF0'), /rem)
       dlam = strcompress(sxpar(shdr, 'COEFF1'), /rem)
       mkhdr, head, 1, npix
       sxaddpar, head, 'COEFF0', wav0
       sxaddpar, head, 'COEFF1', dlam
       wv_template = 10^(wav0 + findgen(npix)*dlam)


;--------------------------------------------------
; loop over all of the spec1d files.

 ; DETERMINE PRIMARY SLITS  
   objname = strcompress(zdata.objname,/rem)
   for ii=0,n_elements(objname)-1 do begin


     zfile = strcompress(data_path + zdata[ii].zspec1d_file,/rem)
     if (strcompress(findfile(zfile),/rem) eq '') then goto, skipslit
     print,zfile
     zspec = mrdfits(zfile,1,hdr)

      ; REBIN DATA ONTO THE SAME WAVELENGTH ARRY AS GALAXY-- NOTE
      ;  THIS IS A SCALE THAN STARS!
        x_specrebin, zspec.lambda, zspec.spec, wv_template, flux
        x_specrebin, zspec.lambda, zspec.ivar, wv_template, ivar


          q=where(NOT float(finite(ivar)))
          if (q[0] ne -1) then ivar[q] = 0 
          if (q[0] ne -1) then flux[q] = 0


       ; REMOVE TELLURIC LINES FROM CROSS CORRELATION
         sky_min1 = [6862, 7588] 
         sky_max1 = [6950, 7700] 

         sky_min2 = [7167, 8210, 8900]
         sky_max2 = [7315, 8325, 9200]


        q1=where(wv_template ge sky_min1[0] and wv_template le sky_max1[0])
        q3=where(wv_template ge sky_min1[1] and wv_template le sky_max1[1])

        q=[q1,q3]
        if (q[0] ne -1) then ivar[q] =ivar[q]*0.0001   ;  THE STRONG A- B-BANDS

        q2=where(wv_template ge sky_min2[0] and wv_template le sky_max2[0])
        q4=where(wv_template ge sky_min2[1] and wv_template le sky_max2[1])
        q5=where(wv_template ge sky_min2[2] and wv_template le sky_max2[2])

        q=[q2,q4,q5]
        if (q[0] ne -1) then ivar[q] =ivar[q]*0.1   ; THERE ARE THE WEAKER ABS LINES
        

; _____________________________________________________________
; MEASURE REDSHIFTS
; _____________________________________________________________

     ; DO STELLAR REDSHIFTS 
        npoly=3
        zmin = 0.0033
;        zmin = 0.0001
        zmax = 3.
        pspace = 1
        nfind = 1


            tmp_zans_gal = mg_zfind(flux,ivar, hdr=head, eigendir=eigendir,$
                                  eigenfile=eigenfile_gal,$
                                  npoly=npoly, zmin=zmin, zmax=zmax, pspace=15,$       
                                   nfind=nfind, width=5*pspace, $
                                  debug=debug, /silent)


            print,tmp_zans_gal.z
             

            if keyword_set(plot) then begin 
;                plot,wv_template,flux,xrange=[6500,8900]
;                synth = synthspec(tmp_zans_gal, loglam=alog10(wv_template),$
;                   eigendir = eigendir)
;                oplot,wv_template,synth,color=180
            endif

            print,tmp_zans_gal.rchi2,zdata[ii].rchi2
            

            galaxy = 0
            if (tmp_zans_gal.rchi2  le  0.95*zdata[ii].rchi2 and $
                tmp_zans_gal.rchi2 ge 0.5) then galaxy=1

            if (galaxy eq 1) then begin
               
               zdata[ii].class = 'GALAXY'
               zdata[ii].subclass = 'GALAXY'
               zdata[ii].tfile = eigenfile_gal
               zdata[ii].z     = tmp_zans_gal.z
               zdata[ii].z_err = tmp_zans_gal.z_err
               zdata[ii].dof   = tmp_zans_gal.dof
               zdata[ii].rchi2 = tmp_zans_gal.rchi2
               zdata[ii].tcolumn= tmp_zans_gal.tcolumn
               zdata[ii].npoly = tmp_zans_gal.npoly
               zdata[ii].theta = tmp_zans_gal.theta



            endif
            zdata[ii].tmpl_rchi2[16] = tmp_zans_gal.rchi2
            zdata[ii].tmpl_z[16]     = tmp_zans_gal.z
            zdata[ii].tmpl_tfile[16] = 'galaxy'

 
   skipslit:a=1

 endfor 

        ; RE-WRITE OVER SPEC1D FILE
          outfile = data_path + '/zresult/mgzresult.'+maskname+'.fits'
          mwrfits, zdata, outfile, /create
          mwrfits, zdata, outfile
          print,'Writing output to ',outfile

end








