;+
;
; NAME
;      mc_mgspec1d
;
; PURPOSE
;      Compute Monte Carlo error bars on velocity and EW.  Add to exisitin mgzspec.
;
; SYNTAX
;      mg_domask_spec1d, [ result, files=files, /doplot, /debug, allresult= ]
;
; INPUTS
;      files = an optional parameter giving a vector of spec1d file
;              names. 
;
; KEYWORDS
;
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
;      remove_telluric
;
; EXAMPLES
;
;
; COMMENTS
;   
;
; HISTORY
;      Created August 7, 2002 by mcc.
;      various revisions 30sep02 MD
;      modified from DEEP reduce1d for stellar spectra MG 3/06
;      modified from mg_spec1d
;-
pro sg_mcerrs, maskname,plot=plot,silent=silent

  c=2.99792e5


   ; SET THIS PATH TO DEIMOS DATA WORKING DIRECTORY
;     loadct,12
     data_path = getenv('DEIMOS_DATA')

   ; GET DIRECTORIES AND NAME OF EIGEN TEMPLATES
     sgidl_path = getenv('SGIDL')
     eigendir = sgidl_path + '/templates/'


   ; READ IN ALLDATA FILE
     file = data_path+'/zresult/zspec.'+maskname+'.fits'
     print,'Reading '+file
     data = mrdfits(file,1)
    nslit = n_elements(data)

     nmc = 500

     mctmp = {objname:' ',$
              mc_vcorr:fltarr(nmc),$
              mc_aband:fltarr(nmc),$ 
              mcerr_vcorr:0.0,$
              mc_v0:  0.0,$     
              mcstrct = replicate(mctmp, nslit)


     eigenfile_star = 'deimos-052715.fits'

     shdr = headfits(djs_filepath(eigenfile_star, root_dir=eigendir))
     nstar = sxpar(shdr, 'NAXIS2') > 1
     tmp_subclass = strarr(nstar)   ;types of stars
     for istar=0, nstar-1 do tmp_subclass[istar] = $
               strtrim( sxpar(shdr, 'NAME'+strtrim(string(istar),2)), 2)
  

    ;; Create the wavelength header    
       npix = strcompress(sxpar(shdr, 'NAXIS1'),/rem)
       wav0 = strcompress(sxpar(shdr, 'COEFF0'), /rem)
       dlam = strcompress(sxpar(shdr, 'COEFF1'), /rem)
       mkhdr, head, 1, npix
       sxaddpar, head, 'COEFF0', wav0
       sxaddpar, head, 'COEFF1', dlam




  ; LOOP THROUGH EACH SPECTRUM
    cnt =0
    for islit = 0L,nslit -1 do begin


       mcstrct[islit].objname = data[islit].objname
       print,islit
       if (data[islit].zquality le 2 or data[islit].z*3e5 ge -100) then goto,skipslit
 
       ; IF THIS SLIT IS MISSING SKIP TO END
         zfile1d = strcompress(data_path + data[islit].zspec1d_file,/remove_all)
         sdata = mrdfits(zfile1d,1,/silent)
         if size(sdata, /tname) ne 'STRUCT' then begin

            help,sdata[islit],/st
            stop            
            goto, skipslit

        endif
      
        zfile1d = strcompress(data_path + data[islit].zspec1d_file,/remove_all)
        for j=0L,nmc do begin

            ; READ 1D SPEC FILE
            sdata = mrdfits(zfile1d,1,/silent)
            wv_template = 10^(wav0 + findgen(npix)*dlam)
          


 
           seed=j+islit
           mc_spec = sdata.spec + randomn(seed, npix)*sqrt(1./sdata.ivar)

           q1=where(NOT Float( Finite(sqrt(1./sdata.ivar))))
           q2=where(sqrt(1./sdata.ivar) ge 1e10)
           q=[q1,q2]
           mc_spec[q]=0.
           mc_ivar = sdata.ivar


           if (j eq 0) then mc_spec = sdata.spec


        ; A-BAND CORRECTION
        zans_aband  = sg_aband(wv_template,mc_spec,mc_ivar,head,plot=plot)


       ; REMOVE TELLURIC LINES FROM CROSS CORRELATION
         sky_min1 = [6862, 7588] * (1. + zans_aband.z)  
         sky_max1 = [6950, 7700] * (1. + zans_aband.z)  

         sky_min2 = [7167, 8210, 8900] * (1. + zans_aband.z)  
         sky_max2 = [7315, 8325, 9200] * (1. + zans_aband.z)  


        q1=where(wv_template ge sky_min1[0] and wv_template le sky_max1[0])
        q3=where(wv_template ge sky_min1[1] and wv_template le sky_max1[1])

        q=[q1,q3]
        if (q[0] ne -1) then mc_ivar[q] =mc_ivar[q]*0.0001   ;  THE STRONG A- B-BANDS

        q2=where(wv_template ge sky_min2[0] and wv_template le sky_max2[0])
        q4=where(wv_template ge sky_min2[1] and wv_template le sky_max2[1])
        q5=where(wv_template ge sky_min2[2] and wv_template le sky_max2[2])

        q=[q2,q4,q5]
        if (q[0] ne -1) then mc_ivar[q] =mc_ivar[q]*0.1   ; THERE ARE THE WEAKER ABS LINES
        




; _____________________________________________________________
; MEASURE REDSHIFTS
; _____________________________________________________________

     ; DO STELLAR REDSHIFTS 
        npoly = data[islit].npoly
        zmin = -0.0033          ; -1000 km/sec
        zmax = 0.0033           ; +1000 km/sec
        pspace = 0.25
        nfind = 1



      ; RUN ERRORS ONLY ON BEST FITTING STAR IN ORIGINAL FIT 
        qst = where(strcompress(data[islit].subclass,/rem) eq tmp_subclass)
        istar = qst

            subclass = strtrim( sxpar(shdr, 'NAME'+strtrim(string(istar),2)), 2)       
             

            zans_star = x_zfind(mc_spec,mc_ivar, hdr=head, eigendir=eigendir,$
                                  eigenfile=eigenfile_star, columns=istar, $
                                  npoly=npoly, zmin=zmin, zmax=zmax, pspace=1, $       
                                   nfind=nfind, width=5*pspace, $
                                   objflux=objflux, objivar=objivar,loglam=loglam, $
                                  debug=debug, /silent)

        mc_z = zans_star.z

        
        if (j eq 0) then begin
           z0     = zans_star.z
           a0     = zans_aband.z
           vcorr0 = c*(z0 - a0)
           mcstrct[islit].mc_v0  = vcorr0
        endif


        if (j ne 0) then begin
           mcv = c*(mc_z - zans_aband.z) 
           mcstrct[islit].mc_vcorr[j-1] = mcv 
           mcstrct[islit].mc_aband[j-1] = zans_aband.z

           print,islit,'  ',j,mcv
;           print,'MC Velocity = '+string(c*(mc_z - zans_aband.z),format='(f7.2)')
;           print,'   Velocity = '+string(vcorr0,format='(f7.2)')
;           print,'   Original = '+string(c*(data[islit].z - data[islit].aband),format='(f7.2)')
         
       endif    


     endfor  ; END MC LOOP ON SLIT


      ; CALCULATE RMS ERRORs
        mcstrct[islit].mcerr_vcorr= sqrt( total((mcstrct[islit].mc_vcorr - $
                                                       vcorr0)^2)/float(nmc))

       
        print,'MC ERROR    = '+string(mcstrct[islit].mcerr_vcorr,format='(f7.2)')
        print,'ZFIND ERROR = '+string(data[islit].z_err*c,format='(f7.2)')


        cnt=cnt+1
        if (cnt gt 10) then begin
          mwrfits, mcstrct, data_path+'/zresult/mcall_'+maskname+'.fits', /create
          cnt=0
        endif

        skipslit:a=1


     endfor  ; END MASK

    spawn,'rm '+ data_path+'/zresult/mcall_'+maskname+'.fits'
    mwrfits, mcstrct, data_path+'/zresult/mcall_'+maskname+'.fits', /create
    spawn,'gzip -f '+ data_path+'/zresult/mcall_'+maskname+'.fits'

end








