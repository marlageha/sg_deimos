;+ 
; NAME:
;    mg_aband
;
; PURPOSE:
;   Calculate A-band correction based on hot star/target star cross
;   correlation.   Outputs a telluric redshift correction (z_telluric)
;   such that:
;    v_corrected = c * (z_observed - z_telluric)
;
; CALLING SEQUENCE:
;   
;     z_telluruc = mg_aband(lambda, flux, ivar, head)
;
; INPUTS:
;     lambda, flux, ivar -- stellar spectrum.  
;
; OPTIONAL INPUTS:
;      create_template -- creates a-band template from multiple stellar spectra
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   
;
; REVISION HISTORY:
;   9/06 MG
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function sg_aband, lambda,flux,ivar,head,create_template = create_template,$
                   plot=plot

    ; GET DIRECTORIES AND NAME OF EIGEN TEMPLATES
     mgidl_path = getenv('SGIDL')
     eigendir = mgidl_path + '/templates'
     aband_star = 'deimos-aband.fits'

    ; CREATE A-BAND TEMPLATE
      if keyword_set(create_template) then begin

          mk_aband_template

     endif
    
     ahead = headfits(djs_filepath(aband_star, root_dir=eigendir))

        q1=where(lambda ge 6866 and lambda le 6912)
        q2=where(lambda ge 7167 and lambda le 7320)
        q3=where(lambda ge 7593 and lambda le 7690)
        q4=where(lambda ge 8110 and lambda le 8320)
;        q5=where(lambda ge 8925 and lambda le 9120)  ; this region has problems
        q=[q1,q2,q3,q4]

        tmp = ivar
        tmp[0:*] = 0
        tmp[q] = 1
        qtmp = where(tmp ne 1)
        
        ivar_sky = ivar
        ivar_sky[qtmp]=ivar[qtmp]*0.0001  ; supress other wavelengths


 
   ; DO STELLAR REDSHIFTS 
        npoly=1
        zmin = -0.0003         ; -100 km/sec
        zmax = 0.0003          ; +100 km/sec
        pspace = 0.25
        nfind = 1

        zans_aband = x_zfind(flux,ivar_sky, hdr=head, eigendir=eigendir,$
                             eigenfile=aband_star,columns=0,$
                             npoly=npoly, zmin=zmin, zmax=zmax, pspace=1, $
                              nfind=nfind, width=5*pspace, $
                              objflux=objflux, objivar=objivar,loglam=loglam,$
                              debug=debug, /silent)


        ;vel_aband =zans_aband.z*2.9979e5
        ;print,vel_aband

  ; PLOT RESULTS
    if keyword_set(plot) then begin
        plot,lambda,flux,xrange=[7500,7800]
        synth = synthspec(zans_aband, loglam=alog10(lambda),$
                          eigendir = eigendir)
        oplot,lambda,synth,color=100
    endif

return,zans_aband

end
