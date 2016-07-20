;--------------------------------------------------------------------
;+
; NAME:
;   mg_zspec
; 
; PURPOSE:
;    interface for DEIMOS redshift quality control.  
;
; CALLING SEQUENCE:
;    mg_zspec, maskname, zspecfile=zspecfile, /deepfile, /seredip
;
;    Need to set path outside idl named 'DEIMOS_DATA' pointing to the
;    directory above the mask files.
;
; INPUTS:
;     maskname  -- name of mask, eg. 'N147_1'
;     zspecfile -- If this keyword is set, return to previous zspec 
;                     session with name eg.,'zspec.N147_1.fits'
;
;
; OUTPUTS:
;       - plots 1D spectrum + model
;       - postage stamp of 2D spectrum +trace       
;       - dialogue box
;
; COMMENTS:
;      Spec1d fits for stellar velocities between -1000 and 1000km/s.  It
;      does not yet fit for galaxies/QSOs.  These need to be flagged by hand.
;
;
; CALLS
;     mg_redo_extract1d
;     mg_oneslit_zfind
;
; MODIFICATION HISTORY:
;   mg 4/06 based on DEEP zspec.pro, mg 4/06
;   mg 2/07 switch to zspec1d files, mg 2/07
;   mg 7/09 changed zquality system, auto ID alignment stars
;           TODO:  add code to determine ra/dec of serendips
;- 
;--------------------------------------------------------------------

PRO deimos_zspec_event, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult
 
	widget_control, ev.id, get_uvalue=value
	if (n_elements(value) eq 0) then value = ''
	name = strmid(tag_names(ev, /structure_name), 7, 4)
	
	case (name) of
		'BUTT': handle_button, ev
		'TEXT': handle_text, ev
    		'DONE': WIDGET_CONTROL, ev.TOP, /DESTROY
  	ENDCASE
END




PRO handle_button, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult

	WIDGET_CONTROL, ev.ID, GET_UVALUE=uval

	case (uval) of
		'Next': show_object, ev, mgzspec_state.ix+1
		'Back': show_object, ev, mgzspec_state.ix-1
		'sky': toggle_twod, ev
		'save': save_file, ev
		'qm2': toggle_qm2, ev
		'qm1': toggle_qm1, ev
		'q0': toggle_q0, ev
		'q1': toggle_q1, ev
		'q2': toggle_q2, ev
		'q3': toggle_q3, ev
		'q4': toggle_q4, ev
		'extr': oned_extract_auto, ev
		'tweek': show_oned, ev
		'zguess_gal': zguess_gal, ev
		'zguess_star': zguess_star, ev

		'DONE': WIDGET_CONTROL, ev.TOP, /DESTROY
	ENDCASE
    END

	
PRO handle_text, ev

	widget_control, ev.id, get_uvalue=uval

	case (uval) of
		'smooth1d': oned_smooth, ev
		'tweek': oned_tweek, ev
		'extr': oned_extract, ev
		'newstr': new_star, ev
		'comment': comment, ev
		'jump_to_slit': jump_to_slit, ev
	end

END

PRO save_file, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult


      ; CALCULATE SERENDIPS RA/DECs
        print,'Getting final RA/DEC of serendip slits'
        mg_y2radec_serendip,zresult

        ;; Write fits file
        mwrfits,zresult,mgzspec_state.filename,/create
        print,'Saved file to ', mgzspec_state.filename
        show_object, ev, mgzspec_state.ix
END




PRO toggle_twod, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult

	mgzspec_state.show_twod = not mgzspec_state.show_twod
	show_twod

    END



PRO toggle_extract, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult

        r = string(mgzspec_state.r1,format='(i2)')+ ', '+$
          string(mgzspec_state.r2,format='(i2)')
	widget_control, ev.id, set_value = r
		

        read_file2d
	show_twod

END





PRO toggle_oned, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult

	mgzspec_state.show_oned = not mgzspec_state.show_oned
	mgzspec_state.show_ktrs_img = 0

	show_object, ev, mgzspec_state.ix
END



PRO show_oned
	COMMON mgzspec_specs_common, mgzspec_state,zresult
	COMMON npk_specs_common, npk_state
	COMMON splot_state,state,graphkeys
	COMMON mgzspec_path, data_path, mgidl_path,mask

	tellam = [6860,7600]

  	 slit = mgzspec_state.ix
         spec = mgzspec_state.spec1d
         wave = mgzspec_state.wave1d
         ivar  = mgzspec_state.ivar1d


     ; REMOVE TELLURIC LINES
;        q1=where(wave ge 6866 and wave le 6912)
;        q3=where(wave ge 7593 and wave le 7690)
;        q=[q1,q3]
;        ivar[q] =ivar[q]*0.0001   ; THERE ARE THE STRONG A- B-BANDS

;        q5=where(wave ge 8850 and wave le 9120) 
;        if (q5[0] ne -1) then ivar[q5] =ivar[q5]*0.1   ; THERE ARE THE WEAKER ABS LINES
;        var = 1./ivar


         zans=zresult[slit]
         eigendir = mgidl_path + '/templates'

      ; HACK TO SKIP BAD SLITS
        if (NOT Float( Finite(zans.rchi2))) then begin
                splot,wave,spec 
                goto,skipz
             endif
        if (Float( Finite(zans.rchi2))) then begin
            if (zans.tcolumn[0] eq -1) then zans.tcolumn[0] = 0
            synth = synthspec(zans, loglam=alog10(wave),eigendir = eigendir)
        endif

	if mgzspec_state.smooth1d[0] eq 1 then begin;
		spec = poly_smooth(spec, mgzspec_state.smooth1d[1])
		synth= poly_smooth(synth, mgzspec_state.smooth1d[1])
                print,'smoothing = ',mgzspec_state.smooth1d[1] 
        endif

      ; set limits not so affected by bad points
        window,0

        slit = mgzspec_state.ix
        n = n_elements(mgzspec_state.linenames)
        zlines = zresult[slit].z * mgzspec_state.linewvl  + $
                   mgzspec_state.linewvl

        !p.multi=[0,1,2]
        plot,wave,spec,xrange=[8470,8570],xstyle=1
        oplot,wave,synth,color=3
        for i = 0, n-1 do begin
                wl = zlines[i]
                wn = mgzspec_state.linenames[i]
                oplot, [wl, wl], [-500, 500000], color=2,linestyle=1
        endfor
;        plot,wave,spec*sqrt(ivar),xrange=[8470,8570],xstyle=1
        plot,wave,spec,xrange=[6540,6600],xstyle=1
        oplot,wave,synth,color=3
        for i = 0, n-1 do begin
                wl = zlines[i]
                wn = mgzspec_state.linenames[i]
                oplot, [wl, wl], [-500, 500000], color=2,linestyle=1
        endfor
        !p.multi=0



        q=where(spec ge 0 and wave ge 6000)
        medspec = median(spec)
        stdspec = stddev(spec)
        maxs = medspec + 3.*stdspec 
        mins = medspec - 2.*stdspec
        mins_new = mins - 0.1*(maxs-mins)
        if (mins ge 0) then mins=0

	if (size(state))[0] eq 0 then begin
		tely = 1.5
		splot, wave,spec,yrange=[mins,maxs],$
                              xrange=[6400,9200] 
	endif else begin
		if mgzspec_state.persistent eq 0 then $
		     splot, wave,spec,yrange=[mins,maxs], $
                              xrange=[6400,9200] $
		else $
		     splot, wave,spec, $
			xrange=[6400,9200], yrange=state.yrange
	endelse
        
        
        soplot,wave,synth,color=4
;        soplot,wave,sqrt(var),color='red'



	tely=state.yrange[1]
	sxyouts, tellam[0], tely, 'B', color=1,charthick=3
	sxyouts, tellam[1], tely, 'A', color=1,charthick=3
	overplot_lines
        skipz:a=1

END


PRO overplot_lines
	COMMON mgzspec_specs_common, mgzspec_state,zresult

        slit = mgzspec_state.ix
        n = n_elements(mgzspec_state.linenames)
        zlines = zresult[slit].z * mgzspec_state.linewvl  + $
                   mgzspec_state.linewvl
        for i = 0, n-1 do begin
                wl = zlines[i]
                wn = mgzspec_state.linenames[i]

                soplot, [wl, wl], [-500, 500000], color=2,linestyle=1
        endfor

END



PRO show_twod
	COMMON mgzspec_specs_common, mgzspec_state,zresult

      ; SET 2D SPECTRUM Sky/NoSKY, NEXP = 1/2
        twod=mgzspec_state.spec2d
        wave2d = mgzspec_state.wave2d

        slit = mgzspec_state.ix
        s=size(twod)
        mx = 2 * mgzspec_state.ycen
        if (mx ge 99) then mx = 99
        spec = twod[*,0:*]

        xatv,spec,/align

 
      ; wavelength array along trace
        s=size(twod)
        wave=fltarr(s(1))
        bwave=fltarr(s(1))
        for j=0,s(1)-1 do begin
            wave[j] = wave2d[j,mgzspec_state.ycen]
            bwave[j] = wave2d[j,mgzspec_state.b_ycen+mgzspec_state.ycen]
        endfor

      ; plot extraction region
        q=where(wave ge 8500 and wave le 8520,n)
        xatvplot,[0,50],[mgzspec_state.r1,mgzspec_state.r1],color=yellow
        xatvplot,[0,50],[mgzspec_state.r2,mgzspec_state.r2],color=yellow

        bs = mgzspec_state.b_ycen
        xatvplot,[0,50],[bs+mgzspec_state.r1,bs+mgzspec_state.r1],color=yellow
        xatvplot,[0,50],[bs+mgzspec_state.r2,bs+mgzspec_state.r2],color=yellow

        
        if (n ne 0) then begin
          xatvplot,[q[0],q[n-1]],[mgzspec_state.r1,mgzspec_state.r1],color=yellow
          xatvplot,[q[0],q[n-1]],[mgzspec_state.r2,mgzspec_state.r2],color=yellow

          xatvplot,[q[0],q[n-1]],[bs+mgzspec_state.r1,bs+mgzspec_state.r1],color=yellow
          xatvplot,[q[0],q[n-1]],[bs+mgzspec_state.r2,bs+mgzspec_state.r2],color=yellow
        endif
      ; plot location of important lines
        zlines = zresult[slit].z * mgzspec_state.linewvl  + $
                                        mgzspec_state.linewvl

;        vactoair,zlines
        tmp =findgen(10)
        xline= fltarr(10)
        yline1 = tmp
        yline2 = tmp+40
        for j=0,n_elements(mgzspec_state.linenames)-1 do begin

           ; RED HALF
            diff = min(abs(wave - zlines[j]),q)
            if (q ne 0) then xatvxyouts,q,50,mgzspec_state.linenames[j],charsize=2,color='yellow',charthick=4

            xline[0:*] = q
            xatvplot,[xline],[yline1],color='yellow'
            xatvplot,[xline],[yline2],color='yellow'

            xatvplot,[xline],[yline1],color='yellow'
            xatvplot,[xline],[yline2],color='yellow'

           ; BLUE HALF
            diff = min(abs(bwave - zlines[j]),q)
            if (q ne 0) then xatvxyouts,q,bs+50,mgzspec_state.linenames[j],charsize=2,color='yellow',charthick=4
            xline[0:*] = q
            xatvplot,[xline],[bs+yline1],color='yellow'
            xatvplot,[xline],[bs+yline2],color='yellow'

            xatvplot,[xline],[bs+yline1],color='yellow'
            xatvplot,[xline],[bs+yline2],color='yellow'
        endfor

        

END



PRO show_object, ev, ix
	COMMON mgzspec_specs_common, mgzspec_state, zresult
	common splot_state, state, graphkeys

	if ix lt 0 or ix gt n_elements(zresult.slitname)-1 then begin
		null = dialog_message(['There are not more objects!'])
		return
	endif


	mgzspec_state.ix = ix
        read_file2d
        read_file1d
	describe_spectrum, ev

	show_oned
	show_twod


END



PRO show_message, ev, msg

        WIDGET_CONTROL, ev.TOP, GET_UVALUE=textwid
        WIDGET_CONTROL, textwid, SET_VALUE=msg[1]

END

PRO describe_spectrum, ev
	COMMON mgzspec_specs_common, mgzspec_state, zresult
        get_info,msg
	show_message, ev, msg
    END



PRO get_info,msg
	COMMON mgzspec_specs_common, mgzspec_state, zresult
        
        c=2.9979e5

        slit = mgzspec_state.ix
        id = zresult[slit].slitname

      ; Mask info
        mask  = strcompress('Mask ID =  ' + zresult[slit].maskname)
        obj= strcompress('Slit Num    = '  + string(zresult[slit].slitname,$
                                                   format='(a10)'))
        if (zresult[slit].pri eq -1) then $
               obj= strcompress(obj  + ' align')

        obj2 =strcompress('Slit Name    = '  + zresult[slit].objname)
        zstr  = strcompress('V = '  + string(zresult[slit].z*c,format='(f8.1)'))
        astr  = strcompress('Aband = '  + string(zresult[slit].aband*c,format='(f8.1)'))
        zqual = strcompress('ZQUAL   = '  + string(zresult[slit].zquality,format='(i4)'))
       
        subc  = strcompress(zresult[slit].subclass)



        msg1 = obj2 +string([10B])+$
               obj +string([10B]) +$
               string([10B])+$
               zstr+string([10B]) +$
               astr+string([10B]) +$
               subc+string([10B]) +$
               zqual

        msg = [mask,msg1]

      ; PRINT TEMPLATE INFORMATION OUT
        objstr= 'slit id    = '  + string(zresult[slit].objname,format='(a10)')
        print,string(mgzspec_state.ix)+' '+objstr
        print

        q=where(zresult[slit].tmpl_rchi2 ne 0,n)
        print,zresult[slit].subclass,zresult[slit].z*3e5
        print
        for jj=0,n-1 do begin
           print,jj,string(zresult[slit].tmpl_rchi2[jj],f='(f7.2)'),$
                 '   ',zresult[slit].tmpl_tfile[jj],$
                 zresult[slit].tmpl_z[jj]*2.9979e5

        endfor



END


 
  ; TOGGLE SPECTRA QUALITY BUTTONS
    PRO toggle_qm2, ev
	COMMON mgzspec_specs_common, mgzspec_state, zresult
	zresult[mgzspec_state.ix].zquality = -2
	describe_spectrum, ev
     END
    PRO toggle_qm1, ev
	COMMON mgzspec_specs_common, mgzspec_state, zresult
	zresult[mgzspec_state.ix].zquality = -1
	describe_spectrum, ev
     END
    PRO toggle_q0, ev
	COMMON mgzspec_specs_common, mgzspec_state, zresult
	zresult[mgzspec_state.ix].zquality = 0
	describe_spectrum, ev
     END
    PRO toggle_q1, ev
	COMMON mgzspec_specs_common, mgzspec_state, zresult
	zresult[mgzspec_state.ix].zquality = 1
	describe_spectrum, ev
    END
    PRO toggle_q2, ev
	COMMON mgzspec_specs_common, mgzspec_state, zresult
	zresult[mgzspec_state.ix].zquality = 2
	describe_spectrum, ev
    END
    PRO toggle_q3, ev
	COMMON mgzspec_specs_common, mgzspec_state, zresult
	zresult[mgzspec_state.ix].zquality = 3
	describe_spectrum, ev
    END
    PRO toggle_q4, ev
	COMMON mgzspec_specs_common, mgzspec_state, zresult
	zresult[mgzspec_state.ix].zquality = 4
	describe_spectrum, ev
    END


PRO oned_smooth, ev
	COMMON mgzspec_specs_common, mgzspec_state

	widget_control, ev.id, get_value=val
	mgzspec_state.smooth1d[1] = fix(val)
	widget_control, ev.id, set_value = $
		strcompress(string(mgzspec_state.smooth1d[1]),/remove_all)
        show_oned

    END

PRO oned_tweek, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult

	widget_control, ev.id, get_value=val
	zresult[mgzspec_state.ix].z = float(val)/2.9979e5
	widget_control, ev.id, set_value = $
		strcompress(string(zresult[mgzspec_state.ix].z*2.9979e5),/remove_all)
        show_oned
        show_twod

    END


; AUTOMATICALLY RE-DEFINE AND EXTRACT 1D SPECTRA
PRO oned_extract_auto, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult
	COMMON mgzspec_path, data_path, mgidl_path,mask

	widget_control, ev.id, get_value=val
        

      ; RE-EXTRACT 
        if (NOT mgzspec_state.deepfile) then begin
          file1d = strcompress(zresult[mgzspec_state.ix].spec1d_file,/remove_all) 
          file2dr = strcompress(zresult[mgzspec_state.ix].spec2d_file,/remove_all) 
          bpos=strpos(file2dr,'R.fits')
          file2db = file2dr
           strput,file2db,'B.fits',bpos
        endif

        if (mgzspec_state.deepfile) then begin
               s = zresult[mgzspec_state.ix].slitname 
               slit = string(zresult[mgzspec_state.ix].slitname) 
               if (s le 99) then slit = '0' + zresult[mgzspec_state.ix].slitname 
               if (s le 9) then slit = '00' + zresult[mgzspec_state.ix].slitname 
               
               rj = getenv('RAJAPATH') 
               date = ''
               if (rj ne '') then $
                    date = '20[0-9][0-9][a-z][a-z][a-z][0-9][0-9]'

               file1d = strcompress( '/'+mask + '/'+date+'/spec1d.' + $
                                     mask + '.' + slit + '.'+$
                                     zresult[mgzspec_state.ix].objname+$
                                     '.fits.gz',$
                                     /remove_all) 

                file2dr = strcompress('/'+mask + '/'+date+'/slit.' + $
                                     mask + '.' + slit + $
                                     'R.fits.gz',$
                                     /remove_all) 
                file2db = strcompress( '/'+mask + '/'+date+'/slit.' + $
                                     mask + '.' + slit + $
                                     'B.fits.gz',$
                                     /remove_all) 
         endif

        print,file1d
        oldfile1d=file1d
        mg_redo_extract1d, oldfile1d, [file2dr,file2db], objpos, /auto
        oldfile1d=file1d
        print,oldfile1d
        mg_oneslit_zfind, zresult,oldfile1d,mgzspec_state.ix

        read_file1d
        read_file2d
        show_oned
        show_twod
        get_info,msg

     END

PRO oned_extract, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult
	COMMON mgzspec_path, data_path, mgidl_path,mask

	widget_control, ev.id, get_value=val
        tmp = strsplit(val,',',/extract)

        mgzspec_state.r1  = float(tmp[0])
        mgzspec_state.r2  = float(tmp[1])
        objpos=[mgzspec_state.r1,mgzspec_state.r2]
       

      ; RE-EXTRACT
      ; RE-EXTRACT 
        if (NOT mgzspec_state.deepfile) then begin
          file1d = strcompress(zresult[mgzspec_state.ix].spec1d_file,/remove_all) 
          file2dr = strcompress(zresult[mgzspec_state.ix].spec2d_file,/remove_all) 
          bpos=strpos(file2dr,'R.fits')
          file2db = file2dr
           strput,file2db,'B.fits',bpos
        endif

        if (mgzspec_state.deepfile) then begin
               s = zresult[mgzspec_state.ix].slitname 
               slit = string(zresult[mgzspec_state.ix].slitname) 
               if (s le 99) then slit = '0' + zresult[mgzspec_state.ix].slitname 
               if (s le 9) then slit = '00' + zresult[mgzspec_state.ix].slitname 
               
               rj = getenv('RAJAPATH') 
               date = ''
               if (rj ne '') then $
                    date = '20[0-9][0-9][a-z][a-z][a-z][0-9][0-9]'


              file1d = strcompress('/' + mask + '/'+date + '/spec1d.' + $
                                     mask + '.' + slit + '.'+$
                                     zresult[mgzspec_state.ix].objname+$
                                     '.fits',$
                                     /remove_all) 

                file2dr = strcompress( '/' + mask + '/'+date+'/slit.' + $
                                     mask + '.' + slit + $
                                     'R.fits.gz',$
                                     /remove_all) 
                file2db = strcompress( '/' + mask + '/'+date+'/slit.' + $
                                     mask + '.' + slit + $
                                     'B.fits.gz',$
                                     /remove_all) 
           endif 

        oldfile1d=file1d
        mg_redo_extract1d, oldfile1d, [file2dr,file2db], objpos
        oldfile1d=file1d
        print,oldfile1d
        mg_oneslit_zfind, zresult,oldfile1d,mgzspec_state.ix

        read_file1d
        read_file2d
        show_oned
        show_twod
        get_info,msg
    END



; AUTOMATICALLY RE-DEFINE AND EXTRACT 1D SPECTRA
PRO new_star, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult
	COMMON mgzspec_path, data_path, mgidl_path,mask
        COMMON templates, tmp_subclass

	widget_control, ev.id, get_value=val
        
        slit = mgzspec_state.ix

        zresult[slit].subclass = zresult[slit].tmpl_tfile[val]
        qst = where(strcompress(zresult[slit].subclass,/rem) eq tmp_subclass)
        if (strcompress(zresult[slit].subclass,/rem) eq 'galaxy') then qst=16
        istar = qst        

        zresult[slit].tcolumn[0] = istar
        zresult[slit].z = zresult[slit].tmpl_z[val]

      ; TOGGLE BETWEEN STAR/GALAXY TEMPLATES
       if (strcompress(zresult[slit].tmpl_tfile[val],/rem) eq 'galaxy') then begin
              zresult[slit].tfile = 'deimos-galaxy.fits'
              zresult[slit].npoly=2
       endif

        if (strcompress(zresult[slit].tmpl_tfile[val],/rem) ne 'galaxy' and $
           strcompress(zresult[slit].tfile,/rem) eq 'deimos-galaxy.fits') then begin

           zresult[slit].tfile = zresult[slit+1].tfile   ; hack
           zresult[slit].npoly=4
        endif

        print,zresult[slit].subclass
        print,zresult[slit].z*3e5

        read_file1d
        show_oned
        describe_spectrum, ev
     END

; AUTOMATICALLY RE-DEFINE AND EXTRACT 1D SPECTRA
PRO comment, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult
	COMMON mgzspec_path, data_path, mgidl_path,mask

	widget_control, ev.id, get_value=val
        
        slit = mgzspec_state.ix
        zresult[slit].comment = val


        get_info,msg
    END




PRO jump_to_slit, ev
	COMMON mgzspec_specs_common, mgzspec_state
	COMMON mgzspec_path, data_path, mgidl_path,mask

	widget_control, ev.id, get_value=val
	mgzspec_state.ix = fix(val)
	widget_control, ev.id, set_value = $
		strcompress(string(mgzspec_state.ix),/remove_all)
        show_object,ev,mgzspec_state.ix

    END

 ; RERUN ZSPEC ON CASES WERE ZFIND FAILS
;    FIT GALAXY TEMPLATE ONLY IN NARROWER RANGE THAN ZFIND
  PRO zguess_gal, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult

	widget_control, ev.id, get_value=val
           if (NOT mgzspec_state.deepfile) then begin
                file1d = strcompress($
                                     zresult[mgzspec_state.ix].spec1d_file,$
                                     /remove_all)    
           endif
           if (mgzspec_state.deepfile) then begin
               s = zresult[mgzspec_state.ix].slitname 
               slit = string(zresult[mgzspec_state.ix].slitname) 
               if (s le 99) then slit = '0' + zresult[mgzspec_state.ix].slitname 
               if (s le 9) then slit = '00' + zresult[mgzspec_state.ix].slitname 

               rj = getenv('RAJAPATH') 
               date = ''
               if (rj ne '') then $
                    date = '20[0-9][0-9][a-z][a-z][a-z][0-9][0-9]'
               
                file1d = strcompress('/' + mask + '/'+date+'/spec1d.' + $
                                     mask + '.' + slit + '.'+$
                                     zresult[mgzspec_state.ix].objname+$
                                     '.fits',$
                                     /remove_all) 
           endif                   

          mg_oneslit_zfind, zresult,file1d,mgzspec_state.ix,/galaxy
          show_oned
          show_twod
          get_info,msg


    END

 ; RERUN ZSPEC ON CASES WERE ZFIND FAILS
;    FIT STAR ONLY
  PRO zguess_star, ev
	COMMON mgzspec_specs_common, mgzspec_state,zresult
	COMMON mgzspec_path, data_path, mgidl_path,mask

;	widget_control, ev.id, get_value=val
           if (NOT mgzspec_state.deepfile) then begin
                file1d = strcompress($
                                     zresult[mgzspec_state.ix].spec1d_file,$
                                     /remove_all)    
           endif
           if (mgzspec_state.deepfile) then begin
               s = zresult[mgzspec_state.ix].slitname 
               slit = string(zresult[mgzspec_state.ix].slitname) 
               if (s le 99) then slit = '0' + zresult[mgzspec_state.ix].slitname 
               if (s le 9) then slit = '00' + zresult[mgzspec_state.ix].slitname 
               
               rj = getenv('RAJAPATH') 
               date = ''
               if (rj ne '') then $
                    date = '20[0-9][0-9][a-z][a-z][a-z][0-9][0-9]'


                file1d = strcompress('/' + mask + '/'+date + '/spec1d.' + $
                                     mask + '.' + slit + '.'+$
                                     zresult[mgzspec_state.ix].objname+$
                                     '.fits',$
                                     /remove_all) 
           endif           

          q=where(zresult.sn ge 2.5 and zresult.z*3e5 le -100. and $
                 zresult.z*3e5 ge -900.)
          zmed = median(zresult[q].z * 2.997e5)
          print,'search +/- 100 km/s around v=',zmed
          mg_oneslit_zfind, zresult,file1d,mgzspec_state.ix,zlimits=zmed + [-200,200]
          show_oned
          show_twod
          get_info,msg
    END


PRO read_file1d
	COMMON mgzspec_specs_common, mgzspec_state,zresult
	COMMON mgzspec_path, data_path, mgidl_path,mask

           if (NOT mgzspec_state.deepfile) then begin
                file1d = strcompress(data_path  +$
                                     zresult[mgzspec_state.ix].zspec1d_file,$
                                     /remove_all)    
           endif
           if (mgzspec_state.deepfile) then begin
               s = zresult[mgzspec_state.ix].slitname 
               slit = string(zresult[mgzspec_state.ix].slitname) 
               if (s le 99) then slit = '0' + zresult[mgzspec_state.ix].slitname 
               if (s le 9) then slit = '00' + zresult[mgzspec_state.ix].slitname 
               

               rj = getenv('RAJAPATH') 
               date = ''
               if (rj ne '') then $
                 date = '20[0-9][0-9][a-z][a-z][a-z][0-9][0-9]'

                file1d = strcompress(data_path + '/' + mask + '/'+date+'/zspec1d.' + $
                                     mask + '.' + slit + '.'+$
                                     zresult[mgzspec_state.ix].objname+$
                                     '.fits.gz',$
                                     /remove_all) 
           endif           

                print, 'reading ',file1d
                data = mrdfits(file1d,1,status=status)

                mgzspec_state.spec1d[0:*] = 0
                if (status ne -1) then begin
                   s=size(data.spec)
                   nx= s(1)-1
                   mgzspec_state.spec1d[0:nx] = data.spec
                   mgzspec_state.ivar1d[0:nx] = data.ivar
                   mgzspec_state.wave1d[0:nx] = data.lambda

                   mgzspec_state.ycen  = avg(data.objpos)
                   mgzspec_state.r1  = data.r1
                   mgzspec_state.r2  = data.r2
                endif


                mgzspec_state.read_1d = 1

        END

PRO read_file2d
	COMMON mgzspec_specs_common, mgzspec_state,zresult
	COMMON mgzspec_path, data_path, mgidl_path, mask


           if (NOT mgzspec_state.deepfile) then begin
                file2d = strcompress(data_path + $
                                     zresult[mgzspec_state.ix].spec2d_file,$
                                     /remove_all) 
           endif
           if (mgzspec_state.deepfile) then begin
               s = zresult[mgzspec_state.ix].slitname 
               slit = string(zresult[mgzspec_state.ix].slitname) 
               if (s le 99) then slit = '0' + zresult[mgzspec_state.ix].slitname 
               if (s le 9) then slit = '00' + zresult[mgzspec_state.ix].slitname 
               

               rj = getenv('RAJAPATH') 
               date = ''
               if (rj ne '') then $
                    date = '20[0-9][0-9][a-z][a-z][a-z][0-9][0-9]'

                file2d = strcompress(data_path + '/' + mask + '/'+date+'/slit.' + $
                                     mask + '.' + slit + $
                                     'R.fits.gz',$
                                     /remove_all) 
           endif

                print,'Reading ',file2d               
                nf =where(strcompress(findfile(file2d),/rem) eq '')
                

                if (nf eq -1) then begin
                 ; READ RED SIDE
                   r_spec2d = mrdfits(file2d,1,status=status) 

                   mgzspec_state.spec2d = 0
                   if (status ne -1) then begin
                      mgzspec_state.spec2d = r_spec2d.flux
                      mgzspec_state.wave2d = lambda_eval(r_spec2d)
                   endif

                 ; READ BLUE SIDE
                   b_file2d = repstr(file2d,'R.fits.gz','B.fits.gz')
                   b_spec2d = mrdfits(b_file2d,1,status=status)                


                   if (status ne -1) then begin
                      s=size(b_spec2d.flux)
                      s1 = s(2) +50
                      print,s1,s(2)

                      mgzspec_state.b_ycen = s1    
                      mgzspec_state.spec2d[*,s1:s1+s(2)-1] = b_spec2d.flux
                      bw = lambda_eval(b_spec2d)
                      mgzspec_state.wave2d[*,s1:s1+s(2)-1] = bw
                   endif

                endif

                if (nf eq 0) then begin
                   mgzspec_state.spec2d = 0
                   mgzspec_state.wave2d = 0
                endif


                mgzspec_state.read_2d = 1
                

END

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------

Pro sg_zspec,maskname,zspecfile = zspecfile, SERENDIP=serendip, $
             deepfile=deepfile,$
             write_only=write_only



	COMMON mgzspec_specs_common, mgzspec_state,zresult
	COMMON mgzspec_path, data_path,mgidl_path,mask
        COMMON templates, tmp_subclass

        !p.font =6


     ; SET THIS PATH TO DEIMOS DATA WORKING DIRECTORY
       mgidl_path = getenv('SGIDL')
       data_path = getenv('DEIMOS_DATA')
       zresultfile = data_path + '/zresult/mgzresult.'+maskname+'.fits'
       output = data_path+'zresult/zspec.'+maskname+'.fits'

       if (keyword_set(zspecfile)) then begin
           zresultfile = data_path+'zresult/'+zspecfile
           print,'WARNING!  overwriting results to file ',zspecfile
           output=zresultfile
           ngd=n_elements(zresult)
      endif


      ; READ ZRESULT OUTPUT
        print,'Reading ',zresultfile
        x=findfile(zresultfile,count=cnt)
        if (cnt eq 0) then begin
             print,zresultfile
             print,'File not found!!'
            return
        endif
        zresult = mrdfits(zresultfile,1)

      ; DETERMINE PRIMARY SLITS  
        serendip=0
        if NOT (keyword_set(zspecfile)) then begin
           objname = strcompress(zresult.objname,/rem)
           gslit = where(objname ne 'serendip1' and $
                      objname ne 'serendip2' and $
                      objname ne 'serendip3' and $
                      objname ne 'serendip4',ngd)
           read,serendip,prompt = 'INCLUDE SERENDIP SLITS? [0=no, 1=yes]: '

         ; DON"T ALLOW MORE THAN 2 SERENDIPS PER SLIT!!
           if (serendip eq 1) then begin
              objname = strcompress(zresult.objname,/rem)
              gslit = where(objname ne 'serendip6' and $
                           objname ne 'serendip7' and $ 
                           objname ne 'serendip8',ngd)   

           endif
           zresult=zresult[gslit]
        endif 


      ; Important lines, convert these to vacuum wvls
        lines = [3726.16,3728.9,$   
                 3933.7,3968.5,$
                 4861.32,4958.9, 5006.84,$
                 6548.1,6562.80, 6583.6,$
                 8183.25, 8194.79,$
                 8498.0,8542.1,8662.2]
        ;airtovac,lines
        nlines= ['[OII]',' ',$
                 ' ','CaH/K',$
                 'H!7b!6',' ',' ',$
                 ' ','H!7a!6',' ','NaI','NaI', 'CaII', 'b', 'c']

	mgzspec_state = {$
		ix: 0, $		; the object # being examined
                nslit:0.0,$             ; number of good slits
                frame: ' ',$            ; frame          
                smooth1d: [1,5], $	; [smooth y/n, smoothing in pixels]
		persistent: 0,$		; Are splot options persistent
		show_twod: 1, $		; Show the twod ATV plot
		show_nexp: 0, $		; Show the nexp 
		read_1d: 0, $		; has nexp1 been read
		read_2d: 0, $		; has nexp2 been read
		linewvl: lines,$	;  list of names to display lines for
		linenames: nlines,$	;  list of names to display lines for
                zquality:0,$               ; quality flag
		filename: output, $	;  
                ycen:0.0,$          
                b_ycen:0.0,$          
                spec2d:fltarr(4096,2500),$
                wave2d:fltarr(4096,2500),$
                spec1d:fltarr(7200),$
                wave1d:fltarr(7200),$
                ivar1d:fltarr(7200),$
                r1:0.0,$
                r2:0.0,$
                serendip:0.0,$
                deepfile:0.0$
	}


        if NOT keyword_set(zspecfile) and (serendip eq 1) then mgzspec_state.serendip=1


        mask=maskname
        mgzspec_state.nslit = ngd
        mgzspec_state.frame = zresultfile
        mgzspec_state.zquality = zresult[mgzspec_state.ix].zquality
        if NOT keyword_set(deepfile) then deepfile=0
        mgzspec_state.deepfile = deepfile

       ; SET ZQUALITY OF ALIGNMENT STARS
       if (NOT keyword_set(zspecfile)) then begin
          qalign = where(zresult.pri eq -1)
          zresult[qalign].zquality = -1                    
       endif

       ; IF CONVERTING zresult TO zspec ONLY
         if (keyword_set(write_only)) then begin
            mwrfits,zresult,mgzspec_state.filename,/create
           return
         endif


      ; READ AND LOAD TEMPLATES
        tfile = strcompress(zresult[0].tfile,/rem)
        if (tfile eq '') then tfile = strcompress(zresult[1].tfile,/rem)
        if (tfile eq '') then tfile = strcompress(zresult[2].tfile,/rem)


         eigenfile = strcompress(mgidl_path + '/templates/'+tfile,/rem)
         shdr = headfits(eigenfile)
         nstar = sxpar(shdr, 'NAXIS2') > 1
         tmp_subclass = strarr(nstar) ;types of stars
         for istar=0, nstar-1 do tmp_subclass[istar] = $
               strtrim( sxpar(shdr, 'NAME'+strtrim(string(istar),2)), 2)


      ; READ 2D SPECTRUM FRAME
        read_file2d
        read_file1d


        get_info,msg
        

	; HIGHEST LEVEL
  	base = WIDGET_BASE(/column)
	WIDGET_CONTROL, /MANAGED, base

	; Level 1

	b1 = WIDGET_BASE(base, /frame, /row)
  	button = WIDGET_BUTTON(b1, VALUE='Back', UVALUE='Back')
  	button = WIDGET_BUTTON(b1, VALUE='Next', UVALUE='Next')
  	text = WIDGET_LABEL(b1, value='jump to slit =')
  	text = WIDGET_TEXT(b1, XSIZE=5,value=string(mgzspec_state.ix),$
                           /EDITABLE,uvalue='jump_to_slit')

	b2 = WIDGET_BASE(base, /frame, /row)



	; Level 2 Middle
	l2_right = WIDGET_BASE(b2, /column, /frame)
	l2_right_top = WIDGET_BASE(l2_right,/row)
        ;button = WIDGET_LABEL(l2_right_top, VALUE='Comment')
        ;text2 = WIDGET_TEXT(l2_right_top, XSIZE=15, $
        ;        value=zresult[mgzspec_state.ix].comment, $
	;	/EDITABLE, uvalue='comment')        
        blah = WIDGET_LABEL(l2_right_top, VALUE='Smoothing:')
        button = WIDGET_TEXT(l2_right_top, XSIZE=5, value='5', /editable, uvalue='smooth1d')

	l2_right_mid = WIDGET_BASE(l2_right, /row)
	button = WIDGET_LABEL(l2_right_mid, value='Tweek vel')
	text2 = WIDGET_TEXT(l2_right_mid, XSIZE=5, value='-240', $
		/EDITABLE, uvalue='tweek')


	l2_right_bot = WIDGET_BASE(l2_right, /row)
	button = WIDGET_BUTTON(l2_right_bot, value='Auto Re-Extract', uvalue='extr')
        r = string(mgzspec_state.r1,format='(i2)')+ ', '+$
          string(mgzspec_state.r2,format='(i3)')

	text2 = WIDGET_TEXT(l2_right_bot, XSIZE=10, value=r, $
		/EDITABLE, uvalue='extr')




	l2_right_bot = WIDGET_BASE(l2_right, /row)
	button = WIDGET_BUTTON(l2_right_bot, value='Save Fits File', uvalue='save')


	b2_r = WIDGET_BASE(b2, /frame, /column)
	b2_rb = WIDGET_BASE(b2, /frame,/column)
	b2_l = WIDGET_BASE(b2, /frame, /column)


      
        ; LEVEL 3

	b1 = WIDGET_BASE(base, /row)
        text = WIDGET_LABEL(b1, value='Quality =')
  	button = WIDGET_BUTTON(b1, VALUE=' -2 ', UVALUE='qm2')
  	button = WIDGET_BUTTON(b1, VALUE=' -1 ', UVALUE='qm1')
  	button = WIDGET_BUTTON(b1, VALUE=' 0 ', UVALUE='q0')
  	button = WIDGET_BUTTON(b1, VALUE=' 1 ', UVALUE='q1')
  	button = WIDGET_BUTTON(b1, VALUE=' 2 ', UVALUE='q2')
  	button = WIDGET_BUTTON(b1, VALUE=' 3 ', UVALUE='q3')
  	button = WIDGET_BUTTON(b1, VALUE=' 4 ', UVALUE='q4')


	l2_right_bot = WIDGET_BASE(b1, /row)
        text = WIDGET_LABEL(b1, value='New Template: ')
	text2 = WIDGET_TEXT(b1, XSIZE=10, value='0' , $
		/EDITABLE, uvalue='newstr')



        text = WIDGET_LABEL(b2_r, value='ZQUALITY FLAGS: ')
        text = WIDGET_LABEL(b2_r, value='-2 = Default/Failed slit',/align_left)
        text = WIDGET_LABEL(b2_r, value='-1 = Alignment Star',/align_left)
        text = WIDGET_LABEL(b2_r, value=' 0 = Galaxy/QSO',/align_left)
        text = WIDGET_LABEL(b2_r, value=' 1 = Velocity by hand',/align_left)
        text = WIDGET_LABEL(b2_r, value=' 2 = Failed star vel.',/align_left)
        text = WIDGET_LABEL(b2_r, value=' 3 = Marginal star vel.',/align_left)
        text = WIDGET_LABEL(b2_r, value=' 4 = Secure star vel.',/align_left)

        text = WIDGET_LABEL(b2_rb, value=msg[0])
        text = WIDGET_LABEL(b2_rb, value=msg[1],/align_left)


  	button4 = WIDGET_BUTTON(base, VALUE='Done', UVALUE='DONE')

        show_oned
        show_twod
        
  	WIDGET_CONTROL, base, SET_UVALUE=text
  	WIDGET_CONTROL, base, /REALIZE
  	XMANAGER, 'deimos_zspec', base



    end






