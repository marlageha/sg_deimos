;--------------------------------------------------------------------
;+
; NAME:
;
; 
; PURPOSE:
;
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; MODIFICATION HISTORY:
;~~~mg
;- 
;--------------------------------------------------------------------
Pro mg_mccheck, maskname,noplot=noplot

  c=2.99792e5
  loadct,12

   ; SET THIS PATH TO DEIMOS DATA WORKING DIRECTORY
     loadct,12
     data_path = getenv('DEIMOS_DATA')

   ; GET DIRECTORIES AND NAME OF EIGEN TEMPLATES
     sgidl_path = getenv('SGIDL')
     eigendir = sgidl_path + '/idl/mgidl/templates/'

   ; READ IN ALLDATA FILE
     file = data_path+'/zresult/zspec.'+maskname+'.fits'
     print,'Reading '+file
     data = mrdfits(file,1)
     nslit = n_elements(data)

     mctmp = {objname:' ',$
              vwidth:0.0,$
              mcerr_vcorr:0.0,$
              mc_v0:  0.0,$     
              mc_vmed:0.0}
              mc = replicate(mctmp, nslit)



  ; READ MC FILE
     mfile = data_path+'/zresult/mcall_'+maskname+'.fits.gz'
     print,'Reading '+mfile
     mdata = mrdfits(mfile,1) 

;     omfile = data_path+'/zresult/old/mc_'+maskname+'.fits'
;     olddata = mrdfits(omfile,1) 

!p.multi=[0,1,2]
for i=0,nslit-1 do begin


    mc[i].objname = mdata[i].objname
    mc[i].mc_v0 =  data[i].z*c

    
    if (data[i].zquality le 2 or data[i].z*3e5 ge -100) then goto,skipslit

    print,i
    print,'zerr mcerr = ',data[i].z_err*c


    ; HISTOGRAM INITIAL DISTRIBUTION AND SET BINNING
      b=0.5
      plothist,mdata[i].mc_vcorr,x,y,bin=b

      if (n_elements(x) le 5) then begin
        b=0.1
        plothist,mdata[i].mc_vcorr,x,y,bin=b
     endif
      if (n_elements(x) le 5) then begin
        b=0.025
        plothist,mdata[i].mc_vcorr,x,y,bin=b
     endif
      if (n_elements(x) le 5) then begin
        b=0.001
        plothist,mdata[i].mc_vcorr,x,y,bin=b
      endif

    ; FIRST ITERATATION TO GET 'TRUE' DISPERSION
;      if (data[i].z*3e5 le 0) then data[i].z =0
      vzero = mdata[i].mc_vcorr - (data[i].z-data[i].aband)*c

      q=where(abs(vzero) le 20,n)
      if (n le 2) then q=where(abs(vzero) le 200,n)
      if (n le 2) then q=where(abs(vzero) le 2000,n)
      
      plothist,mdata[i].mc_vcorr[q],x,y,bin=b,/overplot
      if (n_elements(x) le 4) then begin
        b=0.05
        plothist,mdata[i].mc_vcorr[q],x,y,bin=b
      endif


    ; SECOND ITERATION -- 5-sigma rejection
      g = gauss_fit(x,y,a)

      oplot,x,g,color=30
      sig = abs(a(2))
      if (sig le 0.5) then sig = 10
      q = where(abs(mdata[i].mc_vcorr - a(1)) le sig*5,nm)
      m=a(1)

;   if NOT keyword_set(noplot) then begin
      plothist,mdata[i].mc_vcorr[q],/overplot,bin=b,color=100
      oplot,[mdata[i].mc_v0,mdata[i].mc_v0],[0,1000],linestyle=0
      oplot,[m,m],[0,1000],linestyle=0,color=30




   mc[i].mcerr_vcorr=100
   if (q[0] ne -1) then begin
    mcerr_vcorr= sqrt( total((mdata[i].mc_vcorr[q] - m)^2)/float(nm))

     mc[i].mcerr_vcorr = sqrt(mcerr_vcorr^2 + 2.2^2)  ; add base error

     mc[i].mc_vmed = m

    print,'****'
    print,'z0 = ',(data[i].z - data[i].aband)*c
    print,'V0 = ',mdata[i].mc_v0
    print,'Fitted v = ',m
    print,'Fitted error = ',mc[i].mcerr_vcorr
    print,'****'

  ; BUILT IN STOPS
    redo = 0
    print,sig
    if (abs(m) ge 1e3) then print,'**** REDO=1'
    if (abs(m) ge 1e3) then redo=1
    r=nm*1./n_elements(mdata[i].mc_vcorr)
    if (r le 0.7) then print,'**** REDO=1'
    if (r le 0.7) then redo=1
    if (r le 0.7) then sig=10
    stop
    while (redo eq 1) do begin
       q = where(abs(mdata[i].mc_vcorr - m) le sig*5,nm)
       plothist,mdata[i].mc_vcorr[q],/overplot,bin=b,color=180

       mcerr_vcorr= sqrt( total((mdata[i].mc_vcorr[q] - m)^2)/float(nm))
       mc[i].mcerr_vcorr = sqrt(mcerr_vcorr^2 + 2.2^2) ; add base error
       mc[i].mc_vmed = m

        print,'****'
        print,'z0 = ',(data[i].z- data[i].aband)*2.9979e5
        print,'V0 = ',mdata[i].mc_v0
        print,'Fitted v = ',m
        print,'Fitted error = ',mc[i].mcerr_vcorr
        print,'****'
        print,sig
        stop
     endwhile  
    endif
     mc[i].vwidth = sig



    skipslit:a=0
 endfor

    mwrfits, mc, data_path+'/zresult/mc_'+maskname+'.fits', /create

    !p.multi=[0,1,2]
    plotsym,0,1,/fill
    q=where(data.zquality ge 3)
    plot,data[q].sn,mc[q].mcerr_vcorr,psym=8,xrange=[0,80],yrange=[0,40]
    
    x=findgen(1000)/10.
    oplot,x,3.0+28.*exp(-0.3*x)


stop

end
