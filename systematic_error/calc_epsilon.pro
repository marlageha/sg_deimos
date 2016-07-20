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
Pro calc_epsilon,create_ps=create_ps

  c=2.99792e5
  loadct,12
  !p.multi=0


    if keyword_set(create_ps) then begin
         set_plot,'ps'
         device,xoffset=-0.5,yoffset=2,filename= 'error_fig.eps',/times, $
         xsize=23,ysize=12,/portrait,/encapsulate,/color
    endif


    ; READ REPEATS FILE
     file = 'repeats_oldtemplates.fits'
     x=mrdfits(file,1)

   ; PLOT STRAIGHT FUNCTION  
     norm_err = x.vdiff/sqrt(x.rand_error1^2 + x.rand_error2^2)

   ; TRIM OUT OUTLIERS  
     q = where(abs(norm_err) le 10)
     x=x[q]
     norm_err = x.vdiff/sqrt(x.rand_error1^2 + x.rand_error2^2)
     
     plothist,norm_err,a,b,bin=1.0,xrange=[-50,50],charsize=2,xstyle=1
     stop

     n=100.
     epsilon = findgen(n)/10.
     sigma = fltarr(n)
     for i=0,n-1 do begin
     
        norm_err = x.vdiff/sqrt(x.rand_error1^2 + x.rand_error2^2 + 2.*epsilon[i]^2)
        q = where(abs(norm_err) le 10)
     
        plothist,norm_err[q],a,b,bin=0.5
        gfit = gaussfit(a,b,coeff, NTERMS=3)
        print,coeff
        sigma[i] = coeff[2]


    endfor
    tmp=min(abs(sigma-1.0),q)
    print,sigma[q],epsilon[q]
     stop

for i=0,n-1 do begin
     
        norm_err = x.vdiff/sqrt(x.rand_error1^2 + x.rand_error2^2 + 2.*epsilon[i]^2)
        q = where(abs(norm_err) le 10)

        err=x.vdiff
        err[0:*] = 0.001
        mg_maxlike_sigma,norm_err,err,a2,b2
        

        sigma[i] = b2

    endfor
    tmp=min(abs(sigma-1.0),q)
    print,sigma[q],epsilon[q]


 if keyword_set(create_ps) then begin
   device,/close
   set_plot,'x'
 endif

     
end
