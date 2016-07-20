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
;
;- 
;--------------------------------------------------------------------
Pro sg_spec1dstrct__define

 
tmp = {sg_spec1dstrct, $
         class:'', $
         subclass:'', $
         objname: ' ', $
         slitname: ' ',$
         maskname: ' ',$
         date:' ',     $
         mjd:0.0D,     $
         z:0.0, $
         z_err:0.0, $
         rchi2:0.0, $
         dof:0L, $
         rchi2diff:0.0, $
         tfile:'', $
         tcolumn:lonarr(10) - 1L, $
         tmpl_rchi2:fltarr(35),$
         tmpl_z:fltarr(35),$
         tmpl_tfile:strarr(35),$
         npoly:0L, $
         theta:fltarr(10), $
         vdisp:0.0, $
         vdisp_err:0.0, $
         zquality:0, $
         eigendir:' ', $
         ra:' ',$
         dec:' ',$
         pri:0,$
         imag: 0., $   
         rmag: 0., $   
         spec1d_file:' ',$
         zspec1d_file:' ',$
         spec2d_file:' ',$
         comment:' ',$
         airmass:0.0,$
         SN:0.0,$
         aband:0.0,$
         vhelio:0.0,$
         vcorr:0.0$
     }


 
end
   
