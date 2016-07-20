; SEARCH THROUGH ALL ZSPEC AND MCERR DATA
; FIND REPEATED MEASUREMENTS AND SAVE THEM 
; TO A FITS TABLE
;
Pro add_to_repeats,r,s,i,i1,i2
  r[i].vdiff = (s[i1].z - s[i2].z)*3e5

  r[i].subclass1  = s[i1].subclass
  r[i].subclass2  = s[i2].subclass

  r[i].objname1  = s[i1].objname
  r[i].objname2  = s[i2].objname

  r[i].z1  = s[i1].z
  r[i].z2  = s[i2].z

  r[i].ra1  = s[i1].ra
  r[i].ra2  = s[i2].ra

  r[i].dec1  = s[i1].dec
  r[i].dec2  = s[i2].dec

  r[i].sn1  = s[i1].sn
  r[i].sn2  = s[i2].sn

  r[i].mjd1  = s[i1].mjd
  r[i].mjd2  = s[i2].mjd
  
  r[i].airmass1  = s[i1].airmass
  r[i].airmass2  = s[i2].airmass

  r[i].rand_error1  = s[i1].rand_error
  r[i].rand_error2  = s[i2].rand_error

end
    
      



Pro find_repeats

  path = '~/Dropbox/Projects/zresult/'
  zfiles = findfile(path+'zspec*')
  nfiles = n_elements(zfiles)

 ; CREATE TMP STRUCTURE FOR ALL DATA
  temp = {stmp,  subclass:' ',$
                 objname:' ',$
                 z:0.0,$
                 ra:0.0d,$
                 dec:0.0d,$
                 sn:0.0,$
                 mjd:0.0,$
                 airmass:0.0,$
                 rand_error:0.0}

   stars = replicate({stmp}, nfiles*100)


  n=0
  for i=0,n_elements(zfiles)-1 do begin


     ; DETERMINE ASSOCIATED MCERRS FILES
     s = strsplit(zfiles[i],'.',/extract)
     mcfile = path+'mc_'+s[1]+'.fits'

     ; READ FILES
     x  = mrdfits(zfiles[i],1,h)
     mc = mrdfits(mcfile,1)


     ; KEEP GOOD VELOCITIES, NO SERENDIPS
     q = where(x.zquality ge 3 and strmatch(x.objname,'serendip*') ne 1,nc)
     x=x[q]
     mc=mc[q]

     ; CONVERT RA/DEC TO DECIMAL
     ra = dblarr(nc)
     dec=dblarr(nc)
     for j=0,nc-1 do begin
        print,x[j].ra,x[j].dec
        sradec,x[j].ra,x[j].dec,r,d
        ra[j]=r
        dec[j]=d

     endfor

     
     
     ; CONVERT TO RANDOM ERROR
     rand = sqrt(mc.mcerr_vcorr^2 - 2.2^2)

     stars[n:n+nc-1].subclass = x.subclass
     stars[n:n+nc-1].objname = x.objname
     stars[n:n+nc-1].z = x.z
     stars[n:n+nc-1].sn = x.sn
     stars[n:n+nc-1].mjd = x.mjd
     stars[n:n+nc-1].ra = ra
     stars[n:n+nc-1].dec = dec 
     stars[n:n+nc-1].rand_error = rand


     n=n+nc
  endfor
  q=where(stars.z ne 0,nstars)
  stars=stars[q]
  print,nstars
  print

;******* NOW SEARCH FOR REPEATS  **********
  
 ; CREATE TMP STRUCTURE FOR ALL DATA
  temp = {rtmp,  vdiff:0.0,$
                 subclass1:' ',$
                 objname1:' ',$
                 z1:0.0,$
                 ra1:0.0,$
                 dec1:0.0,$
                 sn1:0.0,$
                 mjd1:0.0,$
                 airmass1:0.0,$
                 rand_error1:0.0,$
                 subclass2:' ',$
                 objname2:' ',$
                 z2:0.0,$
                 ra2:0.0,$
                 dec2:0.0,$
                 sn2:0.0,$
                 mjd2:0.0,$
                 airmass2:0.0,$
                 rand_error2:0.0}
   repeats = replicate({rtmp}, fix(nstars/2.))

   nr = 0
   for i=0,nstars-1 do begin

      if stars[i].ra ne 0 then begin
         spherematch,stars.ra,stars.dec,stars[i].ra,stars[i].dec,1./3600,m1,m2,maxmatch=10


         ; FIND REPEATS AND ADD TO STRUCTURE
         if n_elements(m1) ge 2 then begin
            print,stars[m1[0]].ra,stars[m1[1]].ra
            print,stars[m1[0]].z*3e5,stars[m1[1]].z*3e5
            print
            add_to_repeats,repeats,stars,nr,m1[0],m1[1]
            nr =nr+1
            
            if n_elements(m1) ge 3 then begin
               print,stars[m1[0]].ra,stars[m1[2]].ra
               print,stars[m1[0]].z*3e5,stars[m1[2]].z*3e5
               print
               add_to_repeats,repeats,stars,nr,m1[0],m1[2]
               nr =nr+1
               if n_elements(m1) ge 4 then begin
                  print,stars[m1[0]].ra,stars[m1[3]].ra
                  print,stars[m1[0]].z*3e5,stars[m1[3]].z*3e5
                  print
                  add_to_repeats,repeats,stars,nr,m1[0],m1[3]
                  nr =nr+1
               endif
            endif 
         endif 
         stars[m1].ra = 0  ; REMOVE THESE STARS FROM SEARCH STRUCTURE
      endif 
      
   endfor

; WRITE REPEAT ARRAY
   q=where(repeats.ra1 ne 0)
   repeats=repeats[q]
   plothist,repeats.vdiff


   mwrfits,repeats,'repeats_oldtemplates.fits',/create
   
end
