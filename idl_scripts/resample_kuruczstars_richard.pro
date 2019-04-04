pro resample_kuruczstars_richard
  
; Resamples the Kurucz spectra for easier log_rebin call in
; fit_telluric_richard.pro. Don't need to run, as the resampled
; spectra are included as part of the NIFS pipeline download, but this
; code is included for completeness.
;
; Kurucz spectra from: kurucz.harvard.edu/stars/


libfile = '/Users/jlwalsh/Data/LP_2016/nifs_reduction_info/full_reduce/'+$
          'nifs_pipeline_v4/ref_files/kurucz_all.lib'
; Ang/pix
stepsize = 1.0
; Useful range for NIFS, in Ang. Also sets new wavelength sampling,
; together with stepsize
wlim = [9400.,25000.]

readcol,libfile,libdir,f='a',numline=1,/silent
readcol,libfile,files,f='a',skipline=1,/silent
nfiles = n_elements(files)

for i=0,nfiles-1 do begin
   
   readcol,libdir[0]+files[i],lsun,fsun,/silent
   ; Put into Angstrom
   lsun = 10.0*lsun
   ; Truncate to useful region
   ind = where((lsun ge wlim[0]) AND (lsun le wlim[1]))
   ind = [ind[0]-1,ind,max(ind)+1]
   lsun = lsun[ind]
   fsun = fsun[ind]
   
   ; Rebin to something more reasonable, but still higher than NIFS.
   ; This saves time in log-rebin
   npix = FLOOR(1.0 + (max(wlim)-min(wlim))/stepsize)
   lnew = min(wlim) + FINDGEN(npix)*stepsize
   print,lnew[0],npix,minmax(lnew)
   fnew = interpol(fsun, lsun, lnew, /spline)

   ; Write out the file
   forprint,lnew,fnew,textout=libdir[0]+files[i]+'_resamp',/silent,$
            comment='#Vac. Wavelength (ang), Flux (ergs/cm^2/s/ster/nm)'

endfor

end
