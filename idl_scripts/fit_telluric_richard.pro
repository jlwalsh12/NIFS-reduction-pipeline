pro fit_telluric_richard, workdir, telluric_file, refdir

;Create the 1D telluric correction spectrum using the combined 1D
;spectrum of a (A0 V) telluric star. Remove the absorption line(s) and
;the blackbody shape from the telluric star spectrum. This is done by
;fitting Kurucz models of Vega and A-stars with pPXF.
  
;change the screen size for the PPXF fitting.
r = GET_SCREEN_SIZE()
window, 0, xsize=r[0]*0.4, ysize=r[1]*0.3
  
for i=0, n_elements(telluric_file)-1 do begin

   ;==========
   ; Read in the NIFS telluric standard
   ;=========
   ;
   print, 'Modifying the NIFS spectrum for '+$
          strcompress(string(telluric_file[i]),/remove_all)
   spec0 = mrdfits(workdir+telluric_file[i],0,h0,/silent)
   spec1 = mrdfits(workdir+telluric_file[i],1,h1,/silent)
   spec2 = mrdfits(workdir+telluric_file[i],2,h2,/silent)
   spec3 = mrdfits(workdir+telluric_file[i],3,h3,/silent)
   grating = strmid(sxpar(h0,'GRATING'),0,1)

   ;=========
   ; Template used depends on spectral range. Visual inspection
   ; of fits suggests that Vega alone is good for all gratings
   ; except K, where some small benefit comes from the library of
   ; Astar models. These tend to have other weak lines that create
   ; spurious 'emission' peaks in the final telluric spectrum.
   ;===========
   ;
   case grating of
      'K': libfile = refdir+'kurucz_resamp_all.lib'
      else: libfile = refdir+'kurucz_resamp_vega.lib'
   endcase

   ;==========
   ; Create lambda vector
   ;==========
   ;
   lstar = sxpar(h1,'CRVAL1') + [findgen(sxpar(h1,'NAXIS1'))+1 - $
                                 sxpar(h1,'CRPIX1')] * sxpar(h1,'CD1_1')
   wlim = minmax(lstar)

   ;==========
   ; Oversample NIFS spectrum
   ;==========
   ;
   factor = 2.
   npix = sxpar(h1,'NAXIS1')
   rstar = rebin(spec1,factor*npix)
   rlstar = rebin(lstar,factor*npix)

   ;==========
   ; Log-rebin the telluric standard
   ;==========
   ;
   log_rebin, minmax(rlstar), rstar, temp1, lnrlstar, VELSCALE=velscale
   nlnpix = n_elements(temp1)
   
   ;==========
   ; Get the hi-res reference spectrum and log-rebin.
   ;=========
   ;
   readcol,libfile,temps,f='a',/silent,comment="#"
   ntemp = n_elements(temps)
   
   for j = 0, ntemp-1 do begin
      ; Read resampled spectra, already in Angstrom (vaccuum wavelengths)
      readcol,refdir+'kurucz_stars/'+temps[j],lsun,fsun,/silent
      ; Truncate to useful region
      ind = where((lsun ge wlim[0]) AND (lsun le wlim[1]))
      ind = [ind[0]-1,ind, max(ind)+1]
      lsun = lsun[ind]
      fsun = fsun[ind]
      ; Log-rebin the reference spectrum
      log_rebin, minmax(lsun), fsun, temp2, lnlsun, VELSCALE=velscale
      if j eq 0 then template = fltarr(n_elements(temp2),ntemp)
      template[*,j] = temp2
   endfor
   

   ;==========
   ; Establish 'goodpixels'
   ;==========
   ;
   s = size(temp2)
   goodpixels = [cap_range(30,s[1]-30)]
   case grating of
      "Z":goodpixels = [cap_range(900,1300),cap_range(2600,3000)]+100
      "J":goodpixels = [cap_range(30,s[1]-30)]
      "H":goodpixels = [cap_range(30,s[1]-30)]
      ;"K":goodpixels = [cap_range(1300,1800)]
      "K":goodpixels = [cap_range(1300,1950)]
      else:
   endcase
   
   ;==========
   ; Run pPXF; the template spectra are flux calibrated
   ; and we are fitting over a small wavelength range. Set
   ; the multiplicative polynomial to degree 0 and use no
   ; additive polynomial.
   ;==========
   ;
   mpolynom = 0
   polynom = -1
   nmom = 2
   vdif = 2.99e5 * (min(lnlsun)-min(lnrlstar))

   ppxf, template, temp1, replicate(1.D,nlnpix), velScale, $
         [vdif+100.0, 20.0], sol, mdegree=mpolynom, degree=polynom, $
         moments=nmom, bestfit=bestfit, goodpixels=goodpixels, /plot,$
         weights=weights, bias=0.
   
   ;==========
   ; Take the ratio of the fit. The result is the telluric absorption
   ; spectrum.
   ;==========
   ;
   abs = temp1 / bestfit

   ;==========
   ; Interpolate to go back to linear lambda and to the original binning
   ;==========
   ;
   dw = (max(lstar)-min(lstar))/(n_elements(lstar)-1) 
   abs_org = interpol(abs,exp(lnrlstar),lstar-dw/2.0,/spline)
   fit_org = interpol(bestfit,exp(lnrlstar),lstar-dw/2.0,/spline)
   
   ;=========
   ; Scale the variance
   ;=========
   abs_var = spec2 / fit_org^2

   ;========
   ; Write total weights to the header - useful for preliminary
   ; flux calibration
   ;=========
   ;
   sxaddpar,h0,'WEIGHT',total(weights),'ergs/cm^2/s/ster/nm'

   ;==========
   ; Write the spectrum to disk
   ;==========
   print,'Writing telluric absorption spectrum to disk.'
   tmp = strsplit(telluric_file[i],'.',/extract)
   outfile = 'c'+tmp[0]+'_final.fits'
   mwrfits,spec0,workdir+outfile,h0,/create
   mwrfits,abs_org,workdir+outfile,h1,/silent
   mwrfits,abs_var,workdir+outfile,h2,/silent
   mwrfits,spec3,workdir+outfile,h3,/silent

endfor

end
