pro measure_lsf_lp, workdir, skycube, refdir, galaxy

;this routine measures the line-spread-function using a sky cube. at
;each spatial location, 13 OH lines between 2.0 and 2.25 microns are
;fit with gaussians, and the median FWHM (in Angstroms) of the lines
;for that spaxel is saved to an output fits file.
  
COMMON skyblock, gauss_sep

;read in the merged sky cube
cube=readfits(workdir+skycube,head,exten_no=1,/silent)
cubesize=size(cube,/dim)
;set the wavelength
lambda0=sxpar(head,'CRVAL3')
dlambda=sxpar(head,'CD3_3')
crpix=sxpar(head,'CRPIX3')
lambda=( (dindgen(cubesize[2])+1-crpix)*dlambda ) + lambda0

;conversion between FWHM and dispersion
sigma2fwhm=2.*Sqrt(2.*alog(2.))

;read in the file giving the rest wavelengths (vaccuum) for the OH
;lines
readcol,refdir+'rousselot2000.dat',skylambda,skyint,format='F,F',$
        skip=28, /silent
;these lines all have close doublets except 21955 line and are
;relatively isolated
uselambda=double([20339.,20412.7,20563.5,20729.0,21176.,21802.3,21955.6,$
                  22052.,22125.5,22247.,22312.7,22460.0,22517.9])
nlines=n_elements(uselambda)
reflambda=dblarr(nlines,2)
for i=0,nlines-1 do begin
   ind=where(abs(skylambda-uselambda[i]) LT 2. AND skyint GT 1.,$
             nind)
   if (nind NE 2) then begin
      print, 'Problem finding sky line in the reference list.'
      stop
   endif
   reflambda[i,*]=skylambda[ind]
endfor

;set speed of light and create arrays that will hold the fit results
c=2.99792458D5
intmap=dblarr(cubesize[0],cubesize[1],nlines)
velmap=dblarr(cubesize[0],cubesize[1],nlines)
fwhmmap=dblarr(cubesize[0],cubesize[1],nlines)

minlcont=-30.
maxlcont=-20.
minucont=20.
maxucont=30.
fitwidth=15.

loadrgb
;loop over sky lines
for k=0,nlines-1 do begin

   centerwave=mean(reflambda[k,*])
   gauss_sep=reflambda[k,1]-reflambda[k,0]
   lind=where(lambda GE centerwave+minlcont AND $
              lambda LE centerwave+maxlcont,nlind)
   uind=where(lambda GE centerwave+minucont AND $
              lambda LE centerwave+maxucont,nuind)
   fitind=where(lambda GE centerwave-fitwidth AND $
                lambda LE centerwave+fitwidth,nfitind)

   ;loop over spatial dimensions. i for x and j for y.
   for i=0,cubesize[0]-1 do begin
      for j=0,cubesize[1]-1 do begin

         ;the spectrum for this spaxel. multiply by 1.D18 to avoid
         ;very small numbers from flux calibration.
         spec=cube[i,j,*] * 1.D18
         ;the continuum level (a constant)
         background=median([spec[lind],spec[uind]])
         ;attempt to filter out obviously poor specta where the line
         ;is not present
         get_element, lambda, centerwave, tmp_ind
         if spec[tmp_ind] LE 0. OR spec[tmp_ind-1] LE 0. OR $
            spec[tmp_ind+1] LE 0. OR $
            spec[tmp_ind] LE background then goto, SKIPTHIS
         ;set the initial guesses for the 4 parameters
         parinfo = replicate({value:0.D, fixed:0, limited:[0,0],$
                              limits:[0.D,0]}, 4)
         parinfo(*).value = [max(spec[fitind]),centerwave,$
                             4.2D/sigma2fwhm,background]
         ;assume constant errors for now (until have a correct
         ;variance spectrum)
         specerr = replicate(1.D,nfitind)

         ;fit a double gaussian, where the separation between the
         ;components is known, and the amplitude and dispersion is the
         ;same
         param = mpfitfun('gauss_binary_mpfit',lambda[fitind],$
                          spec[fitind],specerr,parinfo=parinfo,$
                          yfit=yfit,/quiet)

         ;store fit results
         intmap[i,j,k]=param[0]*param[2]
         velmap[i,j,k]=(param[1]-centerwave)/centerwave*c
         fwhmmap[i,j,k]=param[2]*sigma2fwhm

         ;another attempt to filter out poor spectra/fits
         if param[0] LE 0. OR param[2]*sigma2fwhm GE 10. OR $
            ((param[1]-centerwave)/centerwave*c) GE 10. then begin
            intmap[i,j,k]=0.D
            velmap[i,j,k]=0.D
            fwhmmap[i,j,k]=0.D
         endif
         
         ;plot some example fits
         titstring=string(i)+'_'+string(j)+' '
         if (j EQ 30) then begin
            plot,lambda,spec,title=titstring,$
                 xrange=[centerwave+minlcont,centerwave+maxucont],/xsty,psym=10
            oplot,lambda[fitind],yfit,color=1,psym=10
            wait,0.1
         endif

         SKIPTHIS:
      endfor
   endfor

   ;normalize the integrated flux of the line
   intmap[*,*,k]=intmap[*,*,k]/median(intmap[*,*,k])

endfor

;print some diagnostic information to the screen.
med_velmap = median(velmap,dim=3)
med_velmap_trim = med_velmap[3:(size(med_velmap))[1]-3,$
                             3:(size(med_velmap))[2]-3]
print, 'Excluding the edges of the cube, the median velocity '+$
       'shift is: '+ $
       strcompress(string(median(abs(med_velmap_trim))),/remove_all)+$
       ' km/s.'

line_fwhm = dblarr(nlines)
fwhmmap_trim = fwhmmap[3:(size(fwhmmap))[1]-3,3:(size(fwhmmap))[2]-3,*]
for k=0,nlines-1 do line_fwhm[k] = median(fwhmmap_trim[*,*,k])
print, 'Excluding the edges of the cube, there is a range of '+$
       strcompress(string(max(line_fwhm)-min(line_fwhm)),/remove_all)+' ang '+$
       'in the FWHM of different lines.'

;write the output LSF map
med_fwhmmap = median(fwhmmap,dim=3)
writefits,workdir+galaxy+'_fullskycube_fwhm_med.fits',med_fwhmmap

end
