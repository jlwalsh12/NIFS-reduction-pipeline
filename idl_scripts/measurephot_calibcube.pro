pro measurephot_calibcube, teldir, telfile, rootdir, reducedir, $
                           galaxy, obs_setup, telstar_name

;Determines the multiplicative conversion factor to go from counts/s
;to the known flux density at the isophotal wavelength for a telluric
;star (given in
;rootdir+reducedir+'tellurics/tellurics_flux.dat'). This conversion
;factor is then applied to each spectra of the merged galaxy cube and
;the merged sky cube.


;read in the header and the telluric cube
junk=readfits(teldir+telfile,head0,ext=0)
im=readfits(teldir+telfile,head,ext=1)

;set the wavelength vector.
imsize=size(im,/dim)
lambda0=sxpar(head,'CRVAL3')
dlambda=sxpar(head,'CD3_3')
crpix=sxpar(head,'CRPIX3')
lambda=( (dindgen(imsize[2])+1-crpix)*dlambda ) + lambda0
nlambda=n_elements(lambda)

;get the exposure time
exptime = sxpar(head0,'EXPTIME')

;collapse the telluric cube
collim = dblarr(imsize[0],imsize[1])
for i=0, imsize[0]-1 do begin
   for j=0, imsize[1]-1 do begin
      collim[i,j] = total(im[i,j,*],/double)
   endfor
endfor

;fit a 2D Gaussian to the collapsed telluric cube to get the x,y
;center.
;first find the spaxel with the maximum value
maxval = max(collim, maxindex)
ncol = imsize[0]
xmax = maxindex MOD ncol & ymax = maxindex / ncol
;trying to fit the entire img, seems to often give errors, focus near
;the xmax, ymax instead.
xtrimdown = xmax-15 > 0
xtrimup = xmax+15 < imsize[0] - 1
ytrimdown = ymax-15 > 0
ytrimup = ymax+1 < imsize[1] - 1
yfit = mpfit2dpeak(collim[xtrimdown:xtrimup,ytrimdown:ytrimup], $
                   bestparam, /gaussian)
xcen_gauss = bestparam[4]+xtrimdown
ycen_gauss = bestparam[5]+ytrimdown

;extract a 1D spectrum from the telluric cube using a 1.5"
;diameter. divide by exposure time to get counts/s.
spec=dblarr(nlambda)
for i=50, nlambda-50 do begin
   aper,im[*,*,i],xcen_gauss,ycen_gauss,flux,fluxerr,sky,skyerr,$
        1.0,[15.],setskyval=0.0,/exact,/flux,/silent,/nan
    spec[i]=flux[0]
endfor
spec=spec/exptime

;get the flux density of the telluric star at the effective/isophotal
;wavelength.
;isophotal wavelength for Ks, in Angstroms
isolambda=2.159D4 
;get the k band mag from simbad
querysimbad, telstar_name, ra, dec, kmag=kmag
;target flux in Jy
jyflux=10.D0^(-0.4D*kmag)*666.7D
;target flux in in units of ergs/s/cm^2/Angstrom
targetflux=jyflux*1.D-23*2.99792D18/(isolambda)^2

;calculate the ratio between the known flux (ergs/s/cm^2/A) and the
;pixel flux (counts/s) at the isophotal wavelength. need to
;interpolate to get values at the isophotal wavelength.
lowind=where(lambda GT 2.02D4 AND lambda LT 2.14D4)
highind=where(lambda GT 2.19D4 AND lambda LT 2.4D4)
ind=[lowind,highind]
fitpar=poly_fit(lambda[ind],spec[ind],3)
fit=poly(lambda,fitpar)
;plot, lambda,spec
;oplot,lambda,fit,linestyle=2,color=200
;wait, 5.
targetval=interpol(fit,lambda,isolambda)
conversion=targetflux/targetval

print, 'conversion = ', conversion

;plot the calibrated 1D telluric spectrum
print, 'Plotting the calibrated 1D telluric spectrum. Cross shows'+$
       ' the input flux density at the isophotal wavelength.'
plot, lambda, spec*conversion, psym=10
plots,!X.CRANGE,[targetflux,targetflux]
plots,[isolambda,isolambda],!Y.CRANGE
print, ''
;wait, 5.

;now apply the conversion to the merged galaxy and sky cubes.

;start with the galaxy cube

;copy the galaxy infile to the outfile. will modify the outfile.
spawn, 'cp '+rootdir+reducedir+galaxy+'/merged/'+obs_setup+'/'+$
       galaxy+'_combined.fits '+rootdir+reducedir+galaxy+$
       '/merged/'+obs_setup+'/'+galaxy+'_combined_flux.fits'

;read in the main header, the cube, and the variance cube
cube = readfits(rootdir+reducedir+galaxy+'/merged/'+obs_setup+$
                '/'+galaxy+'_combined.fits', exten_no=1, head)
var = readfits(rootdir+reducedir+galaxy+'/merged/'+obs_setup+$
               '/'+galaxy+'_combined.fits', exten_no=2)


;flux calibrate each spectrum of the cube
outcube = cube
outvar = var
for i=0, (size(cube))[1]-1 do begin
   for j=0, (size(cube))[2]-1 do begin
      outcube[i,j,*] = cube[i,j,*] * conversion
      outvar[i,j,*] = var[i,j,*] * conversion^2
   endfor
endfor


;modify the output fits file
modfits, rootdir+reducedir+galaxy+'/merged/'+obs_setup+'/'+$
         galaxy+'_combined_flux.fits', outcube, exten_no=1
modfits, rootdir+reducedir+galaxy+'/merged/'+obs_setup+'/'+$
         galaxy+'_combined_flux.fits', outvar, exten_no=2


;repeat for the sky cube if one exists

exist_test = file_test(rootdir+reducedir+galaxy+'/merged/'+obs_setup+'/'+$
                       galaxy+'_sky_combined.fits')

if exist_test EQ 1 then begin

   ;copy the sky infile to the outfile. will modify the outfile.
   spawn, 'cp '+rootdir+reducedir+galaxy+'/merged/'+obs_setup+'/'+$
          galaxy+'_sky_combined.fits '+rootdir+reducedir+galaxy+$
          '/merged/'+obs_setup+'/'+galaxy+'_sky_combined_flux.fits'

   ;read in the main header, the cube, and the variance cube
   cube = readfits(rootdir+reducedir+galaxy+'/merged/'+obs_setup+$
                   '/'+galaxy+'_sky_combined.fits', exten_no=1)
   var = readfits(rootdir+reducedir+galaxy+'/merged/'+obs_setup+$
                  '/'+galaxy+'_sky_combined.fits', exten_no=2)

   ;flux calibrate each spectrum of the cube
   outcube = cube
   outvar = var
   for i=0, (size(cube))[1]-1 do begin
      for j=0, (size(cube))[2]-1 do begin
         outcube[i,j,*] = cube[i,j,*] * conversion
         outvar[i,j,*] = var[i,j,*] * conversion^2
      endfor
   endfor

   ;modify the output fits file
   modfits, rootdir+reducedir+galaxy+'/merged/'+obs_setup+'/'+$
            galaxy+'_sky_combined_flux.fits', outcube, exten_no=1
   modfits, rootdir+reducedir+galaxy+'/merged/'+obs_setup+'/'+$
            galaxy+'_sky_combined_flux.fits', outvar, exten_no=2
endif

end
