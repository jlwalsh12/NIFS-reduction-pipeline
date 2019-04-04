pro find_offsets_lp, dir, files
  
;uses the collapsed galaxy data cubes to determine the spatial offsets
;between each cube. routine calculates the offsets in 3 ways: by
;finding the location of the maximum flux value of the image, by
;fitting a 2D gaussian to the image to get the central location, and
;by cross-correlating the images.
  
if n_elements(files) GT 1 then begin
   print, ''
   print, 'Now calculating the spatial offsets between the collapsed '+$
          'data cubes, using the first file as the reference, in 3 ways:'
   print, '(1) finding the maximum value of the image'
   print, '(2) fitting a 2D gaussian to the image'
   print, '(3) cross-correlating the images'
   print, ''
endif else begin
   print, ''
   print, 'Only one collapsed data cube. Now:'
   print, '(1) finding the maximum value of the image'
   print, '(2) fitting a 2D gaussian to the image'
   print, ''
endelse

;read in the collapsed data cubes
img_ref = readfits(dir+files[0], head, /silent)
img_all = fltarr((size(img_ref))[1], (size(img_ref))[2], n_elements(files))
for i = 0, n_elements(files)-1 do begin
   img_all[*,*,i] = readfits(dir+files[i], /silent)
endfor

;array that will hold the centers
cen_max_final = fltarr(2, n_elements(files))
cen_gauss_final = fltarr(2, n_elements(files))
cen_crosscorr_final = fltarr(2, n_elements(files))

for i = 0, n_elements(files)-1 do begin

   img_tmp = reform(img_all[*,*,i])

   ;first find the spaxel with the maximum value
   maxval = max(img_tmp, maxindex)
   ncol = (size(img_tmp))[1]
   xmax = maxindex MOD ncol & ymax = maxindex / ncol
   cen_max_final[0,i] = xmax & cen_max_final[1,i] = ymax

   ;now try fitting a 2D Gaussian. trying to fit the entire img, seems
   ;to often give errors, focus near the xmax, ymax instead.
   xtrimdown = xmax-15 > 0
   xtrimup = xmax+15 < (size(img_tmp))[1] - 1
   ytrimdown = ymax-15 > 0
   ytrimup = ymax+1 < (size(img_tmp))[2] - 1
   yfit = mpfit2dpeak(img_tmp[xtrimdown:xtrimup,ytrimdown:ytrimup], $
                      bestparam, /gaussian)
   xcen_gauss = bestparam[4]+xtrimdown & ycen_gauss = bestparam[5]+ytrimdown
   cen_gauss_final[0,i] = xcen_gauss & cen_gauss_final[1,i] = ycen_gauss

end

;calculate the offset relative to the first file
offsets_max_final = fltarr(2, n_elements(files))
offsets_gauss_final = fltarr(2, n_elements(files))
if n_elements(files) GT 1 then begin
   for i = 1, n_elements(files)-1 do begin
      offsets_max_final[0,i] = cen_max_final[0,0] - cen_max_final[0,i]
      offsets_max_final[1,i] = cen_max_final[1,0] - cen_max_final[1,i]
      offsets_gauss_final[0,i] = cen_gauss_final[0,0] - cen_gauss_final[0,i]
      offsets_gauss_final[1,i] = cen_gauss_final[1,0] - cen_gauss_final[1,i]
   endfor
endif else begin
   offsets_max_final[0,0] = cen_max_final[0,0] - cen_max_final[0,0]
   offsets_max_final[1,0] = cen_max_final[1,0] - cen_max_final[1,0]
   offsets_gauss_final[0,0] = cen_gauss_final[0,0] - cen_gauss_final[0,0]
   offsets_gauss_final[1,0] = cen_gauss_final[1,0] - cen_gauss_final[1,0]
endelse

;find the RA and Dec in decimal degrees of the reference image center
;pixel
extast, head, astr, noparams
if noparams EQ -1 then begin
   print, 'Problem with the header of the reference image.'
   stop
endif
;perform one more check to make sure the header is correct.
checkastr = strcompress(string(astr.ctype[0]), /remove_all)
if checkastr NE 'RA---TAN' then begin
   print, 'Problem with the header of the reference image.'
   stop
endif
xy2ad, cen_max_final[0,0], cen_max_final[1,0], astr, ra_deg_max, dec_deg_max
xy2ad, cen_gauss_final[0,0], cen_gauss_final[1,0], astr, ra_deg_gauss, $
       dec_deg_gauss

;calculate the offsets using a cross-correlation
if n_elements(files) GT 1 then begin
   
   offsets_crosscorr_final = fltarr(2, n_elements(files))
   for i = 1, n_elements(files)-1 do begin
   
      correl_optimize, reform(img_all[*,*,0]), reform(img_all[*,*,i]), $
                       xoffset_optimum, yoffset_optimum, $
                       XOFF_INIT = 0,   $
                       YOFF_INIT = 0,   $
                       /NUMPIX, $
                       MAGNIFICATION = 10., $
                       PLATEAU_TRESH = 0.01
   
      offsets_crosscorr_final[0,i] = xoffset_optimum
      offsets_crosscorr_final[1,i] = yoffset_optimum
   
   endfor
   ;don't find the absolute x and y centers, just relative with the
   ;cross-correlation method. so will just adopt ra_deg_gauss,
   ;dec_deg_gauss.
endif

;get the pixel scale to include in comments of offset file.
pixscale = sxpar(head, 'PIXSCALE')

;print the results of each offset method to output files
forprint, offsets_max_final[0,*], offsets_max_final[1,*], format='(2f15.3)', $
          textout=dir+'offset_max.list', /silent, $
          comment='#x,y offsets (pixels, pix='+$
          strcompress(string(pixscale),/remove_all)+'") needed to match '+$
          'first entry'
forprint, offsets_gauss_final[0,*], offsets_gauss_final[1,*], $
          format='(2f15.3)', textout=dir+'offset_gauss.list', /silent, $
          comment='#x,y offsets (pixels, pix='+$
          strcompress(string(pixscale),/remove_all)+'") needed to match '+$
          'first entry'
if n_elements(files) GT 1 then $
   forprint, offsets_crosscorr_final[0,*], offsets_crosscorr_final[1,*], $
             format='(2f15.3)', textout=dir+'offset_crosscorr.list', /silent, $
             comment='#x,y offsets (pixels, pix='+$
             strcompress(string(pixscale),/remove_all)+'") needed to match'+$
             ' first entry'

;print the x and y reference pixels (i.e., crpix1 and crpix2) and the
;ra and dec (i.e., crval1 and crval2). need to conform to the fits
;standard where indexing starts at 1 and not 0 for crpix1 and crpix2.
forprint, cen_max_final[0,0] + 1., cen_max_final[1,0] + 1., $
          ra_deg_max, dec_deg_max, format='(4f15.6)', $
          textout=dir+'headervals_offset_max.dat', $
          comment='#CRPIX1, CRPIX2 (in fits standard, '+$
          'index starts at 1), CRVAL1 (RA deg), CRVAL2 '+$
          '(Dec deg) for reference image', /silent
forprint, cen_gauss_final[0,0] + 1., cen_gauss_final[1,0] + 1., $
          ra_deg_gauss, dec_deg_gauss, format='(4f15.6)', $
          textout=dir+'headervals_offset_gauss.dat', $
          comment='#CRPIX1, CRPIX2 (in fits standard, '+$
          'index starts at 1), CRVAL1 (RA deg), CRVAL2 '+$
          '(Dec deg) for reference image', /silent
if n_elements(files) GT 1 then $
   forprint, cen_gauss_final[0,0] + 1., cen_gauss_final[1,0] + 1., $
             ra_deg_gauss, dec_deg_gauss, format='(4f15.6)', $
             textout=dir+'headervals_offset_crosscorr.dat', $
             comment='#CRPIX1, CRPIX2 (in fits standard, '+$
             'index starts at 1), CRVAL1 (RA deg), CRVAL2 '+$
             '(Dec deg) for reference image', /silent

end
