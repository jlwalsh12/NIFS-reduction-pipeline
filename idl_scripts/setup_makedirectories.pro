pro setup_makedirectories, rootdir, datadir, workdir, dates, $
                           galaxies, tellurics

;Program to create directory structure, sort files into the correct
;directories, create the daycal file lists (i.e., flatlist,
;flatdarklist), galaxy file lists, telluric file lists, and psf file
;lists.


;make main directory structure
;------------------------------

;for workdir
if (file_test(workdir,/dir) EQ 0) then spawn, 'mkdir '+workdir

;for daycals
if (file_test(workdir+'daycals',/dir) EQ 0) then $
   spawn, 'mkdir '+workdir+'daycals'

;for tellurics
if (file_test(workdir+'tellurics',/dir) EQ 0) then $
   spawn, 'mkdir '+workdir+'tellurics'

;for object exposures
for i = 0, n_elements(galaxies)-1 do begin
   if (file_test(workdir+galaxies[i],/dir) EQ 0) then $
      spawn, 'mkdir '+workdir+galaxies[i]
   if (file_test(workdir+galaxies[i]+'/merged',/dir) EQ 0) then $
      spawn, 'mkdir '+workdir+galaxies[i]+'/merged'
endfor

;for PSF star. start by making this directory, and delete later if no
;PSF star observation was made.
if (file_test(workdir+'psfs',/dir) EQ 0) then $
   spawn, 'mkdir '+workdir+'psfs'

;-----------------------

;loop over each night
for i=0, n_elements(dates)-1 do begin

   cd, datadir+dates[i]+'/'

   readcol,'datalist.txt', filestem, object, obstype, obsclass, $
           exptime, filter, centralwave, aperture, xoff, yoff, $
           airmass, pa, medianim, obsid, ut, $
           FORMAT='A,A,A,A,I,A,F,A,F,F,F,F,F,A,A', /silent


   ;make the dated subdirectories
   ;------------------------------
   
   ;for the daycals
   if (file_test(workdir+'daycals/'+dates[i],/dir) EQ 0) then $
      spawn, 'mkdir '+workdir+'daycals/'+dates[i]


   ;for the tellurics
   if (file_test(workdir+'tellurics/'+dates[i],/dir) EQ 0) then $
      spawn, 'mkdir '+workdir+'tellurics/'+dates[i]


   ;for the objects
   for j = 0, n_elements(galaxies)-1 do begin

      ;see if the full galaxy name listed by user matches the object
      ;name in the data file.
      test = strmatch(object, galaxies[j], /fold_case)
      index_tmp = where(test GT 0, count_tmp)
      ;if the above does not work, see if the end of the galaxy name
      ;listed by the user matches some part of the object name in the
      ;data file.
      if count_tmp EQ 0 then begin
         catalog_firsttwoletters = strmid(galaxies[j],0,2)
         case catalog_firsttwoletters of
            'ng': galaxies_end = strsplit(galaxies[j],'ngc',/extract)
            'ug': galaxies_end = strsplit(galaxies[j],'ugc',/extract)
            'mr': galaxies_end = strsplit(galaxies[j],'mrk',/extract)
            'pg': galaxies_end = strsplit(galaxies[j],'pgc',/extract)
            'ic': galaxies_end = strsplit(galaxies[j],'ic',/extract)
            else: begin
               print, 'Cannot match galaxy name given by user to '+$
                      'object name in any fits file headers.'
               stop
            end
         endcase
         test = strmatch(object, '*'+galaxies_end+'*', /fold_case)
      endif
      
      ;if the galaxy files have been located and are science exposures
      ;(not acquisition exposures) then create the subdirectory
      index_galfiles=where(test GT 0 AND obsclass EQ 'science' AND $
                           exptime GE 180., count_galfiles)
      if count_galfiles GT 0 then begin
         if (file_test(workdir+galaxies[j]+'/'+dates[i],/dir) EQ $
             0) then spawn, 'mkdir '+workdir+galaxies[j]+'/'+$
                            dates[i]
      endif
      
   endfor

   ;for the psf. again will delete later if there are no PSF star
   ;observations
   if (file_test(workdir+'psfs/'+dates[i],/dir) EQ 0) then $
      spawn, 'mkdir '+workdir+'psfs/'+dates[i]
   if (file_test(workdir+'psfs/merged',/dir) EQ 0) then $
      spawn, 'mkdir '+workdir+'psfs/merged'

   ;------------------------------

   
   ;loop over different central wavelengths
   uniqcentralwave = centralwave[uniq(centralwave)]
   nuniqcentralwave = n_elements(uniqcentralwave)

   for j=0, nuniqcentralwave-1 do begin

      index = where(centralwave EQ uniqcentralwave[j],nindex)
      if nindex EQ 0 then begin
         print, 'Problem finding files with the same central '+$
                'wavelength.'
         stop
      endif
      filestem_tmp = filestem[index]
      object_tmp = object[index]
      obstype_tmp = obstype[index]
      obsclass_tmp = obsclass[index]
      exptime_tmp = exptime[index]
      filter_tmp = filter[index]
      centralwave_tmp = centralwave[index]
      aperture_tmp = aperture[index]
      xoff_tmp = xoff[index]
      yoff_tmp = yoff[index]
      airmass_tmp = airmass[index]
      pa_tmp = pa[index]
      medianim_tmp = medianim[index]
      obsid_tmp = obsid[index]
      ut_tmp = ut[index]


      ;make the observational setup subdirectories. for the tellurics,
      ;galaxies, and psf, sort files into the subdirectories. make the
      ;object and sky list files for the tellurics, galaxies, and psf.
      ;---------------------------------------------------------------

      ;use filter_tmp[0] because all elements of filter_tmp will be
      ;the same given this particular central wavelength
      setup = strcompress(strlowcase(filter_tmp[0]),/remove_all)+'_'+$
              strcompress(sigfig(uniqcentralwave[j],3),$
                                /remove_all)

      ;for the daycals      
      if (file_test(workdir+'daycals/'+dates[i]+'/'+setup,/dir) EQ $
          0) then spawn, 'mkdir '+workdir+'daycals/'+dates[i]+'/'+$
                         setup

      
      ;for the tellurics
      if (file_test(workdir+'tellurics/'+dates[i]+'/'+setup,/dir) EQ $
          0) then spawn, 'mkdir '+workdir+'tellurics/'+dates[i]+'/'+$
                         setup

      for k = 0, n_elements(tellurics)-1 do begin

         ;see if the full telluric name listed by user matches the
         ;object name in the data file.
         test = strmatch(object_tmp, tellurics[k], /fold_case)
         index_tmp = where(test GT 0, count_tmp)
         ;if the above does not work, see if the end of the telluric
         ;name listed by the user matches some part of the object name
         ;in the data file.
         if count_tmp EQ 0 then begin
            catalog_firsttwoletters = strmid(tellurics[k],0,2)
            case catalog_firsttwoletters of
               'hi':tellurics_end=strsplit(tellurics[k],'hip',/extract)
               'hd':tellurics_end=strsplit(tellurics[k],'hd',/extract)
               'hr':tellurics_end=strsplit(tellurics[k],'hr',/extract)
               else: begin
                  print, 'Cannot match telluric name given by user '+$
                         'to object name in any fits file headers.'
                  stop
               end
            endcase
            test=strmatch(object_tmp,'*'+tellurics_end+'*',/fold_case)
         endif

         ;if the telluric files have been located and are science
         ;exposures (not acquisition exposures) then create the
         ;subdirectory and file lists. note, can't cut on exposure
         ;time for telluric stars, so rely on acquisition exposures
         ;for tellurics to be correctly labeled as acqCal and not
         ;partnerCal.
         index_telfiles= where(test GT 0 AND obsclass_tmp EQ $
                               'partnerCal', count_telfiles)
         if count_telfiles GT 0 then begin

            if (file_test(workdir+'tellurics/'+dates[i]+'/'+setup+$
                          '/'+tellurics[k],/dir) EQ 0) then $
                             spawn, 'mkdir '+workdir+'tellurics/'+$
                                    dates[i]+'/'+setup+'/'+$
                                    tellurics[k]

            ;could have the same telluric observed multiple times
            ;during the same night in the same setup. don't want to
            ;combine all of these exposures together so use the ut
            ;time to distinguish multiple sets of observations.
            wentthru = 'no'
            counter = 1
            telgroup = intarr(n_elements(test))
            for kk = 0, n_elements(test)-1 do begin
               if test[kk] EQ 1 AND obsclass_tmp[kk] EQ $
                  'partnerCal' then begin
                  ut_tel_split = strsplit(ut_tmp[kk],':',/extract)
                  ut_tel = ut_tel_split[0] + ut_tel_split[1]/60.D
                  ;pull out the first telluric file of the first set
                  if wentthru EQ 'no' then ut_tel_ref = ut_tel
                  ;assume that a telluric belongs to a set if the
                  ;observation was taken within 30 mins of the first
                  ;exposure.
                  if ut_tel - ut_tel_ref LE (30.D/60.D) then begin
                     telgroup[kk] = counter
                  endif else begin
                     counter = counter + 1
                     telgroup[kk] = counter
                     ut_tel_ref = ut_tel
                  endelse
                  wentthru = 'yes'
               endif
            endfor
            
            for kk = 1, max(telgroup) do begin

               openw, lun1, workdir+'tellurics/'+dates[i]+'/'+setup+$
                      '/'+tellurics[k]+'/telluriclist_'+$
                      strcompress(string(kk),/remove_all), /get_lun
               openw, lun2, workdir+'tellurics/'+dates[i]+'/'+setup+$
                      '/'+tellurics[k]+'/skylist_'+$
                      strcompress(string(kk),/remove_all), /get_lun
               openw, lun3, workdir+'tellurics/'+dates[i]+'/'+setup+$
                      '/'+tellurics[k]+'/skylistshort_'+$
                      strcompress(string(kk),/remove_all), /get_lun

               for kkk = 0, n_elements(test)-1 do begin
                  if test[kkk] EQ 1 AND obsclass_tmp[kkk] EQ $
                     'partnerCal' AND telgroup[kkk] EQ kk then begin
                     spawn, 'cp '+filestem_tmp[kkk]+'.fits '+workdir+$
                            'tellurics/'+dates[i]+'/'+setup+'/'+$
                            tellurics[k]+'/'
                     ;to determine if this is an object or sky
                     ;exposure of the telluric, look at the x and y
                     ;offsets. must be within a radius of 2" to be an
                     ;object exposure.
                     radius_tmp = Sqrt(xoff_tmp[kkk]^2 + yoff_tmp[kkk]^2)
                     if radius_tmp LE 2.D then $
                        printf, lun1, filestem_tmp[kkk], format='(A)'
                     if radius_tmp GE 2.D then $
                        printf, lun3, filestem_tmp[kkk], format='(A)'
                  endif
               endfor
               close, lun1
               free_lun, lun1
               close, lun3
               free_lun, lun3
               
               ;determine which sky observations go with which
               ;telluric observations. separate the case of just 1 sky
               ;exposure, otherwise use +1/-1 of the file number
               readcol, workdir+'tellurics/'+dates[i]+'/'+setup+'/'+$
                        tellurics[k]+'/telluriclist_'+$
                        strcompress(string(kk),/remove_all), $
                        telluric_file_tmp, format='A', /silent
               readcol, workdir+'tellurics/'+dates[i]+'/'+setup+'/'+$
                        tellurics[k]+'/skylistshort_'+$
                        strcompress(string(kk),/remove_all), $
                        telluricskyshort_file_tmp, format='A', /silent
               if n_elements(telluricskyshort_file_tmp) EQ 1 then begin
                  for kkk=0, n_elements(telluric_file_tmp)-1 do begin
                     printf, lun2, telluricskyshort_file_tmp[0], $
                             format='(A)'
                  endfor
               endif else begin
                  for kkk=0, n_elements(telluric_file_tmp)-1 do begin
                     index_telluricsky = $
                        where(strmid(telluricskyshort_file_tmp,3,$
                                     /reverse_offset) EQ $
                              strmid(telluric_file_tmp[kkk],3,$
                                     /reverse_offset) - 1 OR $
                              strmid(telluricskyshort_file_tmp,3,$
                                     /reverse_offset) EQ $
                              strmid(telluric_file_tmp[kkk],3,$
                                     /reverse_offset) + 1, count_telluricsky)
                     if count_telluricsky EQ 1 then begin
                        printf, lun2, $
                                telluricskyshort_file_tmp[index_telluricsky], $
                                format='(A)'
                     ;if cannot find a sky with a file number +1/-1
                     ;within the telluric file number, just use the
                     ;first sky. will remind the user to check the
                     ;telluric and sky lists.
                     endif else begin
                        printf, lun2, telluricskyshort_file_tmp[0], $
                                format='(A)'
                     endelse
                  endfor
               endelse
               close, lun2
               free_lun, lun2

            endfor
               
         endif
      endfor
      

      ;for the objects
      for k = 0, n_elements(galaxies)-1 do begin

         ;see if the full galaxy name listed by user matches the
         ;object name in the data file.
         test = strmatch(object_tmp, galaxies[k], /fold_case)
         index_tmp = where(test GT 0, count_tmp)
         ;if the above does not work, see if the end of the galaxy
         ;name listed by the user matches some part of the object name
         ;in the data file.
         if count_tmp EQ 0 then begin
            catalog_firsttwoletters = strmid(galaxies[k],0,2)
            case catalog_firsttwoletters of
               'ng':galaxies_end=strsplit(galaxies[k],'ngc',/extract)
               'ug':galaxies_end=strsplit(galaxies[k],'ugc',/extract)
               'mr':galaxies_end=strsplit(galaxies[k],'mrk',/extract)
               'pg':galaxies_end=strsplit(galaxies[k],'pgc',/extract)
               'ic':galaxies_end=strsplit(galaxies[k],'ic',/extract)
               else: begin
                  print, 'Cannot match galaxy name given by user '+$
                         'to object name in any fits file headers.'
                  stop
               end
            endcase
            test=strmatch(object_tmp,'*'+galaxies_end+'*',/fold_case)
         endif

         ;if the galaxy files have been located and are science
         ;exposures (not acquisition exposures) then create the
         ;subdirectory and file lists.
         index_galfiles= where(test GT 0 AND obsclass_tmp EQ $
                               'science' AND exptime_tmp GE 180., $
                               count_galfiles)
         if count_galfiles GT 0 then begin
      
            if (file_test(workdir+galaxies[k]+'/'+dates[i]+'/'+$
                          setup,/dir) EQ 0) then $
                             spawn, 'mkdir '+workdir+galaxies[k]+'/'+$
                                    dates[i]+'/'+setup
            if (file_test(workdir+galaxies[k]+'/merged/'+$
                          setup,/dir) EQ 0) then $
                             spawn, 'mkdir '+workdir+galaxies[k]+$
                                    '/merged/'+setup
            
            ;note, creating one list for all galaxy exposures. the
            ;user will have to modify later depending on which
            ;exposures should be in a single group and reduced with a
            ;single telluric
            openw, lun1, workdir+galaxies[k]+'/'+dates[i]+'/'+setup+$
                   '/gallist', /get_lun
            openw, lun2, workdir+galaxies[k]+'/'+dates[i]+'/'+setup+$
                   '/skylist', /get_lun
            openw, lun3,workdir+galaxies[k]+'/'+dates[i]+'/'+setup+$
                   '/skylistshort', /get_lun

            for kk = 0, n_elements(test)-1 do begin
               if test[kk] EQ 1 AND obsclass_tmp[kk] EQ $
                  'science' AND exptime_tmp[kk] GE 180. then begin
                  spawn, 'cp '+filestem_tmp[kk]+'.fits '+workdir+$
                         galaxies[k]+'/'+dates[i]+'/'+setup+'/'
                  ;to determine if this is an object or sky exposure
                  ;of the telluric, look at the x and y offsets. must
                  ;be within a radius of 2" to be an object exposure.
                  radius_tmp = Sqrt(xoff_tmp[kk]^2 + yoff_tmp[kk]^2)
                  if radius_tmp LE 2.D then $
                     printf, lun1, filestem_tmp[kk], format='(A)'
                  if radius_tmp GE 2.D then $
                     printf, lun3, filestem_tmp[kk], format='(A)'
               endif
            endfor
            close, lun1
            free_lun, lun1
            close, lun3
            free_lun, lun3

            ;determine which sky observations go with which galaxy
            ;observations. separate the case of just 1 sky exposure,
            ;otherwise use +1/-1 of the file number
            readcol, workdir+galaxies[k]+'/'+dates[i]+'/'+setup+$
                     '/gallist', gal_file_tmp, format='A', /silent
            readcol, workdir+galaxies[k]+'/'+dates[i]+'/'+setup+$
                     '/skylistshort', galskyshort_file_tmp, $
                     format='A', /silent
            if n_elements(galskyshort_file_tmp) EQ 1 then begin
               for kk=0, n_elements(gal_file_tmp)-1 do begin
                  printf, lun2, galskyshort_file_tmp[0], $
                          format='(A)'
               endfor
            endif else begin
               for kk=0, n_elements(gal_file_tmp)-1 do begin
                  index_galsky = where(strmid(galskyshort_file_tmp,3,$
                                              /reverse_offset) EQ $
                                       strmid(gal_file_tmp[kk],3,$
                                              /reverse_offset) - 1 OR $
                                       strmid(galskyshort_file_tmp,3,$
                                              /reverse_offset) EQ $
                                       strmid(gal_file_tmp[kk],3,$
                                              /reverse_offset) + 1, $
                                       count_galsky)
                  if count_galsky EQ 1 then begin
                     printf, lun2, galskyshort_file_tmp[index_galsky], $
                             format='(A)'
                  ;if cannot find a sky with a file number +1/-1
                  ;within the galaxy file number, just use the first
                  ;sky. will remind the user to check the galaxy and
                  ;sky lists.
                  endif else begin
                     printf, lun2, galskyshort_file_tmp[0],$
                             format='(A)'
                  endelse
               endfor
            endelse
            close, lun2
            free_lun, lun2

         endif
      endfor


      ;for the psf. this is setup to deal with one PSF star taken
      ;before the galaxy exposures.

      ;create directories to hold psf star observations on this date
      ;in this setup. will delete later if no observation exists.
      if (file_test(workdir+'psfs/'+dates[i]+'/'+setup,/dir) EQ 0) then $
         spawn, 'mkdir '+workdir+'psfs/'+dates[i]+'/'+setup
      if (file_test(workdir+'psfs/merged/'+setup,/dir) EQ 0) then $
         spawn, 'mkdir '+workdir+'psfs/merged/'+setup

      ;find the psf stars for this night. assume that the object name
      ;contains 'star' (e.g., SFOTuningStar, PSFStar, * - NGS Star) or
      ;contains 'tt' (e.g., Off-axis TTGS)
      index1_tmp = strmatch(object_tmp,'*star*',/fold_case)
      index2_tmp = strmatch(object_tmp,'*tt*',/fold_case)
      index_psf = index1_tmp + index2_tmp
      index_tmp = where(index_psf EQ 2, count_tmp)
      if count_tmp GT 0 then index_psf[index_tmp] = 1

      ;only continue if there was an observation of a PSF star on this
      ;night.
      index_tmp = where(index_psf EQ 1, count_tmp)
      if count_tmp GT 0 then begin

         ;assume the psf star observation was taken within ~45 mins
         ;before the first galaxy exposure it corresponds to.
         for k = 0, n_elements(galaxies)-1 do begin
            readcol, workdir+galaxies[k]+'/'+dates[i]+'/'+setup+$
                     '/gallist', file_tmp, format='A', /silent

            index_gal1 = where(filestem_tmp EQ file_tmp[0])
            ut_galexp1_split = strsplit(ut_tmp[index_gal1],':',$
                                        /extract)
            ut_galexp1 = ut_galexp1_split[0] + $
                         ut_galexp1_split[1]/60.D

            ;go through once to know if there is a PSF star
            ;observation for this galaxy
            wentthru = 'no'
            for kk = 0, n_elements(index_psf)-1 do begin
               if index_psf[kk] EQ 1 then begin
                  ut_psf_split = strsplit(ut_tmp[kk],':',/extract)
                  ut_psf = ut_psf_split[0] + ut_psf_split[1]/60.D
                  if ut_psf GE (ut_galexp1 - (45.D/60.D)) AND $
                     ut_psf LE ut_galexp1 then begin
                     wentthru = 'yes'
                  endif
               endif
            endfor

            ;if there is a PSF star on this night for this particular
            ;galaxy, make directory, copy files, and make the object
            ;and sky lists.
            if wentthru EQ 'yes' then begin

               if (file_test(workdir+'psfs/'+dates[i]+'/'+setup+'/psf_'+$
                             galaxies[k],/dir) EQ 0) then $
                                spawn, 'mkdir '+workdir+'psfs/'+dates[i]+$
                                       '/'+setup+'/psf_'+galaxies[k]
               if (file_test(workdir+'psfs/merged/'+setup+'/psf_'+$
                             galaxies[k],/dir) EQ 0) then $
                                spawn, 'mkdir '+workdir+'psfs/merged/'+$
                                       setup+'/psf_'+galaxies[k]
               
               openw, lun1, workdir+'psfs/'+dates[i]+'/'+setup+$
                      '/psf_'+galaxies[k]+'/psflist', /get_lun
               openw, lun2, workdir+'psfs/'+dates[i]+'/'+setup+$
                      '/psf_'+galaxies[k]+'/skylist', /get_lun
               openw, lun3, workdir+'psfs/'+dates[i]+'/'+setup+$
                      '/psf_'+galaxies[k]+'/skylistshort', /get_lun

               for kk = 0, n_elements(index_psf)-1 do begin
                  if index_psf[kk] EQ 1 then begin
                     ut_psf_split = strsplit(ut_tmp[kk],':',/extract)
                     ut_psf = ut_psf_split[0] + ut_psf_split[1]/60.D
                     if ut_psf GE (ut_galexp1 - (45.D/60.D)) AND $
                        ut_psf LE ut_galexp1 then begin
                        spawn, 'cp '+filestem_tmp[kk]+'.fits '+$
                               workdir+'psfs/'+dates[i]+'/'+setup+$
                               '/psf_'+galaxies[k]+'/'
                        ;to determine if this is an object or sky
                        ;exposure of the psf, look at the x and y
                        ;offsets. must be within a radius of 2" to be
                        ;an object exposure.
                        radius_tmp = Sqrt(xoff_tmp[kk]^2 + $
                                          yoff_tmp[kk]^2)
                        if radius_tmp LE 2.D then $
                           printf, lun1, filestem_tmp[kk],format='(A)'
                        if radius_tmp GE 2.D then $
                           printf, lun3, filestem_tmp[kk],format='(A)'
                     endif
                  endif
               endfor
               close, lun1
               free_lun, lun1
               close, lun3
               free_lun, lun3
               
               ;determine which sky observations go with which psf
               ;observations. separate the case of just 1 sky
               ;exposure, otherwise use +1/-1 of the file number
               readcol, workdir+'psfs/'+dates[i]+'/'+setup+'/psf_'+$
                        galaxies[k]+'/psflist', psf_file_tmp, $
                        format='A', /silent
               readcol, workdir+'psfs/'+dates[i]+'/'+setup+'/psf_'+$
                        galaxies[k]+'/skylistshort', $
                        psfskyshort_file_tmp, format='A', /silent
               if n_elements(psfskyshort_file_tmp) EQ 1 then begin
                  for kk=0, n_elements(psf_file_tmp)-1 do begin
                     printf, lun2, psfskyshort_file_tmp[0], $
                             format='(A)'
                  endfor
               endif else begin
                  for kk=0, n_elements(psf_file_tmp)-1 do begin
                     index_psfsky = where(strmid(psfskyshort_file_tmp,3,$
                                                 /reverse_offset) EQ $
                                          strmid(psf_file_tmp[kk],3,$
                                                 /reverse_offset) - 1 OR $
                                          strmid(psfskyshort_file_tmp,3,$
                                                 /reverse_offset) EQ $
                                          strmid(psf_file_tmp[kk],3,$
                                                 /reverse_offset) + 1, $
                                          count_psfsky)
                     if count_psfsky EQ 1 then begin
                        printf, lun2, psfskyshort_file_tmp[index_psfsky], $
                                format='(A)'
                     ;if cannot find a sky with a file number +1/-1
                     ;within the psf file number, just use the first
                     ;sky. will remind the user to check the psf and
                     ;sky lists.
                     endif else begin
                        printf, lun2, psfskyshort_file_tmp[0], $
                                format='(A)'
                     endelse
                  endfor
               endelse
               close, lun2
               free_lun, lun2

            endif

         endfor

      ;otherwise delete the setup directory if there is no PSF star
      ;observation this night in this setup for this galaxy
      endif else begin
         spawn, 'rm -r '+workdir+'psfs/'+dates[i]+'/'+setup
         spawn, 'rm -r '+workdir+'psfs/merged/'+setup
      endelse


      ;sort calibration data files and make lists (e.g., flatlist,
      ;flatdarklist)
      ;-------------------------------------------------------------

      ;for the object darks
      darks_index = where(obstype_tmp EQ 'DARK' AND $
                          aperture_tmp EQ 'Blocked' AND $
                          exptime_tmp GE 180., ndarks)
      if ndarks EQ 0 then begin
         print, 'Cannot find object darks.'
         stop
      endif
      ;make lists of object darks that are all of the same exposure
      ;time. for example, could have 300s darks for an object for this
      ;night, but also 600s for a different object for this same
      ;night.
      exptime_darks = exptime_tmp[darks_index]
      uniq_darks = $
         exptime_darks[uniq(exptime_darks(sort(exptime_darks)))]
      for k=0, n_elements(uniq_darks)-1 do begin
         
         openw, lun, workdir+'daycals/'+dates[i]+'/'+setup+$
                '/darklist_'+strcompress(string(uniq_darks[k]),$
                                         /remove_all), /get_lun
         for kk=0, ndarks-1 do begin
            if exptime_darks[kk] EQ uniq_darks[k] then begin
               spawn, 'cp '+filestem_tmp[darks_index[kk]]+'.fits '+$
                      workdir+'daycals/'+dates[i]+'/'+setup+'/'
               printf, lun, filestem_tmp[darks_index[kk]],format='(A)'
            endif 
         endfor
         
         close, lun
         free_lun, lun
      endfor


      ;for the flats
      flats_index = where(obstype_tmp EQ 'FLAT' AND $
                          aperture_tmp EQ '3.0' AND $
                          medianim_tmp GT 200., nflats)
      if nflats EQ 0 then begin
         print, 'Cannot find flats.'
         stop
      endif
      openw, lun, workdir+'daycals/'+dates[i]+'/'+setup+'/flatlist',$
             /get_lun
      for k=0, nflats-1 do begin
         spawn, 'cp '+filestem_tmp[flats_index[k]]+'.fits '+workdir+$
                'daycals/'+dates[i]+'/'+setup+'/'
         printf, lun, filestem_tmp[flats_index[k]], format='(A)'
      endfor
      close, lun
      free_lun, lun

      
      ;for the flat darks
      flatdarks_index = where(obstype_tmp EQ 'FLAT' AND $
                              aperture_tmp EQ '3.0' AND $
                              medianim_tmp LT 200., nflatdarks)
      if nflatdarks EQ 0 then begin
         print, 'Cannot find flat darks.'
         stop
      endif
      openw, lun, workdir+'daycals/'+dates[i]+'/'+setup+$
             '/flatdarklist', /get_lun
      for k=0, nflatdarks-1 do begin
         spawn, 'cp '+filestem_tmp[flatdarks_index[k]]+'.fits '+$
                workdir+'daycals/'+dates[i]+'/'+setup+'/'
         printf, lun, filestem_tmp[flatdarks_index[k]], format='(A)'
      endfor
      close, lun
      free_lun, lun

      
      ;for the ronchi mask
      ronchi_index = where(obstype_tmp EQ 'FLAT' AND $
                           aperture_tmp EQ 'Ronchi' AND $
                           medianim_tmp GT 200.,nronchi)
      if nronchi EQ 0 then begin
         print, 'Cannot find ronchi mask.'
         stop
      endif
      openw, lun, workdir+'daycals/'+dates[i]+'/'+setup+$
             '/ronchilist', /get_lun
      for k=0, nronchi-1 do begin
         spawn, 'cp '+filestem_tmp[ronchi_index[k]]+'.fits '+workdir+$
                'daycals/'+dates[i]+'/'+setup+'/'
         printf, lun, filestem_tmp[ronchi_index[k]], format='(A)'
      endfor
      close, lun
      free_lun, lun

      ;for the ronchi mask darks
      ronchidarks_index = where(obstype_tmp EQ 'FLAT' AND $
                                aperture_tmp EQ 'Ronchi' AND $
                                medianim_tmp LT 20.,nronchidarks)
      if nronchidarks EQ 0 then begin
         print, 'Cannot find ronchi mask darks.'
         stop
      endif
      openw, lun, workdir+'daycals/'+dates[i]+'/'+setup+$
             '/ronchidarklist', /get_lun
      for k=0, nronchidarks-1 do begin
         spawn, 'cp '+filestem_tmp[ronchidarks_index[k]]+'.fits '+workdir+$
                'daycals/'+dates[i]+'/'+setup+'/'
         printf, lun, filestem_tmp[ronchidarks_index[k]], format='(A)'
      endfor
      close, lun
      free_lun, lun
      

      ;for the arcs.
      arcs_index = where(obstype_tmp EQ 'ARC', narcs)
      if narcs EQ 0 then begin
         print, 'Cannot find arcs.'
         stop
      endif

      ;any arcs taken back-to-back are put into the same file list,
      ;otherwise each arc gets its own file list.
      counter = 1
      for k=0, narcs-1 do begin
         spawn, 'cp '+filestem_tmp[arcs_index[k]]+'.fits '+workdir+$
                'daycals/'+dates[i]+'/'+setup+'/'
         if k EQ 0 then begin
            openw, lun, workdir+'daycals/'+dates[i]+'/'+setup+$
                   '/arclist_'+strcompress(string(counter),$
                                           /remove_all),/get_lun
            printf, lun, filestem_tmp[arcs_index[k]], format='(A)'
         endif else begin
            if strmid(filestem_tmp[arcs_index[k]],3,/reverse_offset) - 1 EQ $
               strmid(filestem_tmp[arcs_index[k-1]],3,$
                      /reverse_offset) then begin
               printf, lun, filestem_tmp[arcs_index[k]], format='(A)'
            endif else begin
               close, lun
               free_lun, lun
               counter = counter + 1
               openw, lun, workdir+'daycals/'+dates[i]+'/'+setup+$
                      '/arclist_'+strcompress(string(counter),$
                                              /remove_all),/get_lun
               printf, lun, filestem_tmp[arcs_index[k]], format='(A)'
            endelse
         endelse
      endfor
      close, lun
      free_lun, lun

      
      ;for the arc darks
      arcdarks_index = where(obstype_tmp EQ 'DARK' AND $
                             aperture_tmp EQ 'Blocked' AND $
                             exptime_tmp LT 60., narcdarks)
      if narcdarks EQ 0 then begin
         print, 'Cannot find arc darks.'
      endif
      openw, lun, workdir+'daycals/'+dates[i]+'/'+setup+$
             '/arcdarklist', /get_lun
      for k=0, narcdarks-1 do begin
         spawn, 'cp '+filestem_tmp[arcdarks_index[k]]+'.fits '+$
                workdir+'daycals/'+dates[i]+'/'+setup+'/'
         printf, lun, filestem_tmp[arcdarks_index[k]], format='(A)'
      endfor
      close, lun
      free_lun, lun

   endfor
   
   ;remove directory if there is no psf star observation on this night
   if file_test(workdir+'psfs/'+dates[i]+'/*') EQ 0 then $
      spawn, 'rm -r '+workdir+'psfs/'+dates[i]
   if file_test(workdir+'psfs/merged/*') EQ 0 then $
      spawn, 'rm -r '+workdir+'psfs/merged'
   
endfor

;remove directory if there was no psf star observation for any of the
;nights being reduced.
if file_test(workdir+'psfs/*') EQ 0 then spawn, 'rm -r '+workdir+'psfs'

end
