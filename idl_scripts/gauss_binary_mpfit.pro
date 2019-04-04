function gauss_binary_mpfit,x,p

; this routine is called by measure_lsf_lp in order to fit a double
; gaussian to OH doublet lines. the separation between the two
; gaussians is known (gauss_sep), and the amplitude and dispersion of
; the gaussians are the same. a constant is added to the double
; gaussian function.
  
COMMON skyblock, gauss_sep
  
n = n_elements(p)
nx = N_ELEMENTS(x)

z1 = (x-p[1]-gauss_sep/2.)/p[2]
z2 = (x-p[1]+gauss_sep/2.)/p[2]
ez1 = exp(-z1^2/2.)
ez2 = exp(-z2^2/2.)

f = p[0]*(ez1+ez2)/2. + p[3]

return, f

end
