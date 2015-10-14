FUNCTION sigmascat, energy

  ; from WD01, referencing DL84
  Cgra = 10^(-25.13) ; units of cm^2.5
  Cslc = 10^(-25.11) ; units of cm^2.5
  ; from MG84, with (0.1 um)^-4 scaled to 1e20 cm^-4
  norm = 6.3e-11 * 1e20
  ; default densities
  rho_gra = 2.2 ; in g/cm^3
  rho_slc = 3.3 ; in g/cm^3
  amax = 0.25e-4 ; in cm 
  amin = 5e-8    ; in cm 
  scaling =  norm * (Cgra + Cslc) * ((rho_gra/3.0)^2 + (rho_slc/3.0)^2) * $
	(1/1.5)*(amax^1.5 - amin^1.5)
  result = scaling/(energy^2)

  return, result
end
