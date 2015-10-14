.run read_rawxscat
.run sigmascat

mrn = read_rawxscat('MRN')
energy = read_energy('MRN')
xpos = [0,0.2,0.5,0.75,0.9,0.95,0.990,0.999] ; 8 values
rad  = [10., 30., 60., 120.]  ; 4 values

!p.font=1
!p.charsize=1.3
!xtitle='Energy (keV)'
!ytitle='Sigma (cm!E2!N)'
plot_io, energy, mrn[0,0,*]

oplot, energy, mrn[0,2,*]
oplot, energy, mrn[1,2,*],linestyle=1
oplot, energy, mrn[2,2,*],linestyle=2
oplot, energy, mrn[3,2,*],linestyle=3
oplot, energy, mrn[0,4,*],linestyle=1,thick=3
;oplot, energy, mrn[0,7,*]

oplot, energy, sigmascat(energy),thick=3

