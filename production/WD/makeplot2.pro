.run read_rawxscat
.run sigmascat
.run getcolor

colors=getcolor(/load)

mrn = read_rawxscat('WD3100AGAL')
energy = read_energy('WD3100AGAL')
xpos = [0,0.2,0.5,0.75,0.9] ; 8 values
rad  = [10., 30., 60., 120.]  ; 4 values

!p.font=1
!p.charsize=1.3
!xtitle='Energy (keV)'
!ytitle='Sigma (cm!E2!N)'
plot_io, energy, mrn[0,0,*],yrange=[1e-25,1e-21]

oplot, energy, mrn[0,0,*]
oplot, energy, mrn[0,1,*],color=colors.red
oplot, energy, mrn[0,2,*],color=colors.blue
oplot, energy, mrn[0,3,*],color=colors.green
oplot, energy, mrn[0,4,*],color=colors.yellow
;oplot, energy, mrn[0,5,*]
;oplot, energy, mrn[0,6,*]
;oplot, energy, mrn[0,7,*]

;oplot, energy, sigmascat(energy),thick=3



