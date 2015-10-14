.run read_rawxscat
.run sigmascat
.run getcolor

colors=getcolor(/load)

mrn = read_rawxscat('ZDABAS')
energy = read_energy('ZDABAS')
xpos = [0,0.2,0.5,0.75,0.9,0.95,0.990,0.999] ; 8 values
rad  = [10., 30., 60., 120.]  ; 4 values

!p.font=1
!p.charsize=1.3
!xtitle='Energy (keV)'
!ytitle='Sigma (cm!E2!N)'
plot_io, energy, mrn[0,0,*],yrange=[1e-24,1e-20]

oplot, energy, mrn[0,2,*]
oplot, energy, mrn[1,2,*],color=colors.red
oplot, energy, mrn[2,2,*],color=colors.blue
oplot, energy, mrn[3,2,*],color=colors.green
oplot, energy, mrn[0,4,*],color=colors.yellow
;oplot, energy, mrn[0,7,*]

oplot, energy, sigmascat(energy),thick=3

oplot, energy, ionabs(energy,/iKeV)

