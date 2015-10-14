.run read_rawxscat
.run sigmascat

read_array,'xs_MRN_0.200_000.dat',a200
read_array,'xs_MRN_0.500_000.dat',a500
read_array,'xs_MRN_0.750_000.dat',a750
read_array,'xs_MRN_0.900_000.dat',a900
read_array,'xs_MRN_0.950_000.dat',a950
read_array,'xs_MRN_0.990_000.dat',a990
read_array,'xs_MRN_0.999_000.dat',a999

;mrn = read_rawxscat('MRN')
;energy = read_energy('MRN')
;xpos = [0,0.2,0.5,0.75,0.9,0.95,0.990,0.999] ; 8 values
;rad  = [10., 30., 60., 120.]  ; 4 values

plot_io,a200[0,*],a200[2,*]
oplot,  a500[0,*],a500[2,*]
oplot,  a750[0,*],a750[2,*]
oplot,  a900[0,*],a900[2,*],linestyle=1
oplot,  a950[0,*],a950[2,*]
oplot,  a990[0,*],a990[2,*]
oplot,  a999[0,*],a990[2,*]
;oplot, energy, mrn[0,4,*]
;oplot, energy, mrn[0,5,*]
;oplot, energy, mrn[0,6,*]
;oplot, energy, mrn[0,7,*]

;oplot, energy, 0.5*sigmascat(energy),psym=1

