read_array,'mie_dat.out',mie, LENGTH=20000

E_cm = 2.0
a = reform(mie[3,*])
Th_cm = reform(mie[0,*])
c = 1.1
fdst_henke = 1.0
rho = 3.3

Isca = c * a^6 * (rho/3.0)^2 * exp(-0.4575d0*a*a*E_cm*E_cm*(Th_cm/60.0)^2)

plot_oo,Th_cm,mie[4,*],psym=1
oplot,Th_cm,Isca,psym=2

