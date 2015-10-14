read_array,'xs_WD3100AGAL_0.000_010.dat',a00
read_array,'xs_WD3100AGAL_0.000_030.dat',a01
read_array,'xs_WD3100AGAL_0.000_060.dat',a02
read_array,'xs_WD3100AGAL_0.000_120.dat',a03
read_array,'xs_WD3100AGAL_0.200_010.dat',a10
read_array,'xs_WD3100AGAL_0.200_030.dat',a11
read_array,'xs_WD3100AGAL_0.200_060.dat',a12
read_array,'xs_WD3100AGAL_0.200_120.dat',a13
read_array,'xs_WD3100AGAL_0.500_010.dat',a20
read_array,'xs_WD3100AGAL_0.500_030.dat',a21
read_array,'xs_WD3100AGAL_0.500_060.dat',a22
read_array,'xs_WD3100AGAL_0.500_120.dat',a23
read_array,'xs_WD3100AGAL_0.750_010.dat',a30
read_array,'xs_WD3100AGAL_0.750_030.dat',a31
read_array,'xs_WD3100AGAL_0.750_060.dat',a32
read_array,'xs_WD3100AGAL_0.750_120.dat',a33
read_array,'xs_WD3100AGAL_0.900_010.dat',a40
read_array,'xs_WD3100AGAL_0.900_030.dat',a41
read_array,'xs_WD3100AGAL_0.900_060.dat',a42
read_array,'xs_WD3100AGAL_0.900_120.dat',a43
read_array,'xs_WD3100AGAL_0.950_010.dat',a50
read_array,'xs_WD3100AGAL_0.950_030.dat',a51
read_array,'xs_WD3100AGAL_0.950_060.dat',a52
read_array,'xs_WD3100AGAL_0.950_120.dat',a53
read_array,'xs_WD3100AGAL_0.990_010.dat',a60
read_array,'xs_WD3100AGAL_0.990_030.dat',a61
read_array,'xs_WD3100AGAL_0.990_060.dat',a62
read_array,'xs_WD3100AGAL_0.990_120.dat',a63

!p.multi=[0,2,2]
plot_io,a00[0,*],a00[1,*],yrange=[1e-25,1e-21]
oplot,a00[0,*],a00[2,*]

plot_io,a01[0,*],a01[1,*],yrange=[1e-25,1e-21]
oplot,a00[0,*],a01[2,*]

plot_io,a02[0,*],a02[1,*],yrange=[1e-25,1e-21]
oplot,a00[0,*],a02[2,*]

plot_io,a03[0,*],a03[1,*],yrange=[1e-25,1e-21]
oplot,a00[0,*],a03[2,*]


!p.multi=[0,2,2]
plot_io,a10[0,*],a10[1,*],yrange=[1e-25,1e-21]
oplot,a10[0,*],a10[2,*]

plot_io,a11[0,*],a11[1,*],yrange=[1e-25,1e-21]
oplot,a10[0,*],a11[2,*]

plot_io,a12[0,*],a12[1,*],yrange=[1e-25,1e-21]
oplot,a10[0,*],a12[2,*]

plot_io,a13[0,*],a13[1,*],yrange=[1e-25,1e-21]
oplot,a10[0,*],a13[2,*]


!p.multi=[0,2,2]
plot_io,a20[0,*],a20[1,*]
oplot,a20[0,*],a20[2,*]

plot_io,a21[0,*],a21[1,*]
oplot,a20[0,*],a21[2,*]

plot_io,a22[0,*],a22[1,*]
oplot,a20[0,*],a22[2,*]

plot_io,a23[0,*],a23[1,*]
oplot,a20[0,*],a23[2,*]

!p.multi=[0,2,2]
plot_io,a30[0,*],a30[1,*]
oplot,a30[0,*],a30[2,*]

plot_io,a31[0,*],a31[1,*]
oplot,a30[0,*],a31[2,*]

plot_io,a32[0,*],a32[1,*]
oplot,a30[0,*],a32[2,*]

plot_io,a33[0,*],a33[1,*]
oplot,a30[0,*],a33[2,*]

!p.multi=[0,2,2]
plot_io,a40[0,*],a40[1,*]
oplot,a40[0,*],a40[2,*]

plot_io,a41[0,*],a41[1,*]
oplot,a40[0,*],a41[2,*]

plot_io,a42[0,*],a42[1,*]
oplot,a40[0,*],a42[2,*]

plot_io,a43[0,*],a43[1,*]
oplot,a40[0,*],a43[2,*]

!p.multi=[0,2,2]
plot_io,a50[0,*],a50[1,*]
oplot,a50[0,*],a50[2,*]

plot_io,a51[0,*],a51[1,*]
oplot,a50[0,*],a51[2,*]

plot_io,a52[0,*],a52[1,*]
oplot,a50[0,*],a52[2,*]

plot_io,a53[0,*],a53[1,*]
oplot,a50[0,*],a53[2,*]

!p.multi=[0,2,2]
plot_io,a60[0,*],a60[1,*]
oplot,a60[0,*],a60[2,*]

plot_io,a61[0,*],a61[1,*]
oplot,a60[0,*],a61[2,*]

plot_io,a62[0,*],a62[1,*]
oplot,a60[0,*],a62[2,*]

plot_io,a63[0,*],a63[1,*]
oplot,a60[0,*],a63[2,*]


