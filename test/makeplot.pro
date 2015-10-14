.run getcolor
colors=getcolor(/load)

read_array,'MRNx000.dat',m000
read_array,'MRNx001.dat',m001
read_array,'MRNx002.dat',m002
read_array,'MRNx500.dat',m500
read_array,'MRNx750.dat',m750
read_array,'MRNx950.dat',m950
read_array,'MRNx990.dat',m990
read_array,'MRNx999.dat',m999
;read_array,'WD3100AGALx999.dat',wd999
;read_array,'ZDACAGx999.dat',z999


read_array,'Test1-sherpa.dat',t1s
read_array,'Test1.dat',t1
read_array,'Test2.dat',t2
read_array,'Test3.dat',t3
read_array,'Test4.dat',t4
read_array,'TestWD.dat',tWD
read_array,'TestZDACAF.dat',tZ


plot_io,m000[0,*],m000[3,*],yrange=[1e-28,1e-22],thick=3
oplot,m500[0,*],m500[3,*],color=colors.red
oplot,m750[0,*],m750[3,*],color=colors.blue

oplot,m950[0,*],m950[3,*],color=colors.yellow
oplot,m990[0,*],m990[3,*],color=colors.pink
oplot,m999[0,*],m999[3,*],color=colors.green


