XSCATCODE=../src/xscat
export XSCAT=..
export PFILES=/tmp;../syspfiles

test: Test1.dat TestWD.dat TestZDACAF.dat Test2.dat Test3.dat Test4.dat

testx: MRNx000.dat MRNx001.dat MRNx002.dat MRNx950.dat MRNx990.dat \
	MRNx999.dat MRNx500.dat MRNx750.dat WD3100AGALx999.dat ZDACAFx999.dat

oldtest: Test1a.dat Test1b.dat

clean:
	rm *.dat
#
# Basic tests 
#
Mie.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=2.0 deltaE=0.2 NumberOfEnergies=1 \
	ExtractRadius=0.0  DustPosition=0.5 \
	DustModelName=MRN DustModel=-1 OutputFileName=Mie >& mie_dat.out 

Test1.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.2 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.5 \
	DustModelName=MRN DustModel=-1 OutputFileName=Test1 

Test2.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.2 NumberOfEnergies=20 \
	ExtractRadius=0.0  DustPosition=0.5 \
	DustModelName=MRN DustModel=-1 OutputFileName=Test2 

Test3.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.2 NumberOfEnergies=20 \
	ExtractRadius=120.0  DustPosition=0.5 \
	DustModelName=MRN DustModel=-1 OutputFileName=Test3 

Test4.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=0.2 deltaE=0.2 NumberOfEnergies=50 \
	ExtractRadius=10.0  DustPosition=0.5 \
	DustModelName=MRN DustModel=-1 OutputFileName=Test4 

Test1a.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.2 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.5 \
	DustModelName=MRN DustModel=0 OutputFileName=Test1a

Test1b.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.2 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.5 \
	DustModelName=MRN DustModel=1 OutputFileName=Test1b

TestWD.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.5 NumberOfEnergies=3 \
	ExtractRadius=10.0  DustPosition=0.5 \
	DustModelName=WD3100AGAL DustModel=-1 OutputFileName=TestWD

TestZDACAF.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.5 NumberOfEnergies=3 \
	ExtractRadius=10.0  DustPosition=0.5 \
	DustModelName=ZDACAF DustModel=-1 OutputFileName=TestZDACAF

MRNx000.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.25 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.0 \
	DustModelName=MRN DustModel=-1 OutputFileName=MRNx000

MRNx001.dat: 
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.25 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.001 \
	DustModelName=MRN DustModel=-1 OutputFileName=MRNx001

MRNx002.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.25 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.002 \
	DustModelName=MRN DustModel=-1 OutputFileName=MRNx002

MRNx500.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.25 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.500 \
	DustModelName=MRN DustModel=-1 OutputFileName=MRNx500

MRNx750.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.25 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.75 \
	DustModelName=MRN DustModel=-1 OutputFileName=MRNx750

MRNx950.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.25 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.95 \
	DustModelName=MRN DustModel=-1 OutputFileName=MRNx950

MRNx990.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.25 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.99 \
	DustModelName=MRN DustModel=-1 OutputFileName=MRNx990

MRNx999.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.25 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.999 \
	DustModelName=MRN DustModel=-1 OutputFileName=MRNx999

WD3100AGALx999.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.25 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.999 \
	DustModelName=WD3100AGAL DustModel=-1 OutputFileName=WD3100AGALx999

ZDACAFx999.dat:
	$(XSCATCODE) mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.25 NumberOfEnergies=20 \
	ExtractRadius=10.0  DustPosition=0.999 \
	DustModelName=ZDACAF DustModel=-1 OutputFileName=ZDACAFx999
