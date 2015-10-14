#setenv XSCATCODE ../../src/xscat
setenv XSCAT ../..
#setenv PFILES "/tmp;../../syspfiles"

foreach x (0.200 0.500 0.750 0.900 0.950 0.990 0.999)
   foreach dm ( WD3100AGAL )
      ../../src/xscat mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=0.1 deltaE=0.1 NumberOfEnergies=49 \
	ExtractRadius=60.0  DustPosition=${x} \
	DustModelName=${dm} DustModel=-1 OutputFileName=xs_${dm}_${x}_060 \
	>& log_xs_${dm}_${x}_060.txt &

      ../../src/xscat mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=0.1 deltaE=0.1 NumberOfEnergies=49 \
	ExtractRadius=120.0  DustPosition=${x} \
	DustModelName=${dm} DustModel=-1 OutputFileName=xs_${dm}_${x}_120 \
	>& log_xs_${dm}_${x}_120.txt &
   end 
end

