setenv XSCAT ../..
setenv PFILES "/tmp;../../syspfiles"

foreach x (0.200 0.500 0.750 0.900 0.950 0.990 0.999)
   foreach dm ( MRN )
      ../../src/xscat mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
	Emin=1.0 deltaE=0.1 NumberOfEnergies=1 \
	ExtractRadius=0.0  DustPosition=${x} \
	DustModelName=${dm} DustModel=-1 OutputFileName=xs_${dm}_${x}_000 \
	>& log_xs_${dm}_${x}_000.txt &
   end 
end

