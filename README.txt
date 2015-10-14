XSCAT 1.0.0 
Author: Randall Smith
Email: rsmith@cfa.harvard.edu


** Introduction **

The purpose of `xscat' is to calculate the X-ray scattering cross section of
a population of dust grains as a function of energy, given a specific
dust grain position and an extraction region.  As the grain position
gets closer to the source of the X-rays for a constant extraction
region, the cross section drops as more X-rays remain within the
extraction circle.  For constant dust grain position, the cross
section decreases with increasing extraction region as well, again
because more photons remain within it.  

`xscat' refers to two separate programs.  The main code is written in
a combination of C and Fortran.  It uses the astronomical standard
`IRAF' user interface to set all the basic parameters and then outputs
a text file with the requested values.  The second code is a
multiplicative XSPEC model
(https://heasarc.gsfc.nasa.gov/xanadu/xspec/) that is designed to be
used in combination with an absorption model such as `tbabs' or
`phabs.'  Most users will only wish to use the XSPEC model; the main
code is only required to create new dust models for the XSPEC model, or to
test new dust models.

** The Main Code **

Compilation: The code uses the standard GNU interface, and has been
tested on both Linux and Mac systems using gcc and gfortran.  To
install, use the following commands:

unix% tar zxf xscat-1.0.0.tar.gz
unix% cd xscat-1.0.0
unix% ./configure --prefix=`pwd`
unix% make
unix% make install

If you need help with the configure step, try

unix% ./configure --help

Running: To run the main code, you need to first set the XSCAT
environmental variable to the installed directory of xscat, as well as
the PFILES variable to indicate where the parameter files are installed.  In
general, this will be the same directory as it was compiled into.
Assuming you have just followed the above instructions to compile
`xscat', you need only do:

unix% setenv XSCAT `pwd`
unix% setenv PFILES "/tmp;$XSCAT/syspfiles"

Once this is complete, you can see all the parameters used by `xscat'
using the command:

unix% plist xscat
(OutputFileName = TestZDACAF)      Output Event file stem
(DustModelName = ZDACAF)          Dust Model Name (see documentation)
   (DustModel = -1)              Dust Model Number (see documentation)
        (Emin = 1.0)             Minimum Energy of X-ray to consider (keV)
      (deltaE = 0.5)             Delta Energy step of next X-ray to consider (keV)
(NumberOfEnergies = 3)               Number of Energies to calculate
(ExtractRadius = 10.0)            Source Extraction Radius (arcsec)
(DustPosition = 0.5)             "Relative Dust Position (obs=0.0; src=1.0)
 (Interpolate = yes)             Use faster interpolation when calculating scattering
     (Epsilon = 1e-3)            Accuracy to achieve in numerical integration
       (Drude = no)              Use Drude approximation to optical constants (not recommended)
     (clobber = no)              Overwrite existing output file?
        (mode = hl)              mode for parameter file

These are largely self-documenting, except for the DustModelName and
DustModel values.  The DustModel number is a lookup into the dust
model array, and not recommended as it's easy to mix up.  Setting this
to -1 means it is ignored.  Easier to use is the DustModelName.
The existing coded dust models include the MRN77, WD01, and
ZDA04-based models.  See the Smith, Valencic & Corrales (2015) paper
in this directory for details.  The allowed names are:

MRN77: "MRN"
WD01:  "WD3100AGAL","WD3110AGAL","WD3120AGAL","WD3130AGAL","WD3140AGAL",
       "WD3150AGAL","WD3160AGAL","WD4000AGAL","WD4010AGAL","WD4020AGAL",
       "WD4030AGAL","WD4040AGAL","WD5500AGAL","WD5510AGAL","WD5520AGAL",
       "WD5530AGAL","WD4000BGAL","WD4010BGAL","WD4020BGAL","WD4030BGAL",
       "WD4040BGAL","WD5500BGAL","WD5510BGAL","WD5520BGAL","WD5530BGAL",
       "WD2600ALMC","WD2610ALMC","WD2620ALMC","WD2600BLMC","WD2605BLMC",
       "WD2610BLMC","WD2900ASMC"

ZDA04: "ZDABGS","ZDABGF","ZDABGB","ZDACGS","ZDACGF",
       "ZDACGB","ZDABAS","ZDABAF","ZDABAB","ZDACAS",
       "ZDACAF","ZDACAB","ZDACNS","ZDACNF","ZDACNB"

The code comes by default with the MRN, WD3100AGAL, and ZDABAS models
calculated.  Not all have been run simply because they require
approximately a week per, and some require fine-tuning to get the
optical constants to work properly.  

The WD01 model names are WDXXYYZSRC, where XX is the R_V value, YY is
the b_C value in the WD01 paper, Z is A or B depending upon the model
type, and SRC is the source of the model, either Galactic, LMC, or
SMC.

A sample run using the default values used for the XSPEC models is:

unix% xscat mode=hl clobber=no Epsilon=1e-3 Drude=no Interpolate=yes \
      Emin=0.1 deltaE=0.002 NumberOfEnergies=1450 \
      ExtractRadius=10.0  DustPosition=0.5 \
      DustModelName=MRN DustModel=-1 OutputFileName=xs_MRN_0.500_010 >& \
          log_xs_MRN_0.500_010.txt

** The XSPEC Module **

Compilation: Compiling the XSPEC model requires at least configuring the main
code, and having an installed version of FTOOLS/XSPEC that has been
built from source.  Note that a binary installation is not adequate.
A sample run is shown here.

unix% ...install FTOOLS using whatever your command is...
unix% tar zxf xscat-1.0.0.tar.gz
unix% cd xscat-1.0.0
unix% ./configure --prefix=`pwd`
unix% cd xspec
unix% xspec
XSPEC> initpackage xscat model.dat .

Using xscat in XSPEC: After compiling the model code, it must be
installed into XSPEC each time XSPEC is run.  This is done using the
`lmod` command: 

XSPEC> lmod xscat /path/to/xscat/xspec

where /path/to/xscat is the path to the directory the XSCAT xspec
model is installed in.  Note that the xspec model is in the xspec
directory, so this will be $XSCAT/xspec in general.  One installed,
the model is simply called `xscat' and can be used as a multiplicative
model, usually in combination with an absorption model such as:

XSPEC> model tbabs*xscat*pow

Model TBabs<1>*xscat<2>*powerlaw<3> Source No.: 1   Active/Off
Model Model Component  Parameter  Unit     Value
 par  comp
   1    1   TBabs      nH         10^22    1.00000      +/-  0.0          
   2    2   xscat      NH         10^22    1.00000      +/-  0.0          
   3    2   xscat      Xpos                0.500000     +/-  0.0          
   4    2   xscat      Rext       arcsec   10.0000      frozen
   5    2   xscat      DustModel           1            frozen
   6    3   powerlaw   PhoIndex            1.00000      +/-  0.0          
   7    3   powerlaw   norm                1.00000      +/-  0.0          

Note that xscat has four parameters, two of which should always be
frozen.  The NH value is the interstellar hydrogen column density --
NOT necessarily the total column density, which could include
absorption near the source that does not contribute to scattering
losses.  See Smith, Valencic & Corrales (2015) for all of the
details. 
