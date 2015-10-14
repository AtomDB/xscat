# The sub-directories are built before the current directory.
# In order to change this, include "." in the list of SUBDIRS.
SUBDIRS=. ape
DIST_SUBDIRS=ape
EXTRA_DIST=cfitsio mie

# We assume that the configure-script is executed by the main configure.ac

all-local:
	make all-cfitsio
	make all-mie

install-data-local: 
	make install-cfitsio
	make install-mie

clean-local:
	make clean-cfitsio
	make clean-mie

##### CFITSIO #####

cfitsio_dir=cfitsio
all-cfitsio:
	cd ${cfitsio_dir} &&  make shared stand_alone

# need to create and linke an adjusted Makefile
install-cfitsio:
	cd ${cfitsio_dir} && make install

clean-cfitsio:
	cd ${cfitsio_dir} && make clean

###### MIE #######
mie_dir=mie
all-mie:
	cd ${mie_dir} &&  make

install-mie:
	cd ${mie_dir} && make install

clean-mie:
	cd ${mie_dir} && make clean


###### APE  #######
ape_dir=ape
all-ape:
	cd ${ape_dir} &&  make 

install-ape:
	cd ${ape_dir} && make install

clean-ape:
	cd ${ape_dir} && make clean

