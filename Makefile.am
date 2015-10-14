AUTOMAKE_OPTIONS = foreign
SUBDIRS = extlib src 
EXTRA_DIST = xspec inputs production doc test pfiles CHANGELOG

AM_CFLAGS = $(ARCH_FLAGS)
AM_FFLAGS = $(ARCH_FLAGS)

clean-local:
	-rm -rf bin
	-rm -rf lib
	-rm -rf include
	-rm -rf syspfiles
