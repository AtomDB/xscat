lib_LIBRARIES = libmie.a

AM_CFLAGS = $(ARCH_FLAGS)
AM_FFLAGS = $(ARCH_FLAGS)

libmie_a_SOURCES = MIEV0.f ErrPack.f DMiLay.f

install-exec-hook:
	mkdir -p $(top_builddir)/lib
	cp -p libmie.a $(top_builddir)/lib
