XSCAT_TOP = ..
AM_CFLAGS = -Wall -I$(top_builddir)/src -I$(XSCAT_TOP)/extlib/cfitsio -I$(XSCAT_TOP)/extlib/ape/include $(ARCH_FLAGS)
AM_FFLAGS = $(ARCH_FLAGS)

bin_PROGRAMS = xscat

bin_SCRIPTS = xscat-setup.csh xscat-setup.sh

xscat_SOURCES = xscat.c get_input_data.c messages.c calc_dustmodel.c \
	wrap_send_sigmaint.c expand_env.c quadrature_dustsize.c \
	load_dustmodel.c read_ZDA_fits.c qtrap_dm.c trapzd_dm.c \
	wrap_receive_sigmaint.f get_opt_constant.f \
	mie_func.f calc_theta_sca.f mie_interp.f interp_mie_func.f \
	interp_function.f adapt.f qtrap.f r_err.f sort3.f indexx.f \
	messages.h xscat.h xscat.inc 

xscat_LDADD = ../extlib/cfitsio/libcfitsio.a ../extlib/ape/src/libape.a \
	../extlib/mie/libmie.a $(FLIBS)

install-exec-hook:
	cp -p xscat $(top_builddir)/bin/.
	chmod a+x xscat-setup.csh
	chmod a+x xscat-setup.sh
	cp -p xscat-setup.csh $(top_builddir)/bin/.
	cp -p xscat-setup.sh $(top_builddir)/bin/.
	mkdir -p $(top_builddir)/syspfiles
	cp -p $(top_builddir)/pfiles/xscat.par $(top_builddir)/syspfiles/.
