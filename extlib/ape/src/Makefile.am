AM_CFLAGS=-I@top_srcdir@/extlib/ape/include

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LIBRARIES=libape.a

# Sources:
libape_a_SOURCES=ape_trad.c ape_util.c ape_test.c ape_session.c	\
                 ape_par_group.c ape_par.c ape_msg.c ape_list.c	\
                 ape_error.c pil_error.c pil.c ape_io.c ape_binary.c

install-exec-hook:
	mkdir -p $(top_builddir)/lib
	cp -p libape.a $(top_builddir)/lib
	mkdir -p $(top_builddir)/include
	cp -p ../include/*.h $(top_builddir)/include
	mkdir -p $(top_builddir)/include/ape
	cp -p ../include/ape/*.h $(top_builddir)/include/ape
