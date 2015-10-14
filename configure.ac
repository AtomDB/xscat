#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([xscat], [1.0.0], [rsmith@cfa.harvard.edu])
AC_CONFIG_SRCDIR([src/xscat.c])
AC_CONFIG_HEADERS([config.h])
#AC_CONFIG_MACRO_DIRS([m4])

# Checks for programs.
AC_PROG_CC
AC_PROG_F77([gfortran])
AC_F77_LIBRARY_LDFLAGS
AC_PROG_INSTALL
AC_PROG_RANLIB

#Init automake
AM_INIT_AUTOMAKE
#AM_SILENT_RULES([yes])

# Determine name of OS.
ARCH_FLAGS=""
UNAME_SYSTEM=`(uname -s) 2>/dev/null` || UNAME_SYSTEM=unknown
if test "${UNAME_SYSTEM}" = "Darwin"; then
   ARCH_FLAGS="-arch x86_64"
fi
AC_SUBST([ARCH_FLAGS])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h string.h strings.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt strchr])

AC_CONFIG_FILES([Makefile
	         extlib/Makefile
		 extlib/ape/Makefile
		 extlib/ape/src/Makefile
		 extlib/ape/include/Makefile
		 extlib/ape/include/ape/Makefile
                 extlib/mie/Makefile
		 xspec/xscat.h
                 src/Makefile
		 src/xscat-setup.csh
		 src/xscat-setup.sh])

AC_CONFIG_SUBDIRS([extlib/cfitsio])

AC_OUTPUT
