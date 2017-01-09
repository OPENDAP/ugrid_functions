# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Configure script for OPeNDAP Gateway data handler

AC_PREREQ(2.59)
AC_INIT([ugrid_functions], [1.0.8], [opendap-tech@opendap.org])
AC_CONFIG_AUX_DIR([conf])
AC_CONFIG_MACRO_DIR([conf])

AM_INIT_AUTOMAKE(1.10)

AM_CONDITIONAL([DAP_MODULES], [false])

# The handler is a shared object library, so we don't normally need
# the static lib. Override this default by specifying `--enable-static',
# which can make running unit tests with valgrind far more pleasant. 
AM_DISABLE_STATIC

AC_CONFIG_TESTDIR(tests, [.])

DODS_GCOV_VALGRIND

PACKAGE_MAJOR_VERSION=`echo $PACKAGE_VERSION | sed 's@^\([[0-9]]\)*\.\([[0-9]]*\)\.\([[0-9]]*\)$@\1@'`
PACKAGE_MINOR_VERSION=`echo $PACKAGE_VERSION | sed 's@^\([[0-9]]\)*\.\([[0-9]]*\)\.\([[0-9]]*\)$@\2@'`
PACKAGE_SUBMINOR_VERSION=`echo $PACKAGE_VERSION | sed 's@^\([[0-9]]\)*\.\([[0-9]]*\)\.\([[0-9]]*\)$@\3@'`
AC_MSG_NOTICE(Package Major Version:     $PACKAGE_MAJOR_VERSION)
AC_MSG_NOTICE(Package Minor Version:     $PACKAGE_MINOR_VERSION)
AC_MSG_NOTICE(Package SubMinor Version:  $PACKAGE_SUBMINOR_VERSION)
AC_SUBST(PACKAGE_MAJOR_VERSION)
AC_SUBST(PACKAGE_MINOR_VERSION)
AC_SUBST(PACKAGE_SUBMINOR_VERSION)

AC_CONFIG_SRCDIR([ugrid_utils.h])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_CHECK_FUNCS([atexit strchr])

dnl Checks for specific libraries
AC_CHECK_LIBDAP([3.13.0], 
	[ LIBS="$LIBS $DAP_LIBS"  CPPFLAGS="$CPPFLAGS $DAP_CFLAGS"],
	[ AC_MSG_ERROR([Cannot find libdap]) ])

AC_CHECK_BES([3.13.0], 
	[LIBS="$LIBS $BES_LIBS"  CPPFLAGS="$CPPFLAGS $BES_CPPFLAGS"], 
	[ AC_MSG_ERROR([Could not find bes library and headers]) ])


# Is there CPPUNIT?
AM_PATH_CPPUNIT(1.12.0,
	[AM_CONDITIONAL([CPPUNIT], [true])],
	[AM_CONDITIONAL([CPPUNIT], [false])])


dnl When pkgconfig support is added, use this
libdap_libgridfields_module='libgridfields >= 1.0.3'
libgf_needed="1.0.3"

libgf_ok="no"

AC_ARG_WITH([gridfields],
    AC_HELP_STRING([--with-gridfields=path], [Use the gridfields library at this location.]),
    [if test "$withval" = "yes"
     then
         AC_MSG_ERROR([You must supply a path when using --with-gridfields; try --enable-gridfields.])
     else
         AC_MSG_NOTICE([Using $withval as the GridFields prefix directory.])
         GF_LIBS="-L$withval/lib -lgridfields"
         GF_CFLAGS="-I$withval/include"
	 GF_PATH=$withval
         AC_DEFINE([GRIDFIELDS],[1],[define if gridfields lib is present])
         libgf_ok="yes"
     fi],
    [])

AC_ARG_ENABLE([gridfields],
    AC_HELP_STRING([--enable-gridfields], [Enable gridfields (default is YES)]),
    [],
    [enable_gridfields=yes])

if test "$libgf_ok" = "no" -a "$enable_gridfields" = "yes"
then
     
AC_MSG_CHECKING([for libgridfields])
if gridfields-config --version > /dev/null 2>&1
then
    version_libgf=`gridfields-config --version | sed 's@gridfields \(.*\)@\1@'`
    AX_COMPARE_VERSION([$version_libgf], [ge], [$libgf_needed], [libgf_ok='yes'], [libgf_ok='no'])

    if test "$libgf_ok"="no"
    then
        AC_MSG_ERROR([must have libgf $libgf_needed or greater, found $version_libgf])
    fi

    GF_LIBS="`gridfields-config --libs`"
    GF_CFLAGS="`gridfields-config --cflags`"
    GF_PATH="`gridfields-config --prefix`"

    AC_DEFINE([GRDIFIELDS],[1],[define if gridfields lib is present])
    AC_MSG_RESULT([yes; used gridfields-config])
else
    AC_MSG_ERROR([I could not find libgridfields])
fi

fi

AC_SUBST([GF_LIBS])
AC_SUBST([GF_CFLAGS])
AC_SUBST([GF_PATH])

if test "$libgf_ok" = "yes"
then
    echo "Set gridfields CFLAGS to $GF_CFLAGS"
    echo "Set gridfields LIBS to $GF_LIBS"
    AM_CONDITIONAL([USING_GRIDFIELDS],[true])
else
    AC_MSG_ERROR([You must provide the location of the GridFields library. Use --with-gridfields to do so...])
fi

dnl OPENDAP_DEBUG_OPTION

## This sets up the special symbols used to build the bes.conf file for
## the tests. Previously, this was done in the Makefile.am down in the
## test directory, but building the bes.conf here works better with
## autoconf overall and esp. with the 'modules' build. jhrg 11/19/14
modulesdir=$BES_MODULE_DIR
AC_SUBST(modulesdir)

AC_CONFIG_FILES([Makefile 
	unit-tests/Makefile
	unit-tests/test_config.h
	tests/Makefile 
	tests/package.m4 
	tests/atlocal
	tests/bes.conf])

AC_CONFIG_FILES([tests/generate_data_baseline.sh], 
	[chmod 755 tests/generate_data_baseline.sh])
AC_CONFIG_FILES([tests/generate_metadata_baseline.sh],
	[chmod 755 tests/generate_metadata_baseline.sh])

AC_OUTPUT
