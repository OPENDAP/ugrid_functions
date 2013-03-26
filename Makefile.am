
# Build libugridFunctions, part of libdap.

AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I conf

# AM_CPPFLAGS = 

# These are not used by automake but are often useful for certain types of
# debugging. The best way to use these is to run configure as:
#     export CXXFLAGS='...'; ./configure --disable-shared
# the --disable-shared is not required, but it seems to help with debuggers.
CXXFLAGS_DEBUG = -g3 -O0 -fno-defer-pop -Wall -W -Wcast-align -Werror
TEST_COV_FLAGS = -ftest-coverage -fprofile-arcs



# SUBDIRS = 
# DIST_SUBDIRS = 

lib_besdir=$(libdir)/bes
lib_bes_LTLIBRARIES = libugrid_functions.la

#noinst_LTLIBRARIES = libugrid_functions.la

SRCS = UgridFunctions.cc \
	ugrid_utils.cc \
	MeshDataVariable.cc \
	TwoDMeshTopology.cc  \
	ugrid_restrict.cc  \
	gf3.cc

#        gf2.cc 

# NC_GOES_Dataset.cpp NC_GOES_Dataset.h

HDRS = UgridFunctions.h\
	LocationType.h \
	ugrid_utils.h \
	MeshDataVariable.h  \
	TwoDMeshTopology.h \
	ugrid_restrict.h \
	gf3.h

#       gf2.h 

libugrid_functions_la_SOURCES = $(SRCS) $(HDRS)
libugrid_functions_la_CPPFLAGS = $(GF_CFLAGS) $(XML2_CFLAGS)
libugrid_functions_la_LDFLAGS = -avoid-version -module 
libugrid_functions_la_LIBADD = $(GF_LIBS) 


install-data-local: ugrid_functions.conf
	@currdate=`date +"%y%m%d%H%M"`; \
	test -d $(DESTDIR)$(sysconfdir)/bes || $(MKDIR_P) $(DESTDIR)$(sysconfdir)/bes; \
	test -d $(DESTDIR)$(sysconfdir)/bes/modules || $(MKDIR_P) $(DESTDIR)$(sysconfdir)/bes/modules; \
	test -f $(DESTDIR)$(sysconfdir)/bes/modules/ugrid_functions.conf && mv -f $(DESTDIR)$(sysconfdir)/bes/modules/ugrid_functions.conf $(DESTDIR)$(sysconfdir)/bes/modules/ugrid_functions.conf.$$currdate; \
	$(INSTALL_DATA) ugrid_functions.conf $(DESTDIR)$(sysconfdir)/bes/modules/ugrid_functions.conf

uninstall-local:
	rm -f $(DESTDIR)$(sysconfdir)/bes/modules/ugrid_functions.conf

ugrid_functions.conf: ugrid_functions.conf.in ../config.status
	sed -e "s%[@]pkgdatadir[@]%${pkgdatadir}%" \
                -e "s%[@]sysconfdir[@]%${sysconfdir}%" \
                -e "s%[@]pkglibdir[@]%${pkglibdir}%" \
                -e "s%[@]datarootdir[@]%${datarootdir}%" \
                -e "s%[@]datadir[@]%${datadir}%" \
                -e "s%[@]bindir[@]%${bindir}%" \
                -e "s%[@]bes_modules_dir[@]%${lib_besdir}%" $< > ugrid_functions.conf
