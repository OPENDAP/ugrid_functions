
# Build libugridFunctions, part of libdap.

AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I conf

# These are not used by automake but are often useful for certain types of
# debugging. The best way to use these is to run configure as:
#     export CXXFLAGS='...'; ./configure --disable-shared
# the --disable-shared is not required, but it seems to help with debuggers.
CXXFLAGS_DEBUG = -g3 -O0  -Wall -Wcast-align
TEST_COV_FLAGS = -ftest-coverage -fprofile-arcs

SUBDIRS =  . unit-tests tests
DIST_SUBDIRS = unit-tests tests

lib_besdir=$(libdir)/bes
lib_bes_LTLIBRARIES = libugrid_functions.la

SRCS = UgridFunctions.cc \
	ugrid_utils.cc \
	MeshDataVariable.cc \
	TwoDMeshTopology.cc  \
	ugr5.cc  \
	NDimensionalArray.cc 

HDRS = UgridFunctions.h\
	LocationType.h \
	ugrid_utils.h \
	MeshDataVariable.h  \
	TwoDMeshTopology.h \
	ugr5.h \
	NDimensionalArray.h 

libugrid_functions_la_SOURCES = $(SRCS) $(HDRS)
libugrid_functions_la_CPPFLAGS = $(GF_CFLAGS) $(XML2_CFLAGS)
libugrid_functions_la_LDFLAGS = -avoid-version -module 
libugrid_functions_la_LIBADD = $(GF_LIBS) 

EXTRA_DIST = data COPYING ugrid_functions.conf.in 

package_version=$(PACKAGE_MAJOR_VERSION).$(PACKAGE_MINOR_VERSION).$(PACKAGE_SUBMINOR_VERSION)

CLEANFILES = *~ ugrid_functions.conf ugrid_functions-$(package_version).tar.gz

# ugrid_functions.spec OSX_Resources

install-data-local: ugrid_functions.conf
	@currdate=`date +"%y%m%d%H%M"`; \
	test -d $(DESTDIR)$(sysconfdir)/bes || $(MKDIR_P) $(DESTDIR)$(sysconfdir)/bes; \
	test -d $(DESTDIR)$(sysconfdir)/bes/modules || $(MKDIR_P) $(DESTDIR)$(sysconfdir)/bes/modules; \
	test -f $(DESTDIR)$(sysconfdir)/bes/modules/ugrid_functions.conf && mv -f $(DESTDIR)$(sysconfdir)/bes/modules/ugrid_functions.conf $(DESTDIR)$(sysconfdir)/bes/modules/ugrid_functions.conf.$$currdate; \
	$(INSTALL_DATA) ugrid_functions.conf $(DESTDIR)$(sysconfdir)/bes/modules/ugrid_functions.conf

uninstall-local:
	rm -f $(DESTDIR)$(sysconfdir)/bes/modules/ugrid_functions.conf


ugrid_functions.conf: ugrid_functions.conf.in $(top_builddir)/config.status
	sed -e "s%[@]pkgdatadir[@]%${pkgdatadir}%" \
                -e "s%[@]sysconfdir[@]%${sysconfdir}%" \
                -e "s%[@]pkglibdir[@]%${pkglibdir}%" \
                -e "s%[@]datarootdir[@]%${datarootdir}%" \
                -e "s%[@]datadir[@]%${datadir}%" \
                -e "s%[@]bindir[@]%${bindir}%" \
                -e "s%[@]bes_modules_dir[@]%${lib_besdir}%" $< > ugrid_functions.conf
