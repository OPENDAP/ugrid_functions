
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

noinst_LTLIBRARIES = libugrid_functions.la

libugrid_functions_la_SOURCES = $(SRCS) $(HDRS)
libugrid_functions_la_CPPFLAGS = $(GF_CFLAGS) $(XML2_CFLAGS)
libugrid_functions_la_LDFLAGS = -avoid-version -module 
libugrid_functions_la_LIBADD = $(GF_LIBS) 

SRCS = \
	gf2.cc \
	ugrid_utils.cc \
	MeshDataVariable.cc \
	TwoDMeshTopology.cc  \
	gf3.cc

# NC_GOES_Dataset.cpp NC_GOES_Dataset.h

HDRS = \
	LocationType.h \
	ugrid_utils.h \
	MeshDataVariable.h  \
	TwoDMeshTopology.h \
	gf2.h \
	gf3.h
