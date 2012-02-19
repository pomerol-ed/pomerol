MAINTAINERCLEANFILES = Makefile.in
AM_CXXFLAGS+= ${hdf5_CFLAGS} ${eigen3_CFLAGS} -I../jsoncpp/include/ -I../src
#AM_LDFLAGS=../src/libpomerol.a ../jsoncpp/src/libjsoncpp.a ${hdf5_LIBS} -lhdf5 -lhdf5_cpp

check_PROGRAMS=gamma4 green
TESTS=gamma4 green 

#gamma4
gamma4_SOURCES=gamma4.cpp
green_SOURCES=green.cpp

gamma4_DEPENDENCIES = $(top_builddir)/jsoncpp/src/libjsoncpp.a $(top_builddir)/src/libpomerol.a
green_DEPENDENCIES = $(top_builddir)/jsoncpp/src/libjsoncpp.a $(top_builddir)/src/libpomerol.a

all-clean:
	$(RM) *.o $(TESTS_programs)
