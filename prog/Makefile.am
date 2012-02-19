MAINTAINERCLEANFILES = Makefile.in
AM_CXXFLAGS+= ${hdf5_CFLAGS} ${eigen3_CFLAGS} -I../jsoncpp/include/

bin_PROGRAMS=pomerolDiag 

# pomerolDiag
pomerolDiag_CPPFLAGS=-I../src
pomerolDiag_SOURCES=DevelopersExample.cpp 
pomerolDiag_LDADD=../src/libpomerol.a ../jsoncpp/src/libjsoncpp.a ${hdf5_LIBS} -lhdf5 -lhdf5_cpp

pomerolDiag_DEPENDENCIES = $(top_builddir)/jsoncpp/src/libjsoncpp.a $(top_builddir)/src/libpomerol.a

all-clean:
	$(RM) *.o $(bin_PROGRAMS)
