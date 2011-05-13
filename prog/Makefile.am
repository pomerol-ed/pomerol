MAINTAINERCLEANFILES = Makefile.in
AM_CXXFLAGS+= ${hdf5_CFLAGS} ${eigen3_CFLAGS} -I../jsoncpp/include/

bin_PROGRAMS=pomerolDiag

pomerolDiag_CPPFLAGS=-I../lib
pomerolDiag_SOURCES=main.cpp 
pomerolDiag_LDADD=../lib/libpomerol.a ../jsoncpp/src/libjsoncpp.a ${hdf5_LIBS} -lhdf5 -lhdf5_cpp

pomerolDiag_DEPENDENCIES = $(top_builddir)/jsoncpp/src/libjsoncpp.a $(top_builddir)/lib/libpomerol.a

all-clean:
	$(RM) *.o pomerolDiag
