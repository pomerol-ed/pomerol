MAINTAINERCLEANFILES = Makefile.in
AM_CXXFLAGS+= ${hdf5_CFLAGS} ${eigen3_CFLAGS} -I../jsoncpp/include/ -I../src
LDADD=../src/libpomerol.a ../jsoncpp/src/libjsoncpp.a ${hdf5_LIBS} -lhdf5 -lhdf5_cpp

check_PROGRAMS=IndexTerm IndexPermutation_test CCdagOperatorTest NOperatorTest SzOperatorTest HamiltonianPartTest HamiltonianTest FieldOperatorPartTest \
FieldOperatorTest
TESTS=${check_PROGRAMS}

#gamma4_SOURCES=gamma4.cpp
#green_SOURCES=green.cpp
IndexTerm_SOURCES=IndexTerm.cpp
IndexPermutation_test_SOURCES=IndexPermutation_test.cpp
CCdagOperatorTest_SOURCES=CCdagOperatorTest.cpp
NOperatorTest_SOURCES=NOperatorTest.cpp
SzOperatorTest_SOURCES=SzOperatorTest.cpp
HamiltonianPartTest_SOURCES=HamiltonianPartTest.cpp
HamiltonianTest_SOURCES=HamiltonianTest.cpp
FieldOperatorPartTest_SOURCES=FieldOperatorPartTest.cpp
FieldOperatorTest_SOURCES=FieldOperatorTest.cpp

#gamma4_DEPENDENCIES = $(top_builddir)/jsoncpp/src/libjsoncpp.a $(top_builddir)/src/libpomerol.a
#green_DEPENDENCIES = $(top_builddir)/jsoncpp/src/libjsoncpp.a $(top_builddir)/src/libpomerol.a
TESTS_DEPENDENCIES = $(top_builddir)/jsoncpp/src/libjsoncpp.a $(top_builddir)/src/libpomerol.a

all-clean:
	$(RM) *.o $(check_programs)
