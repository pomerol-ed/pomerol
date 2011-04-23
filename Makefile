EIGEN_INCLUDE=$(shell pkg-config --cflags-only-I eigen3) -I$(shell pwd)/eigen
INCLUDES= -I../jsoncpp/include $(EIGEN_INCLUDE) $(shell pkg-config --cflags-only-I hdf5)
POMEROL_LIBS = $(shell pkg-config --libs hdf5)
CFLAGS=$(INCLUDES) -DHRD -DpomerolHDF5
CXXFLAGS=$(CFLAGS)

export CFLAGS
export CXXFLAGS
export POMEROL_LIBS

all : jsoncpp libpomerol

lib:
	@mkdir -p lib

.PHONY: jsoncpp
jsoncpp: lib
	@echo "Bulding jsoncpp"
	cd jsoncpp && $(MAKE)
	@mv jsoncpp/libs/lib_json.a ./lib
	@echo

.PHONY: libpomerol
libpomerol:	lib
	@echo "Building libpomerol"
	cd src && $(MAKE)
	cp src/libpomerol.a ./lib
	@echo

clean:
	cd src && $(MAKE) clean
	cd jsoncpp && $(MAKE) clean

.PHONY: doc-html
doc-html:
	doxygen doc/Doxyfile
