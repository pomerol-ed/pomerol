EIGEN_INCLUDE=-I/opt/local/include/eigen3 $(shell pkg-config --cflags-only-I eigen3) -I$(shell pwd)/eigen
INCLUDES= -I../iniparser/src -I../iniconfig -I../jsoncpp/include $(EIGEN_INCLUDE)
CFLAGS=$(INCLUDES) -DHRD # -DpomerolHDF5
CXXFLAGS=$(CFLAGS)

export CFLAGS
export CXXFLAGS

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

devel:
	cd src && $(MAKE) devel

.PHONY: doc-html
doc-html:
	doxygen doc/Doxyfile
