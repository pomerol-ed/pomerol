default: src

.PHONY : src
src : 
	mkdir -p libs
	cd src && $(MAKE) 
	mv src/lib_json.a ./libs

clean:
	cd src && $(MAKE) clean
