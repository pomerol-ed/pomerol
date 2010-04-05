default: src

.PHONY : src
src : 
	mkdir -p libs
	cd src/lib_json && $(MAKE) 
	mv src/lib_json/lib_json.a ./libs

clean:
	cd src/lib_json && $(MAKE) clean
