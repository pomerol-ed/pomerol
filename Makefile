all : iniconfig jsoncpp meat

.PHONY: iniconfig
iniconfig:
	@echo "Creating iniconfig"
	@mkdir -p lib
	cd iniconfig && $(MAKE)
	@mv iniconfig/libiniconfig.a ./lib
	@echo

.PHONY: jsoncpp
jsoncpp:
	@echo "Creating jsoncpp"
	@mkdir -p lib
	cd jsoncpp && $(MAKE)
	@mv jsoncpp/libs/lib_json.a ./lib
	@echo

.PHONY: meat
meat:
	@echo "Creating libmeat"
	@mkdir -p lib
	cd src && $(MAKE)
	cp src/libmeat.a ./lib
	@echo

clean:
	cd src && $(MAKE) clean
	cd iniconfig && $(MAKE) clean
	cd jsoncpp && $(MAKE) clean

