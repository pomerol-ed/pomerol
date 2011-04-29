JSON_INCLUDE=-I../include/ 
default : lib_json.a

json_objects=json_reader.o json_writer.o json_value.o

json_reader.o: json_reader.cpp
	$(CXX) $(JSON_INCLUDE) -c json_reader.cpp

json_writer.o: json_writer.cpp
	$(CXX) $(JSON_INCLUDE) -c json_writer.cpp

json_value.o: json_value.cpp
	$(CXX) $(JSON_INCLUDE) -c json_value.cpp

lib_json.a: $(json_objects)
	ar -cvr lib_json.a $(json_objects) 

clean: 
	rm -fv *.o *.a 
	
