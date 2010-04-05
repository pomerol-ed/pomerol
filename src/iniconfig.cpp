#include<cstring>
#include"iniconfig.h"

#include<sstream>

// IniConfig::exCouldNotOpenIniFile class
IniConfig::exCouldNotOpenIniFile::exCouldNotOpenIniFile(std::string Name) throw() : 
	FileName(Name), FormatString("Could not open '%s'.")
{
	Message = new char[strlen(FormatString) - 2 + strlen(FileName.c_str()) + 1];
}
IniConfig::exCouldNotOpenIniFile::~exCouldNotOpenIniFile() throw() {delete[] Message;}
const char* IniConfig::exCouldNotOpenIniFile::what() const throw()
{
	::sprintf(Message,FormatString,FileName.c_str());
	return Message;
}

// IniConfig::exKeyNotFound class
IniConfig::exKeyNotFound::exKeyNotFound(std::string K) throw() : Key(K), FormatString("Key '%s' not found.") 
{
	Message = new char[strlen(FormatString) - 2 + strlen(Key.c_str()) + 1];
}

IniConfig::exKeyNotFound::~exKeyNotFound() throw() {delete[] Message;}
const char* IniConfig::exKeyNotFound::what() const throw()
{
	sprintf(Message,FormatString,Key.c_str());
	return Message;
}

// IniConfig::exValueTypeMismatch class
IniConfig::exValueTypeMismatch::exValueTypeMismatch(std::string V, std::string Requested) throw() :
	RequestedType(Requested), Value(V), FormatString("Unable to cast the value '%s' to type '%s'.")
{
	Message = new char[	strlen(FormatString) - 4 + 
				strlen(Value.c_str()) + 
				strlen(RequestedType.c_str()) +1];
}
IniConfig::exValueTypeMismatch::~exValueTypeMismatch() throw() {delete[] Message;}
const char* IniConfig::exValueTypeMismatch::what() const throw()
{
	sprintf(Message,FormatString,Value.c_str(),RequestedType.c_str());
	return Message;
}

// IniValue class
IniValue::IniValue(const char* V) : Value(V) {}

IniValue::operator bool() throw(IniConfig::exValueTypeMismatch)
{
	std::stringstream ss(Value);
	ss.setf(std::ios_base::boolalpha);
	bool b; ss >> b;
	if(ss.fail()) throw(IniConfig::exValueTypeMismatch(Value,"bool"));
	return b;
}

IniValue::operator int() throw(IniConfig::exValueTypeMismatch)
{
	std::stringstream ss(Value);
	signed int i; ss >> i;
	if(ss.fail()) throw(Value,IniConfig::exValueTypeMismatch(Value,"int"));
	return i;
}

IniValue::operator double() throw(IniConfig::exValueTypeMismatch) 
{
	std::stringstream ss(Value); 
	double d; ss >> d; 
	if(ss.fail()) throw(IniConfig::exValueTypeMismatch(ss.str(),"double"));
	return d;
}

IniValue::operator std::complex<double>() throw(IniConfig::exValueTypeMismatch) 
{
	std::stringstream ss(Value); 
	std::complex<double> c; ss >> c; 
	if(ss.fail()) throw(IniConfig::exValueTypeMismatch(ss.str(),"std::complex<double>"));
	return c;
}

IniValue::operator std::string() throw(IniConfig::exValueTypeMismatch) { return Value; }

// IniConfig class
IniConfig::IniConfig(std::string IniFileName)
	throw(IniConfig::exCouldNotOpenIniFile)
{
	if((this->ini = iniparser_load(IniFileName.c_str())) == NULL)
		throw(IniConfig::exCouldNotOpenIniFile(IniFileName));
}

IniConfig::~IniConfig()
{
	iniparser_freedict(this->ini);	
}

IniValue& IniConfig::operator[](std::string Key) throw(IniConfig::exKeyNotFound)
{
	const char* entry = Key.c_str();
	if(!iniparser_find_entry(ini,const_cast<char*>(entry)))
		throw(IniConfig::exKeyNotFound(Key));

	IniValue *pValue = new IniValue(iniparser_getstring(ini,entry,NULL));	
	return *pValue;
}

