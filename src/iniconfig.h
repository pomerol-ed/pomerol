#ifndef __INICONFIG__
#define __INICONFIG__ 

/** \file iniconfig.h
**	\brief Declaration of IniConfig and IniValue classes.
**	A brief description of the ini-file format is here: http://ndevilla.free.fr/iniparser/html/index.html#inidef
** 
** \author	Igor Krivenko (igor@shg.ru)
*/

#include<string>
#include<complex>

// These classes are actually C++ wrappers around the iniparser library.
extern "C"{
	#include<iniparser.h>
}


class IniValue;

/**
  This class represents an ini-file
*/
class IniConfig{

	/** Handle to the ini-file. */
	dictionary* ini;	

public:

	//
	// Exceptions
	//

	/** Exception: could not open an ini-file file. */
	class exCouldNotOpenIniFile : public std::exception{
		char* Message;
		const char* FormatString;
		std::string FileName;
	public:
		exCouldNotOpenIniFile(std::string Name) throw();
		~exCouldNotOpenIniFile() throw();
		const char* what() const throw();
	};

	/** Exception: key not found. */
	class exKeyNotFound : public std::exception{
		char* Message;
		const char* FormatString;
		std::string Key;
	public:
		exKeyNotFound(std::string Key) throw();
		~exKeyNotFound() throw();
		const char* what() const throw();
	};

	/** Exception: value type mismatch (unable to cast the value to the requested type). */
	class exValueTypeMismatch : public std::exception{
		char* Message;
		const char* FormatString;
		std::string RequestedType;
		std::string Value;
	public:
		exValueTypeMismatch(std::string V, std::string Requested) throw();
		~exValueTypeMismatch() throw();
		const char* what() const throw();
	};

	/** Open an ini-file. */
	IniConfig(std::string IniFileName)
		throw(IniConfig::exCouldNotOpenIniFile);
	
	/** Close the ini-file and free all associated resources. */
	~IniConfig();

	/** This operator queries the ini-file for a value. */
	IniValue& operator [](std::string Key) throw(IniConfig::exKeyNotFound); 
};

/** IniInput::operator[] returns an instance of this class.
    This generic object can convert itself to a number of integral types.
*/ 
class IniValue{
	friend class IniConfig;
	
	/** Value as a string.*/
	std::string Value;
	
	/** Construct IniValue from a string. */
	IniValue(const char* V);
public:
	/** Cast IniValue to an integer. */
	operator int() throw(IniConfig::exValueTypeMismatch);
	/** Cast IniValue to a boolean value. 
	Supported values in the ini-file are 'true' and 'false'. */
	operator bool() throw(IniConfig::exValueTypeMismatch);
	/** Cast IniValue to a double-precision number. */
	operator double() throw(IniConfig::exValueTypeMismatch);
	/** Cast IniValue to a string. */
	operator std::string() throw(IniConfig::exValueTypeMismatch);
	/** Cast IniValue to a complex number. In the ini-file complex numbers are defined as '(re,im)'. */
	operator std::complex<double>() throw(IniConfig::exValueTypeMismatch);
};

 #endif // endif :: ___OUTPATH_H___