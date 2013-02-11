/*
 *		An event-driven parser for command-line arguments.
 *  
 *		Copyright (c) 2004-2005 by N.Okazaki
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions (known as zlib license):
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 *
 * Naoaki Okazaki <okazaki at chokkan.org>
 *
 */

/* $Id: optparse.h 2 2006-10-31 00:57:57Z naoaki $ */

/*
 * Class 'optparse' implements a parser for GNU-style command-line arguments.
 * Inherit this class to define your own option variables and to implement an
 * option handler with macros, BEGIN_OPTION_MAP, ON_OPTION(_WITH_ARG), and
 * END_OPTION_MAP. Consult the sample program attached at the bottom of this
 * source code.
 *
 * This code was comfirmed to be compiled with MCVC++ 2003 and gcc 3.3.
 * Define _BUILD_NCL_SAMPLE if you want to build a sample program.
 *	$ g++ -D_BUILD_NCL_SAMPLE -xc++ optparse.h
 */

#ifndef	__INCLUDE_OPTIONPARSER_H
#define	__INCLUDE_OPTIONPARSER_H

#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>

#ifdef	USE_NCL_NAMESPACE
namespace ncl {
#endif/*USE_NCL_NAMESPACE*/


/**
 * An event-driven parser for command-line arguments.
 *	@author	Naoaki Okazaki
 */
class optparse {
public:
	/**
	 * Exception class for unrecognized options.
	 */
	class unrecognized_option : public std::invalid_argument {
	public:
		unrecognized_option(char shortopt)
			: std::invalid_argument(std::string("-") + shortopt) {}
		unrecognized_option(const std::string& longopt)
			: std::invalid_argument(std::string("--") + longopt) {}
	};

	/**
	 * Exception class for invalid values.
	 */
	class invalid_value : public std::invalid_argument {
	public:
		std::string optionstr;

		invalid_value(const std::string& message)
			: std::invalid_argument(message) {}
		invalid_value(char shortopt, const char *longopt, const std::string& message) :
			std::invalid_argument(message),
			optionstr(
				shortopt ?
					(std::string("-") + shortopt) : 
					(longopt ? (std::string("--") + longopt) : std::string(""))
				) {}

        ~invalid_value() throw (){};
		const std::string& option() const {return optionstr; }
	};

public:
	/** Construct. */
	optparse() {}
	/** Destruct. */
	virtual ~optparse() {}

	/**
	 * Parse options.
	 *	@param	argv		array of null-terminated strings to be parsed
	 *	@param	num_argv	specifies the number, in strings, of the array
	 *	@return				the number of used arguments
	 *	@throws				optparse_exception
	 */
	int parse(char * const argv[], int num_argv)
	{
		int i;
		for (i = 0;i < num_argv;++i) {
			const char *token = argv[i];
			if (*token++ == '-') {
				const char *next_token = (i+1 < num_argv) ? argv[i+1] : "";
				if (!*token) {
					break;	// only '-' was found.
				} else if (*token == '-') {
					const char *arg = std::strchr(++token, '=');
					if (arg) {
						arg++;
					} else {
						arg = next_token;
					}
					int ret = handle_option(0, token, arg);
					if (ret < 0) {
						throw unrecognized_option(token);
					}
					if (arg == next_token) {
						i += ret;
					}
				} else {
					char c;
					while ((c = *token++) != '\0') {
						const char *arg = *token ? token : next_token;
						int ret = handle_option(c, token, arg);
						if (ret < 0) {
							throw unrecognized_option(c);
						}
						if (ret > 0) {
							if (arg == token) {
								token = "";
							} else {
								i++;
							}
						}
					} // while
				} // else (*token == '-') 
			} else {
				break;	// a non-option argument was fonud.
			} 
		} // for (i)

		return i;
	}

protected:
	/**
	 * Option handler
	 *	This function should be overridden by inheritance class.
	 *	@param	c			short option character, 0 for long option
	 *	@param	longname	long option name
	 *	@param	arg			an argument for the option
	 *	@return				0 (success);
							1 (success with use of an argument);
							-1 (failed, unrecognized option)
	 *	@throws				option_parser_exception
	 */
	virtual int handle_option(char c, const char *longname, const char *arg)
	{
		return 0;
	}

	int __optstrcmp(const char *option, const char *longname)
	{
		const char *p = std::strchr(option, '=');
		return p ?
			std::strncmp(option, longname, p-option) :
			std::strcmp(option, longname);
	}
};


/** The begin of inline option map. */
#define	BEGIN_OPTION_MAP_INLINE() \
	virtual int handle_option(char __c, const char *__longname, const char *arg) \
	{ \
		int used_args = 0; \
		if (0) { \

/** Define of option map. */
#define	DEFINE_OPTION_MAP() \
	virtual int handle_option(char __c, const char *__longname, const char *arg);

/** Begin of option map implimentation. */
#define	BEGIN_OPTION_MAP(_Class) \
	int _Class::handle_option(char __c, const char *__longname, const char *arg) \
	{ \
		int used_args = 0; \
		if (0) { \

/** An entry of option map */
#define	ON_OPTION(test) \
			return used_args; \
		} else if (test) { \
			used_args = 0; \

#define	ON_OPTION_WITH_ARG(test) \
			return used_args; \
		} else if (test) { \
			used_args = 1; \

/** The end of option map implementation */
#define	END_OPTION_MAP() \
			return used_args; \
		} \
		return -1; \
	} \

/** A predicator for short options */
#define	SHORTOPT(x)		(__c == x)
/** A predicator for long options */
#define	LONGOPT(x)		(!__c && __optstrcmp(__longname, x) == 0)


#ifdef	USE_NCL_NAMESPACE
};
#endif/*USE_NCL_NAMESPACE*/

#include"Misc.h"

/**
 * A class to store parameters specified by command-line arguments
 */
class pomerolOptionParser : public optparse {
public:
	Pomerol::RealType beta ;
	unsigned long NumberOfMatsubaras;
	std::string LatticeFile;
    bool calculateGF;
    bool calculate2PGF;
    bool savePlaintext;
	std::string help;

	pomerolOptionParser() : 
        beta(10), 
        NumberOfMatsubaras(60), 
        LatticeFile("Lattice.json"), 
        calculateGF(false), 
        calculate2PGF(false),
        savePlaintext(false),
        help("") 
     {}

	BEGIN_OPTION_MAP_INLINE()
		ON_OPTION(SHORTOPT('b') || LONGOPT("beta"))
			beta = std::atof(arg);
            //NumberOfMatsubaras = (int)std::fabs(beta+0.5);
			used_args = 1;	// Notify the parser of a consumption of argument.

		ON_OPTION_WITH_ARG(SHORTOPT('m') || LONGOPT("matsubaras"))
			NumberOfMatsubaras = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

		ON_OPTION_WITH_ARG(SHORTOPT('l') || LONGOPT("lattice") || LONGOPT("Lattice"))
			LatticeFile = arg;
			used_args = 1;	// Notify the parser of a consumption of argument.

		ON_OPTION(LONGOPT("calcgf") || LONGOPT("calculategf"))
            calculateGF=true;

		ON_OPTION(LONGOPT("calc2pgf") || LONGOPT("calculate2pgf"))
            calculate2PGF = true;
            calculateGF = true;

        ON_OPTION(LONGOPT("plaintext"))
            savePlaintext = true;

        ON_OPTION(SHORTOPT('h') || LONGOPT("help"))
            std::cout << "pomerolDiag - an ED code, which provides one- and two- particle Greens functions and irreducible vertex part in Matsubara domain" << std::endl;
            std::cout << "Usage: pomerolDiag [options]" << std::endl;
            std::cout << "Options: " << std::endl;
            std::cout << "-b     --beta        : The value of inverse temperature. Default: " << beta << std::endl;
            std::cout << "-m     --matsubaras  : Amount of Matsubara frequencies. Default: " << NumberOfMatsubaras<< std::endl;
            std::cout << "-l     --lattice     : A file with the lattice. Default : " << LatticeFile << std::endl;
            std::cout << "-h     --help        : Show this help message" << std::endl;
            std::cout << "--calculategf        : Defines whether the program will calculate a Green's function. Default: false." << std::endl;
            std::cout << "--calculate2pgf      : Defines whether the program will calculate a vertex. Default: false." << std::endl;
            exit(0);

	END_OPTION_MAP()
};

#endif // endif :: #ifndef __INCLUDE_OPTIONPARSER_H
