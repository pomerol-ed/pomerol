// pomerol/trunk/BitClassification.h
// This file is a part of pomerol diagonalization code

/** \file BitClassification.h
**  \brief Declaration of BitInfo & BitClassification classes.
** 
**  \author	Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __BIT_CLASSIFICATION__
#define __BIT_CLASSIFICATION__

#include "config.h"
#include <json/json.h>
#include <map>
#include <string>


//! BitInfo - a class to reproduce full information about the given bit
/**
 * BitInfo class is an abstract class to handle all info about current bit.  
 */ 

class BitInfo
{
public:
  	unsigned short site;
	unsigned short spin;
	string type;
	unsigned short bitNumber;
	void setBitNumber(const unsigned short &in);
	virtual void print_to_screen()=0;
};

class sBitInfo : public BitInfo
{
public:
	RealType U;
	sBitInfo(unsigned short site_, string &type_, unsigned short spin_, RealType U_):U(U_){site=site_;type=type_;spin=spin_;};
	friend std::ostream& operator<<(std::ostream& output, const sBitInfo& out);
	void print_to_screen(){cout << *this << endl;};
};

class pBitInfo : public BitInfo
{
public:
	RealType U;
	RealType J;
	string basis;
	short index; 
	pBitInfo(unsigned short site_, string &type_, unsigned short spin_, short index_, const string &basis_, RealType U_, RealType J_):U(U_),J(J_){site=site_;type=type_;spin=spin_;index=index_;basis=basis_;};
	friend std::ostream& operator<<(std::ostream& output, const pBitInfo& out);
	void print_to_screen(){cout << *this << endl;};
};



class BitClassification
{
  	Json::Value root;   // will contains the root value after parsing.
	int N_bit;
	RealMatrixType HoppingMatrix;
	vector<BitInfo*> BitInfoList;
public:
	BitClassification();
	int readin();
	void printBitInfoList();
	void printHoppingMatrix();
	RealMatrixType& getHoppingMatrix();
	const int& getBitSize() const;
private:
	enum OrbitalValue {s=0, p=1, d=2, f=3};
	std::map<std::string, OrbitalValue> mapOrbitalValue;
	void defineBits();
	void defineHopping();
	vector<unsigned short>& findBits(const unsigned short &site);
	unsigned short findBit(const unsigned short &site,const unsigned short &spin);
};

#endif

