//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


#ifndef __INCLUDE_LATTICEANALYSIS_H
#define __INCLUDE_LATTICEANALYSIS_H

#include"Misc.h"
#include<json/json.h>

namespace Pomerol{

struct SiteHoppingElement
 {
    unsigned short From;
    unsigned short OrbitalFrom;
    unsigned short To;
    unsigned short OrbitalTo;
    RealType Value;
    SiteHoppingElement(unsigned short From, unsigned short OrbitalFrom, unsigned short To, unsigned short OrbitalTo, RealType Value) : From(From),OrbitalFrom(OrbitalFrom),To(To),OrbitalTo(OrbitalTo), Value(Value){};
    friend std::ostream& operator<<(std::ostream& output, const SiteHoppingElement& out);
 };

class LatticeSite
{
public:
     unsigned short type;
     unsigned short number;
    RealType LocalMu;
    std::list<SiteHoppingElement*> HoppingList;
    bool operator == (const LatticeSite& right);
    bool is_equivalent (const LatticeSite& right);
};

class sLatticeSite : public LatticeSite
{
public:
    RealType U;
    sLatticeSite(unsigned short type_, RealType LocalMu_, unsigned short number_, RealType U);
    friend std::ostream& operator<<(std::ostream& output, const sLatticeSite& out);
};

class pLatticeSite : public LatticeSite
{
public:
    RealType U;
    RealType J;
    std::string basis;
    pLatticeSite(unsigned short type_, RealType LocalMu_, unsigned short number_, RealType U, RealType J, std::string &basis);
    friend std::ostream& operator<<(std::ostream& output, const pLatticeSite& out);
};

class LatticeAnalysis
{
    Json::Value *root;
    std::map<std::string, std::list<std::list<unsigned short> > > SitesPermutations;
    std::vector<LatticeSite*> SitesList;

    std::map<std::string, OrbitalValue> mapOrbitalValue;     //!< The map between string and value in the previous enum

    void findSitesPermutations();
    void classifySites();
    void enterHoppingListForCurrentSite(unsigned short CurrentSite, Json::Value &Hopping, std::list<SiteHoppingElement*> &HoppingList);

public:
    LatticeAnalysis ();
    int readin(std::string &LatticeFile);
    const std::vector<LatticeSite*>& getSitesList();
    std::stringstream& printSitesList();
};

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_LATTICEANALYSIS_H
