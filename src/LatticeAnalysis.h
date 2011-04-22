#ifndef __INCLUDE_LATTICEANALYSIS_H
#define __INCLUDE_LATTICEANALYSIS_H

#include"Misc.h"
#include<json/json.h>

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

#endif // endif :: #ifndef __INCLUDE_LATTICEANALYSIS_H
