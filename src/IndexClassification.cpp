//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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


#include "IndexClassification.h"
#include <boost/functional/hash.hpp>

namespace Pomerol{

//
//IndexClassification::IndexInfo 
//

IndexClassification::IndexInfo::IndexInfo( const std::string &SiteLabel, const unsigned short Orbital, const unsigned short Spin):
    SiteLabel(SiteLabel), Orbital(Orbital), Spin(Spin)
{
    boost::hash<std::string> string_hash;
    SiteLabelHash = string_hash(SiteLabel);
};

bool IndexClassification::IndexInfo::operator<(const IndexClassification::IndexInfo& rhs) const 
{
    if (SiteLabelHash != rhs.SiteLabelHash) return (SiteLabelHash < rhs.SiteLabelHash);
    if (Orbital != rhs.Orbital) return (Orbital < rhs.Orbital);
    return (Spin < rhs.Spin);
}

std::ostream& operator<<(std::ostream& output, const IndexClassification::IndexInfo& out){
    output << "(" << out.SiteLabel << "," << out.Orbital << "," << out.Spin << ")" ;
    return output;
};

//
//IndexClassification
//

IndexClassification::IndexClassification ( const Lattice::SiteMap &Sites ) : IndexSize(0),Sites(Sites)
{
};

void IndexClassification::prepare(bool order_spins)
{
    unsigned int MaxSpinSize=0;
    for (Lattice::SiteMap::const_iterator it1 = Sites.begin(); it1!=Sites.end();++it1) { // first run : determine IndexSpace size & calculate number of spins on each site.
        IndexSize+= (*(it1->second)).OrbitalSize*(*(it1->second)).SpinSize;
        MaxSpinSize=((*(it1->second)).SpinSize>MaxSpinSize)?(*(it1->second)).SpinSize:MaxSpinSize;
        };

    // Split different spins in one group - useful for spin-symmetric cases

    ParticleIndex currentIndex=0;
    IndicesToInfo.resize(IndexSize);

    if (order_spins) {
        for (unsigned int z=0; z<MaxSpinSize; ++z) {
            for (Lattice::SiteMap::const_iterator it1 = Sites.begin(); it1!=Sites.end();++it1) {
                if (z>=(*(it1->second)).SpinSize) break;
                for (unsigned int i=0; i<(*(it1->second)).OrbitalSize; ++i) {
                        IndicesToInfo[currentIndex] = new IndexInfo( it1->first, i, z);
                        currentIndex++;
                        }; // end of orbital loop
                    }; // end of Lattice::SiteMap loop
                }; // end of spin loop
        } // end of order_spins == true
    else {
        for (Lattice::SiteMap::const_iterator it1 = Sites.begin(); it1!=Sites.end();++it1) {
            for (unsigned int i=0; i<(*(it1->second)).OrbitalSize; ++i) {
                for (size_t z=0; z<(*(it1->second)).SpinSize; ++z) {
                        IndicesToInfo[currentIndex] = new IndexInfo( it1->first, i, z);
                        currentIndex++;
                        }; // end of spin loop
                    }; // end of orbital loop
                }; // end of Lattice::SiteMap loop
        }; // end of order_spins == false

    for (ParticleIndex i=0; i<IndexSize; ++i) InfoToIndices[(*(IndicesToInfo[i]))]=i;
}

const ParticleIndex IndexClassification::getIndexSize() const
{
    return IndexSize;
}

bool IndexClassification::checkIndex(ParticleIndex in)
{
    return (in<IndexSize);
}


void IndexClassification::printIndices()
{
    for (ParticleIndex i=0; i<IndexSize; ++i) INFO("Index " << i << " = " << *(IndicesToInfo[i]));
}

ParticleIndex IndexClassification::getIndex(const std::string &Site, const unsigned short &Orbital, const unsigned short &Spin) const
{
    return getIndex(IndexInfo(Site,Orbital,Spin));
}

ParticleIndex IndexClassification::getIndex(const IndexClassification::IndexInfo &in) const
{
    std::map<IndexInfo, ParticleIndex>::const_iterator it=InfoToIndices.find(in); 
    if (it!=InfoToIndices.end()) return (*it).second;
    else return IndexSize;
}

/*
boost::tuple<std::string, unsigned short, unsigned short> IndexClassification::getInfo(ParticleIndex in) const 
{
    if (in >= IndexSize) throw (exWrongIndex());
    IndexInfo *out = IndicesToInfo[in];
    return boost::make_tuple(out->SiteLabel, out->Orbital, out->Spin);
}*/

IndexClassification::IndexInfo IndexClassification::getInfo(ParticleIndex in) const
{
    if (in >= IndexSize) throw (exWrongIndex());
    IndexInfo *out = IndicesToInfo[in];
    return *out;
}

const char* IndexClassification::exWrongIndex::what() const throw(){
    return "Wrong index";
};

} // end of namespace Pomerol
