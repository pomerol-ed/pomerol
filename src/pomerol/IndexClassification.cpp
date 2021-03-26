#include "pomerol/IndexClassification.h"
#include <boost/functional/hash.hpp>

namespace Pomerol{

//
// IndexInfo
//

IndexInfo::IndexInfo( const std::string &SiteLabel, const unsigned short Orbital, const unsigned short Spin):
    SiteLabel(SiteLabel), Orbital(Orbital), Spin(Spin)
{
    boost::hash<std::string> string_hash;
    SiteLabelHash = string_hash(SiteLabel);
};

bool IndexInfo::operator<(const IndexInfo& rhs) const
{
    if (SiteLabelHash != rhs.SiteLabelHash) return (SiteLabelHash < rhs.SiteLabelHash);
    if (Orbital != rhs.Orbital) return (Orbital < rhs.Orbital);
    return (Spin < rhs.Spin);
}

std::ostream& operator<<(std::ostream& output, const IndexInfo& out){
    output << "(" << out.SiteLabel << "," << out.Orbital << "," << out.Spin << ")" ;
    return output;
};

//
// IndexClassification
//

template<bool Complex>
IndexClassification<Complex>::IndexClassification (const typename Lattice<Complex>::SiteMap &Sites) : IndexSize(0),Sites(Sites)
{
};


template<bool Complex>
void IndexClassification<Complex>::prepare(bool order_spins)
{
    unsigned int MaxSpinSize=0;
    for (auto it1 = Sites.begin(); it1!=Sites.end();++it1) { // first run : determine IndexSpace size & calculate number of spins on each site.
        IndexSize+= (*(it1->second)).OrbitalSize*(*(it1->second)).SpinSize;
        MaxSpinSize=((*(it1->second)).SpinSize>MaxSpinSize)?(*(it1->second)).SpinSize:MaxSpinSize;
        };

    // Split different spins in one group - useful for spin-symmetric cases

    ParticleIndex currentIndex=0;
    IndicesToInfo.resize(IndexSize);

    if (order_spins) {
        for (unsigned int z=0; z<MaxSpinSize; ++z) {
            for (auto it1 = Sites.begin(); it1!=Sites.end();++it1) {
                if (z>=(*(it1->second)).SpinSize) break;
                for (unsigned int i=0; i<(*(it1->second)).OrbitalSize; ++i) {
                        IndicesToInfo[currentIndex] = new IndexInfo( it1->first, i, z);
                        currentIndex++;
                        }; // end of orbital loop
                    }; // end of Lattice::SiteMap loop
                }; // end of spin loop
        } // end of order_spins == true
    else {
        for (auto it1 = Sites.begin(); it1!=Sites.end();++it1) {
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

template<bool Complex>
const ParticleIndex IndexClassification<Complex>::getIndexSize() const
{
    return IndexSize;
}

template<bool Complex>
bool IndexClassification<Complex>::checkIndex(ParticleIndex in)
{
    return (in<IndexSize);
}

template<bool Complex>
void IndexClassification<Complex>::printIndices()
{
    for (ParticleIndex i=0; i<IndexSize; ++i) INFO("Index " << i << " = " << *(IndicesToInfo[i]));
}

template<bool Complex>
ParticleIndex IndexClassification<Complex>::getIndex(const std::string &Site, const unsigned short &Orbital, const unsigned short &Spin) const
{
    return getIndex(IndexInfo(Site,Orbital,Spin));
}

template<bool Complex>
ParticleIndex IndexClassification<Complex>::getIndex(const IndexInfo &in) const
{
    std::map<IndexInfo, ParticleIndex>::const_iterator it=InfoToIndices.find(in);
    if (it!=InfoToIndices.end()) return (*it).second;
    else return IndexSize;
}

template<bool Complex>
IndexInfo IndexClassification<Complex>::getInfo(ParticleIndex in) const
{
    if (in >= IndexSize) throw (exWrongIndex());
    IndexInfo *out = IndicesToInfo[in];
    return *out;
}

template<bool Complex>
const char* IndexClassification<Complex>::exWrongIndex::what() const throw(){
    return "Wrong index";
};

} // end of namespace Pomerol
