#include "pomerol/EnsembleAverage.h"

namespace Pomerol{

EnsembleAverage::EnsembleAverage(const StatesClassification& S, const Hamiltonian& H,
                                 const QuadraticOperator& A, const DensityMatrix& DM) :
    Thermal(DM.beta), ComputableObject(), S(S), H(H), A(A), DM(DM), Vanishing(true)
{
}

EnsembleAverage::EnsembleAverage(const EnsembleAverage& EA) :
    Thermal(EA.beta), ComputableObject(EA), S(EA.S), H(EA.H), A(EA.A), DM(EA.DM), Vanishing(EA.Vanishing)
{
    for(std::list<EnsembleAveragePart*>::const_iterator iter = EA.parts.begin(); iter != EA.parts.end(); iter++)
        parts.push_back(new EnsembleAveragePart(**iter));
}

EnsembleAverage::~EnsembleAverage()
{
    for(std::list<EnsembleAveragePart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
        delete *iter;
}

void EnsembleAverage::prepare(void)
{
    if(Status>=Prepared) return;

    // Find out non-trivial blocks of A.
    FieldOperator::BlocksBimap const& ANontrivialBlocks = A.getBlockMapping();

    for(FieldOperator::BlocksBimap::left_const_iterator Aiter = ANontrivialBlocks.left.begin(); Aiter != ANontrivialBlocks.left.end(); Aiter++){
        // <Aleft|A|Aright>
        BlockNumber Aleft = Aiter->first;
        BlockNumber Aright = Aiter->second;

        // Only diagonal blocks
        if(Aleft == Aright){
            DEBUG(S.getQuantumNumbers(Aleft) << "|" << S.getQuantumNumbers(Aright) );
            // check if retained blocks are included. If not, do not push.
            if ( DM.isRetained(Aleft) )
                parts.push_back(new EnsembleAveragePart((QuadraticOperatorPart&)A.getPartFromLeftIndex(Aleft),
                                                        H.getPart(Aleft), DM.getPart(Aleft)));
        }
    }
    if (parts.size() > 0) Vanishing = false;

    Status = Prepared;
}

void EnsembleAverage::compute()
{
    if(Status>=Computed) return;
    if(Status<Prepared) prepare();

    if(Status<Computed){
        for(std::list<EnsembleAveragePart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
            (*iter)->compute();
    }
    Status = Computed;
}

bool EnsembleAverage::isVanishing(void) const
{
    return Vanishing;
}

} // end of namespace Pomerol
