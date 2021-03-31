#include "pomerol/TwoParticleGFPart.h"

#include <cstddef>
#include <mutex>

std::mutex NonResonantTerm_mpi_datatype_mutex;
std::mutex ResonantTerm_mpi_datatype_mutex;

namespace Pomerol{

// Make the lagging index catch up or outrun the leading index.
template<bool Complex>
inline bool chaseIndices(typename RowMajorMatrixType<Complex>::InnerIterator& index1_iter,
                         typename ColMajorMatrixType<Complex>::InnerIterator& index2_iter)
{
    InnerQuantumState index1 = index1_iter.index();
    InnerQuantumState index2 = index2_iter.index();

    if(index1 == index2) return true;

    if(index1 < index2)
        for(;InnerQuantumState(index1_iter.index())<index2 && index1_iter; ++index1_iter);
    else
        for(;InnerQuantumState(index2_iter.index())<index1 && index2_iter; ++index2_iter);

    return false;
}

//
// TwoParticleGFPart::NonResonantTerm
//
template<bool Complex>
auto TwoParticleGFPart<Complex>::NonResonantTerm::operator+=(
                    const NonResonantTerm& AnotherTerm) -> NonResonantTerm&
{
    long combinedWeight=Weight + AnotherTerm.Weight;
    for (unsigned short p=0; p<3; ++p) Poles[p]= (Weight*Poles[p] + AnotherTerm.Weight*AnotherTerm.Poles[p])/combinedWeight;
    Weight=combinedWeight;
    Coeff += AnotherTerm.Coeff;
    return *this;
}

//
// TwoParticleGFPart::ResonantTerm
//
template<bool Complex>
auto TwoParticleGFPart<Complex>::ResonantTerm::operator+=(const ResonantTerm& AnotherTerm) -> ResonantTerm&
{
    long combinedWeight=Weight + AnotherTerm.Weight;
    for (unsigned short p=0; p<3; ++p) Poles[p]= (Weight*Poles[p] + AnotherTerm.Weight*AnotherTerm.Poles[p])/combinedWeight;
    Weight=combinedWeight;
    ResCoeff += AnotherTerm.ResCoeff;
    NonResCoeff += AnotherTerm.NonResCoeff;
    return *this;
}

//
// TwoParticleGFPart
//
template<bool Complex>
TwoParticleGFPart<Complex>::TwoParticleGFPart(
                const FieldOperatorPart<Complex>& O1, const FieldOperatorPart<Complex>& O2,
                const FieldOperatorPart<Complex>& O3, const CreationOperatorPart<Complex>& CX4,
                const HamiltonianPart<Complex>& Hpart1, const HamiltonianPart<Complex>& Hpart2,
                const HamiltonianPart<Complex>& Hpart3, const HamiltonianPart<Complex>& Hpart4,
                const DensityMatrixPart<Complex>& DMpart1, const DensityMatrixPart<Complex>& DMpart2,
                const DensityMatrixPart<Complex>& DMpart3, const DensityMatrixPart<Complex>& DMpart4,
                Permutation3 Permutation) :
    Thermal(DMpart1),
    ComputableObject(),
    NonResonantTerms(typename NonResonantTerm::Compare(1e-8), typename NonResonantTerm::IsNegligible(1e-16)),
    ResonantTerms(typename ResonantTerm::Compare(1e-8), typename ResonantTerm::IsNegligible(1e-16)),
    O1(O1), O2(O2), O3(O3), CX4(CX4),
    Hpart1(Hpart1), Hpart2(Hpart2), Hpart3(Hpart3), Hpart4(Hpart4),
    DMpart1(DMpart1), DMpart2(DMpart2), DMpart3(DMpart3), DMpart4(DMpart4),
    Permutation(Permutation),
    ReduceResonanceTolerance(1e-8),
    CoefficientTolerance (1e-16),
    MultiTermCoefficientTolerance (1e-5)
{}

template<bool Complex>
void TwoParticleGFPart<Complex>::compute()
{
    NonResonantTerms.clear();
    ResonantTerms.clear();

    RealType beta = DMpart1.beta;
    // I don't have any pen now, so I'm writing here:
    // <1 | O1 | 2> <2 | O2 | 3> <3 | O3 |4> <4| CX4 |1>
    // Iterate over all values of |1><1| and |3><3|
    // Chase indices |2> and <2|, |4> and <4|.
    const RowMajorMatrixType<Complex>& O1matrix = O1.getRowMajorValue();
    const ColMajorMatrixType<Complex>& O2matrix = O2.getColMajorValue();
    const RowMajorMatrixType<Complex>& O3matrix = O3.getRowMajorValue();
    const ColMajorMatrixType<Complex>& CX4matrix = CX4.getColMajorValue();

    InnerQuantumState index1;
    InnerQuantumState index1Max = CX4matrix.outerSize(); // One can not make a cutoff in external index for evaluating 2PGF

    InnerQuantumState index3;
    InnerQuantumState index3Max = O2matrix.outerSize();

    std::vector<InnerQuantumState> Index4List;
    Index4List.reserve(index1Max*index3Max);

    for(index1=0; index1<index1Max; ++index1)
    for(index3=0; index3<index3Max; ++index3){
        typename ColMajorMatrixType<Complex>::InnerIterator index4bra_iter(CX4matrix,index1);
        typename RowMajorMatrixType<Complex>::InnerIterator index4ket_iter(O3matrix,index3);
        Index4List.clear();
        while (index4bra_iter && index4ket_iter){
            if(chaseIndices<Complex>(index4ket_iter,index4bra_iter)){
                Index4List.push_back(index4bra_iter.index());
                ++index4bra_iter;
                ++index4ket_iter;
            }
        };

        if (!Index4List.empty())
        {
            RealType E1 = Hpart1.getEigenValue(index1);
            RealType E3 = Hpart3.getEigenValue(index3);
            RealType weight1 = DMpart1.getWeight(index1);
            RealType weight3 = DMpart3.getWeight(index3);

            typename ColMajorMatrixType<Complex>::InnerIterator index2bra_iter(O2matrix,index3);
            typename RowMajorMatrixType<Complex>::InnerIterator index2ket_iter(O1matrix,index1);
            while (index2bra_iter && index2ket_iter){
                if (chaseIndices<Complex>(index2ket_iter,index2bra_iter)){

                    InnerQuantumState index2 = index2ket_iter.index();
                    RealType E2 = Hpart2.getEigenValue(index2);
                    RealType weight2 = DMpart2.getWeight(index2);

                    for (unsigned long p4 = 0; p4 < Index4List.size(); ++p4)
                    {
                        InnerQuantumState index4 = Index4List[p4];//*pIndex4;
                        RealType E4 = Hpart4.getEigenValue(index4);
                        RealType weight4 = DMpart4.getWeight(index4);
                        if (weight1 + weight2 + weight3 + weight4 >= CoefficientTolerance) {
                            ComplexType MatrixElement = index2ket_iter.value()*
                                                        index2bra_iter.value()*
                                                        O3matrix.coeff(index3,index4)*
                                                        CX4matrix.coeff(index4,index1);

                            MatrixElement *= Permutation.sign;

                            addMultiterm(MatrixElement,beta,E1,E2,E3,E4,weight1,weight2,weight3,weight4);
                        }
                    }
                    ++index2bra_iter;
                    ++index2ket_iter;
                };
            }
        };
    }

    std::cout << "Total " << NonResonantTerms.size() << "+" << ResonantTerms.size() << "="
              << NonResonantTerms.size() + ResonantTerms.size() << " terms" << std::endl << std::flush;

    assert(NonResonantTerms.check_terms());
    assert(ResonantTerms.check_terms());

    Status = Computed;
}

template<bool Complex>
inline void TwoParticleGFPart<Complex>::addMultiterm(ComplexType Coeff, RealType beta,
                      RealType Ei, RealType Ej, RealType Ek, RealType El,
                      RealType Wi, RealType Wj, RealType Wk, RealType Wl)
{
    RealType P1 = Ej - Ei;
    RealType P2 = Ek - Ej;
    RealType P3 = El - Ek;

    // Non-resonant part of the multiterm
    ComplexType CoeffZ2 = -Coeff*(Wj + Wk);
    if(abs(CoeffZ2) > CoefficientTolerance)
        NonResonantTerms.add_term(
            NonResonantTerm(CoeffZ2,P1,P2,P3,false));
    ComplexType CoeffZ4 = Coeff*(Wi + Wl);
    if(abs(CoeffZ4) > CoefficientTolerance)
        NonResonantTerms.add_term(
            NonResonantTerm(CoeffZ4,P1,P2,P3,true));

    // Resonant part of the multiterm
    ComplexType CoeffZ1Z2Res = Coeff*beta*Wi;
    ComplexType CoeffZ1Z2NonRes = Coeff*(Wk - Wi);
    if(abs(CoeffZ1Z2Res) > CoefficientTolerance || abs(CoeffZ1Z2NonRes) > CoefficientTolerance)
        ResonantTerms.add_term(
            ResonantTerm(CoeffZ1Z2Res,CoeffZ1Z2NonRes,P1,P2,P3,true));
    ComplexType CoeffZ2Z3Res = -Coeff*beta*Wj;
    ComplexType CoeffZ2Z3NonRes = Coeff*(Wj - Wl);
    if(abs(CoeffZ2Z3Res) > CoefficientTolerance || abs(CoeffZ2Z3NonRes) > CoefficientTolerance)
        ResonantTerms.add_term(
            ResonantTerm(CoeffZ2Z3Res,CoeffZ2Z3NonRes,P1,P2,P3,false));
}

template<bool Complex>
size_t TwoParticleGFPart<Complex>::getNumNonResonantTerms() const
{
    return NonResonantTerms.size();
}

template<bool Complex>
size_t TwoParticleGFPart<Complex>::getNumResonantTerms() const
{
    return ResonantTerms.size();
}

template<bool Complex>
const Permutation3& TwoParticleGFPart<Complex>::getPermutation() const
{
    return Permutation;
}

template<bool Complex>
ComplexType TwoParticleGFPart<Complex>::operator()(long MatsubaraNumber1, long MatsubaraNumber2, long MatsubaraNumber3) const
{
    long MatsubaraNumberOdd1 = 2*MatsubaraNumber1 + 1;
    long MatsubaraNumberOdd2 = 2*MatsubaraNumber2 + 1;
    long MatsubaraNumberOdd3 = 2*MatsubaraNumber3 + 1;
    return (*this)(MatsubaraSpacing * RealType(MatsubaraNumberOdd1),
                   MatsubaraSpacing * RealType(MatsubaraNumberOdd2),
                   MatsubaraSpacing * RealType(MatsubaraNumberOdd3));
}

template<bool Complex>
ComplexType TwoParticleGFPart<Complex>::operator()(ComplexType z1, ComplexType z2, ComplexType z3) const
{
    ComplexType Frequencies[3] = {  z1, z2, -z3 };

    z1 = Frequencies[Permutation.perm[0]];
    z2 = Frequencies[Permutation.perm[1]];
    z3 = Frequencies[Permutation.perm[2]];

    if (Status != Computed) {
        throw std::logic_error("2PGFPart : Calling operator() on uncomputed container, did you purge all the terms when called compute()");
    }

    return NonResonantTerms(z1, z2, z3) + ResonantTerms(z1, z2, z3, ReduceResonanceTolerance);
}

template<bool Complex>
auto TwoParticleGFPart<Complex>::getResonantTerms() const -> const TermList<ResonantTerm>&
{
    return ResonantTerms;
}

template<bool Complex>
auto TwoParticleGFPart<Complex>::getNonResonantTerms() const -> const TermList<NonResonantTerm>&
{
    return NonResonantTerms;
}

template<bool Complex>
void TwoParticleGFPart<Complex>::clear()
{
    NonResonantTerms.clear();
    ResonantTerms.clear();
    Status = Constructed;
}

template class TwoParticleGFPart<false>;
template class TwoParticleGFPart<true>;

} // end of namespace Pomerol

namespace pMPI {

// When called for the first time, this function creates an MPI structure datatype
// describing TwoParticleGFPart::NonResonantTerm and registers it using MPI_Type_commit().
// The registered datatype is stored in a static variable and is immediately returned
// upon successive calls to mpi_datatype().
//
// About MPI structure datatypes:
// https://pages.tacc.utexas.edu/~eijkhout/pcse/html/mpi-data.html#Structtype
template<bool Complex>
MPI_Datatype mpi_datatype_tpgfp_non_resonant_term() {
    static MPI_Datatype dt;

    // Since we are using static variables here, we have to make sure the code
    // is thread-safe.
    const std::lock_guard<std::mutex> lock(NonResonantTerm_mpi_datatype_mutex);

    // Create and commit datatype only once
    static bool type_committed = false;
    if(!type_committed) {
        int blocklengths[] = {1,3,1,1};
        using TermT = typename Pomerol::TwoParticleGFPart<Complex>::NonResonantTerm;
        MPI_Aint displacements[] = {offsetof(TermT, Coeff),
                                    offsetof(TermT, Poles),
                                    offsetof(TermT, isz4),
                                    offsetof(TermT, Weight)
                                  };
        MPI_Datatype types[] = {MPI_CXX_DOUBLE_COMPLEX, // ComplexType Coeff
                                MPI_DOUBLE,             // RealType Poles[3]
                                MPI_CXX_BOOL,           // bool isz4
                                MPI_LONG                // long Weight
                              };
        MPI_Type_create_struct(4, blocklengths, displacements, types, &dt);
        MPI_Type_commit(&dt);
        type_committed = true;
    }
    return dt;
}

MPI_Datatype mpi_datatype<Pomerol::TwoParticleGFPart<false>::NonResonantTerm>::get() {
  return mpi_datatype_tpgfp_non_resonant_term<false>();
}
MPI_Datatype mpi_datatype<Pomerol::TwoParticleGFPart<true>::NonResonantTerm>::get() {
  return mpi_datatype_tpgfp_non_resonant_term<true>();
}

// When called for the first time, this function creates an MPI structure datatype
// describing TwoParticleGFPart::ResonantTerm and registers it using MPI_Type_commit().
// The registered datatype is stored in a static variable and is immediately returned
// upon successive calls to mpi_datatype().
//
// About MPI structure datatypes:
// https://pages.tacc.utexas.edu/~eijkhout/pcse/html/mpi-data.html#Structtype
template<bool Complex>
MPI_Datatype mpi_datatype_tpgfp_resonant_term() {
    static MPI_Datatype dt;

    // Since we are using static variables here, we have to make sure the code
    // is thread-safe.
    const std::lock_guard<std::mutex> lock(ResonantTerm_mpi_datatype_mutex);

    // Create and commit datatype only once
    static bool type_committed = false;
    if(!type_committed) {
        int blocklengths[] = {1,1,3,1,1};
        using TermT = typename Pomerol::TwoParticleGFPart<Complex>::ResonantTerm;
        MPI_Aint displacements[] = {offsetof(TermT, ResCoeff),
                                    offsetof(TermT, NonResCoeff),
                                    offsetof(TermT, Poles),
                                    offsetof(TermT, isz1z2),
                                    offsetof(TermT, Weight)
                                   };
        MPI_Datatype types[] = {MPI_CXX_DOUBLE_COMPLEX, // ComplexType ResCoeff
                                MPI_CXX_DOUBLE_COMPLEX, // ComplexType NonResCoeff
                                MPI_DOUBLE,             // RealType Poles[3]
                                MPI_CXX_BOOL,           // bool isz1z2
                                MPI_LONG                // long Weight
                               };
        MPI_Type_create_struct(5, blocklengths, displacements, types, &dt);
        MPI_Type_commit(&dt);
        type_committed = true;
    }
    return dt;
}

MPI_Datatype mpi_datatype<Pomerol::TwoParticleGFPart<false>::ResonantTerm>::get() {
  return mpi_datatype_tpgfp_resonant_term<false>();
}
MPI_Datatype mpi_datatype<Pomerol::TwoParticleGFPart<true>::ResonantTerm>::get() {
  return mpi_datatype_tpgfp_resonant_term<true>();
}

} // end of namespace pMPI
