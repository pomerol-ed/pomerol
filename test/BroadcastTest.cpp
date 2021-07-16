#include <mpi_dispatcher/misc.hpp>

#include <pomerol/TwoParticleGFPart.hpp>

#include "catch2/catch-pomerol.hpp"

#include <algorithm>
#include <vector>

namespace Pomerol {

bool operator==(const TwoParticleGFPart::NonResonantTerm &t1,
                const TwoParticleGFPart::NonResonantTerm &t2) {
    return t1.Coeff == t2.Coeff &&
           std::equal(t1.Poles, t1.Poles + 3, t2.Poles) &&
           t1.isz4 == t2.isz4 &&
           t1.Weight == t2.Weight;
}
bool operator!=(const TwoParticleGFPart::NonResonantTerm &t1,
                const TwoParticleGFPart::NonResonantTerm &t2) {
    return !operator==(t1, t2);
}

bool operator==(const TwoParticleGFPart::ResonantTerm &t1,
                const TwoParticleGFPart::ResonantTerm &t2) {
    return t1.ResCoeff == t2.ResCoeff &&
           t1.NonResCoeff == t2.NonResCoeff &&
           std::equal(t1.Poles, t1.Poles + 3, t2.Poles) &&
           t1.isz1z2 == t2.isz1z2 &&
           t1.Weight == t2.Weight;
}
bool operator!=(const TwoParticleGFPart::ResonantTerm &t1,
                const TwoParticleGFPart::ResonantTerm &t2) {
    return !operator==(t1, t2);
}

}

using namespace Pomerol;

TEST_CASE("broadcast() methods of various objects", "[broadcast]") {
    int rank = pMPI::rank(MPI_COMM_WORLD);

    // Reference objects
    TwoParticleGFPart::NonResonantTerm tnr_ref(ComplexType(4.0, 3.0),
                                               -0.1, 0.2, 0.3, true
    );
    tnr_ref.Weight = 100;
    TwoParticleGFPart::ResonantTerm tr_ref(ComplexType(4.0, 3.0),
                                            ComplexType(5.0, 6.0),
                                          -0.1, 0.2, 0.3, true
    );
    tr_ref.Weight = 100;

    SECTION("TwoParticleGFPart::NonResonantTerm::mpi_datatype()") {
        TwoParticleGFPart::NonResonantTerm tnr(0, 0, 0, 0, false);
        if(rank == 0)
            tnr = tnr_ref;

        std::vector<TwoParticleGFPart::NonResonantTerm> tnr_v(10, tnr);

        MPI_Bcast(tnr_v.data(), tnr_v.size(),
                  TwoParticleGFPart::NonResonantTerm::mpi_datatype(),
                  0, MPI_COMM_WORLD);

        for(auto const& t : tnr_v)
            REQUIRE(t == tnr_ref);
    }

    SECTION("TwoParticleGFPart::ResonantTerm::mpi_datatype()") {
        TwoParticleGFPart::ResonantTerm tr(0, 0, 0, 0, 0, false);

        if(rank == 0)
            tr = tr_ref;

        std::vector<TwoParticleGFPart::ResonantTerm> tr_v(10, tr);
        MPI_Bcast(tr_v.data(), tr_v.size(),
                  TwoParticleGFPart::ResonantTerm::mpi_datatype(),
                  0, MPI_COMM_WORLD);

        for(auto const& t : tr_v)
            REQUIRE(t == tr_ref);
    }

    SECTION("TermList::broadcast()") {
        TermList<TwoParticleGFPart::NonResonantTerm> tl(
            TwoParticleGFPart::NonResonantTerm::Compare(1.0/1024),
            TwoParticleGFPart::NonResonantTerm::IsNegligible(1.0/1024)
        );
        tl.add_term(tnr_ref);

        TermList<TwoParticleGFPart::NonResonantTerm> tl_ref(
            TwoParticleGFPart::NonResonantTerm::Compare(1.0/2048),
            TwoParticleGFPart::NonResonantTerm::IsNegligible(1.0/2048)
        );
        tl_ref.add_term(TwoParticleGFPart::NonResonantTerm(ComplexType(1.0, 2.0),
                                                       -0.1, 0.2, 0.4, true));
        tl_ref.add_term(TwoParticleGFPart::NonResonantTerm(ComplexType(1.0, 8.0),
                                                       -0.4, 0.2, 0.4, false));
        tl_ref.add_term(TwoParticleGFPart::NonResonantTerm(ComplexType(7.0, 2.0),
                                                       -0.6, 0.2, 0.4, true));

        if(rank == 0)
            tl = tl_ref;

        tl.broadcast(MPI_COMM_WORLD, 0);

        REQUIRE(tl.get_is_negligible().Tolerance == tl_ref.get_is_negligible().Tolerance);
        REQUIRE(tl.as_set().key_comp().Tolerance == tl_ref.as_set().key_comp().Tolerance);
        REQUIRE(tl.as_set() == tl_ref.as_set());
    }
}
