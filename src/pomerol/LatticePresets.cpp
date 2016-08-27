//
// This file is a part of pomerol - a scientific ED code for obtaining
// properties of a Hubbard model on a finite-size lattice
//
// Copyright (C) 2010-2011 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2011 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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

/** \file src/Lattice.h
** \brief A lattice handler. Reads and stores the lattice from a JSON file. Designed to be extended to other formats.
**
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include "pomerol/LatticePresets.h"

namespace Pomerol {

//
// Lattice::Term::Presets
//

// 2-index presets

Lattice::Term* Lattice::Term::Presets::Hopping ( const std::string& Label1, const std::string& Label2, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2 )
{
    Lattice::Term *T = new Term(2);
    bool OperatorSequence[2] = { 1, 0 };
    std::string Labels[2] = { Label1, Label2 };
    unsigned short Orbitals[2] = {orbital1, orbital2 };
    unsigned short Spins[2] = { spin1, spin2 };
    T->OperatorSequence.assign(OperatorSequence,OperatorSequence+2);
    T->SiteLabels.assign(Labels,Labels+2);
    T->Orbitals.assign(Orbitals,Orbitals+2);
    T->Spins.assign(Spins,Spins+2);
    T->Value = Value;
    return T;
};

Lattice::Term* Lattice::Term::Presets::Hopping ( const std::string& Label1, const std::string& Label2, MelemType Value, unsigned short orbital, unsigned short spin)
{
    return Hopping(Label1, Label2, Value, orbital, orbital, spin, spin);
}

Lattice::Term* Lattice::Term::Presets::Level ( const std::string& Label, MelemType Value, unsigned short orbital, unsigned short spin)
{
    Lattice::Term *T = new Term(2);
    bool OperatorSequence[2] = { 1, 0 };
    std::string Labels[2] = { Label, Label };
    unsigned short Orbitals[2] = {orbital, orbital };
    unsigned short Spins[2] = { spin, spin };
    T->OperatorSequence.assign(OperatorSequence,OperatorSequence+2);
    T->SiteLabels.assign(Labels,Labels+2);
    T->Orbitals.assign(Orbitals,Orbitals+2);
    T->Spins.assign(Spins,Spins+2);
    T->Value = Value;
    return T;
};

// 4-index presets

Lattice::Term* Lattice::Term::Presets::NupNdown ( const std::string& Label1, const std::string& Label2,  MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2 )
{
    if (Label1 == Label2 && spin1 == spin2 && orbital1 == orbital2) {
        ERROR("NupNdown terms should have different sites or spin or orbital components. Returning simple Level term.");
        return Lattice::Term::Presets::Level(Label1, Value, orbital1, spin1);
        };
    Lattice::Term *T = new Term(4);
    bool order[4]              =      { 1,     0,     1,     0     };
    unsigned short Spins[4]    =      { spin1, spin1, spin2, spin2  };
    std::string Labels[4]      =      { Label1, Label1, Label2, Label2 };
    unsigned short Orbitals[4] =      { orbital1, orbital1, orbital2, orbital2 };
    T->OperatorSequence.assign(order,order+4);
    T->SiteLabels.assign(Labels,Labels+4);
    T->Orbitals.assign(Orbitals,Orbitals+4);
    T->Spins.assign(Spins,Spins+4);
    T->Value = Value;
    return T;
}

Lattice::Term* Lattice::Term::Presets::NupNdown ( const std::string& Label, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2 )
{
    return NupNdown(Label, Label, Value, orbital1, orbital2, spin1, spin2);
}

Lattice::Term* Lattice::Term::Presets::NupNdown ( const std::string& Label, MelemType Value, unsigned short orbital1, unsigned short orbital2 )
{
    return NupNdown(Label, Label, Value, orbital1, orbital2, up, down);
};

Lattice::Term* Lattice::Term::Presets::NupNdown ( const std::string& Label, MelemType Value, unsigned short Orbital, unsigned short spin1, unsigned short spin2 )
{
    return NupNdown(Label, Label, Value, Orbital, Orbital, spin1, spin2);
};

Lattice::Term* Lattice::Term::Presets::Spinflip ( const std::string& Label, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2 )
{
    if (orbital1 == orbital2 || spin1 == spin2) { ERROR("Spinflips should have different spins and orbitals"); throw (exWrongIndices()); }
    Lattice::Term *T = new Term(4);
    bool order[4]              =      { 1,     1,     0,     0 };
    unsigned short Spins[4]    =      { spin1, spin2, spin1, spin2  };
    std::string Labels[4]      =      { Label, Label, Label, Label };
    unsigned short Orbitals[4] =      { orbital1, orbital2, orbital2, orbital1 };
    T->OperatorSequence.assign(order,order+4);
    T->Spins.assign(Spins,Spins+4);
    T->SiteLabels.assign(Labels,Labels+4);
    T->Orbitals.assign(Orbitals,Orbitals+4);
    T->Value = Value;
    return T;
}

Lattice::Term* Lattice::Term::Presets::PairHopping ( const std::string& Label, MelemType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2 )
{
    if (orbital1 == orbital2 || spin1 == spin2) { ERROR("Pair hopping terms should have different spins and orbitals"); throw (exWrongIndices()); }
    Lattice::Term *T = new Term(4);
    bool order[4]              =      { 1,     1,     0,     0 };
    unsigned short Spins[4]    =      { spin1, spin2, spin1, spin2  };
    std::string Labels[4]      =      { Label, Label, Label, Label };
    unsigned short Orbitals[4] =      { orbital1, orbital1, orbital2, orbital2 };
    T->OperatorSequence.assign(order,order+4);
    T->Spins.assign(Spins,Spins+4);
    T->SiteLabels.assign(Labels,Labels+4);
    T->Orbitals.assign(Orbitals,Orbitals+4);
    T->Value = Value;
    return T;
}

Lattice::Term* Lattice::Term::Presets::SplusSminus ( const std::string& Label1, const std::string& Label2, MelemType Value, unsigned short orbital)
{
    Term *T = new Term(4);
    bool order[4]              =      { 1,     0,     1,     0 };
    unsigned short Spins[4]    =      { up, down, down, up  };
    std::string Labels[4]      =      { Label1, Label1, Label2, Label2 };
    unsigned short Orbitals[4] =      { orbital, orbital, orbital, orbital };
    T->OperatorSequence.assign(order,order+4);
    T->Spins.assign(Spins,Spins+4);
    T->SiteLabels.assign(Labels,Labels+4);
    T->Orbitals.assign(Orbitals,Orbitals+4);
    T->Value = Value;
    return T;
}

Lattice::Term* Lattice::Term::Presets::SminusSplus ( const std::string& Label1, const std::string& Label2, MelemType Value, unsigned short orbital)
{
    Term *T = SplusSminus(Label1, Label2, Value, orbital);
    unsigned short Spins[4]    =      { down, up, up, down  };
    T->Spins.assign(Spins,Spins+4);
    T->Value = Value;
    return T;
}


const char* Lattice::Term::Presets::exWrongIndices::what() const throw(){
    return "Wrong indices for selected term";
};

//
// LatticePresets
//

void LatticePresets::addCoulombS(Lattice *L, const std::string& Label, MelemType U, MelemType Level)
{
    if (L->Sites.find(Label)==L->Sites.end()) throw (Lattice::exWrongLabel());
    unsigned short Orbitals = L->Sites[Label]->OrbitalSize;
    unsigned short Spins = L->Sites[Label]->SpinSize;
    for (unsigned short i=0; i<Orbitals; ++i)
        for (unsigned short z1=0; z1<Spins; ++z1){
            if (std::abs(Level)) L->Terms->addTerm(Lattice::Term::Presets::Level(Label, Level, i, z1 ));
            for (unsigned short z2=0; z2<z1; ++z2) {
                if (std::abs(U)) L->Terms->addTerm(Lattice::Term::Presets::NupNdown(Label, U, i, i, z1, z2));
                };
       };
};

void LatticePresets::addCoulombP(Lattice *L, const std::string& Label, MelemType U, MelemType U_p, MelemType J, MelemType Level)
{
    if (L->Sites.find(Label)==L->Sites.end()) throw (Lattice::exWrongLabel());
    unsigned short Orbitals = L->Sites[Label]->OrbitalSize;
    unsigned short Spins = L->Sites[Label]->SpinSize;
    if (Orbitals<=1 || Spins<=1) { ERROR("Cannot add multiorbital interaction to a site with 1 orbital or 1 spin"); throw (Lattice::Term::Presets::exWrongIndices()); };
    for (unsigned short i=0; i<Orbitals; ++i){
        for (unsigned short z1=0; z1<Spins; ++z1){
            if (std::abs(Level)) L->Terms->addTerm(Lattice::Term::Presets::Level(Label, Level, i, z1 ));
            for (unsigned short j=0; j<Orbitals; ++j) if (i!=j) L->Terms->addTerm(Lattice::Term::Presets::NupNdown(Label, (U_p-J)/2., i, j, z1, z1));
            for (unsigned short z2=0; z2<z1; ++z2){
                if (std::abs(U)) L->Terms->addTerm(Lattice::Term::Presets::NupNdown(Label, U, i, i, z1, z2));
                for (unsigned short j=0; j<Orbitals; ++j){
                        if (i!=j) {
                            if (std::abs(U_p)) L->Terms->addTerm(Lattice::Term::Presets::NupNdown(Label, U_p, i, j, z1, z2));
                            if (std::abs(J)) {
                                L->Terms->addTerm(Lattice::Term::Presets::Spinflip(Label, -J, i, j, z1, z2));
                                L->Terms->addTerm(Lattice::Term::Presets::PairHopping(Label, -J, i, j, z1, z2));
                                };
                            };
                        };
                    };
            };
        };
};

void LatticePresets::addCoulombP(Lattice *L, const std::string& Label, MelemType U, MelemType J, MelemType Level)
{
    addCoulombP(L, Label, U, U-2.0*J, J, Level);
}

void LatticePresets::addLevel(Lattice *L, const std::string& Label, MelemType Level)
{
    if (L->Sites.find(Label)==L->Sites.end()) throw (Lattice::exWrongLabel());
    unsigned short Orbitals = L->Sites[Label]->OrbitalSize;
    unsigned short Spins = L->Sites[Label]->SpinSize;
    for (unsigned short i=0; i<Orbitals; ++i)
        for (unsigned short z=0; z<Spins; ++z) {
            if (std::abs(Level)) L->Terms->addTerm(Lattice::Term::Presets::Level(Label, Level, i, z ));
            };
}

void LatticePresets::addMagnetization(Lattice *L, const std::string& Label, MelemType Magnetization)
{
    if (L->Sites.find(Label)==L->Sites.end()) { ERROR("No site" << Label << "found."); throw (Lattice::exWrongLabel()); };
    unsigned short Orbitals = L->Sites[Label]->OrbitalSize;
    unsigned short Spins = L->Sites[Label]->SpinSize;
    if (Spins!=2) { ERROR("addSzSz doesn't work not for 2 spins."); throw (Lattice::Term::Presets::exWrongIndices()); };
    for (unsigned short i=0; i<Orbitals; ++i) {
            L->Terms->addTerm(Lattice::Term::Presets::Level(Label, Magnetization, i, Pomerol::up ));
            L->Terms->addTerm(Lattice::Term::Presets::Level(Label, -Magnetization, i, Pomerol::down ));
        };
}

void LatticePresets::addSzSz ( Lattice *L, const std::string& Label1, const std::string& Label2, MelemType ExchJ)
{

    if (L->Sites.find(Label1)==L->Sites.end()) { ERROR("No site" << Label1 << "found."); throw (Lattice::exWrongLabel()); };
    if (L->Sites.find(Label2)==L->Sites.end()) { ERROR("No site" << Label2 << "found."); throw (Lattice::exWrongLabel()); };
    unsigned short Orbitals = L->Sites[Label1]->OrbitalSize;
    unsigned short Spins = L->Sites[Label1]->SpinSize;
    if ( Orbitals != L->Sites[Label2]->OrbitalSize || Spins != L->Sites[Label1]->SpinSize ) { ERROR("Adjacent sites spin and orbital sizes do not match. Use Lattice::addTerm for specific interaction."); throw (Lattice::Term::Presets::exWrongIndices()); };
    if (Spins!=2) { ERROR("addSzSz doesn't work not for 2 spins."); throw (Lattice::exWrongLabel()); };
    for (unsigned short i=0; i<Orbitals; ++i) {
            L->Terms->addTerm(Lattice::Term::Presets::NupNdown(Label1, Label2, -ExchJ/4., i, i, up, down));
            L->Terms->addTerm(Lattice::Term::Presets::NupNdown(Label1, Label2, -ExchJ/4., i, i, down, up));
            if ( Label1 != Label2) {
                L->Terms->addTerm(Lattice::Term::Presets::NupNdown(Label1, Label2, ExchJ/4., i, i, up, up));
                L->Terms->addTerm(Lattice::Term::Presets::NupNdown(Label1, Label2, ExchJ/4., i, i, down, down));
                }
            else {
                L->Terms->addTerm(Lattice::Term::Presets::Level(Label1, ExchJ/4., i, up));
                L->Terms->addTerm(Lattice::Term::Presets::Level(Label1, ExchJ/4., i, down));
                };
        };
}

void LatticePresets::addSS ( Lattice *L, const std::string& Label1, const std::string& Label2, MelemType ExchJ)
{
    if (L->Sites.find(Label1)==L->Sites.end()) { ERROR("No site" << Label1 << "found."); throw (Lattice::exWrongLabel()); };
    if (L->Sites.find(Label2)==L->Sites.end()) { ERROR("No site" << Label2 << "found."); throw (Lattice::exWrongLabel()); };
    unsigned short Orbitals = L->Sites[Label1]->OrbitalSize;
    unsigned short Spins = L->Sites[Label1]->SpinSize;
    if ( Orbitals != L->Sites[Label2]->OrbitalSize || Spins != L->Sites[Label1]->SpinSize ) { ERROR("Adjacent sites spin and orbital sizes do not match. Use Lattice::addTerm for specific interaction."); throw (Lattice::Term::Presets::exWrongIndices()); };
    if (Spins!=2) { ERROR("addSS doesn't work not for 2 spins."); throw (Lattice::exWrongLabel()); };

    addSzSz(L, Label1, Label2, ExchJ);
    for (unsigned short i=0; i<Orbitals; ++i) {
        L->Terms->addTerm(Lattice::Term::Presets::SplusSminus(Label1, Label2, ExchJ/2., i));
        L->Terms->addTerm(Lattice::Term::Presets::SminusSplus(Label1, Label2, ExchJ/2., i));
        };
}

void LatticePresets::addHopping ( Lattice *L, const std::string& Label1, const std::string& Label2, MelemType t, unsigned short Orbital1, unsigned short Orbital2, unsigned short Spin1, unsigned short Spin2)
{
    if (L->Sites.find(Label1)==L->Sites.end()) { ERROR("No site" << Label1 << "found."); throw (Lattice::exWrongLabel()); };
    if (L->Sites.find(Label2)==L->Sites.end()) { ERROR("No site" << Label2 << "found."); throw (Lattice::exWrongLabel()); };
    if (Orbital1 >= L->Sites[Label1]->OrbitalSize || Orbital2 >= L->Sites[Label2]->OrbitalSize || Spin1 >= L->Sites[Label1]->SpinSize || Spin2 >= L->Sites[Label2]->SpinSize ) {
        ERROR("Orbital or Spin index mismatch"); throw ( Lattice::Term::Presets::exWrongIndices() );
        };
    L->addTerm(Lattice::Term::Presets::Hopping(Label1, Label2, t, Orbital1, Orbital2, Spin1, Spin2));
    #ifdef POMEROL_COMPLEX_MATRIX_ELEMENTS
    L->addTerm(Lattice::Term::Presets::Hopping(Label2, Label1, conj(t), Orbital2, Orbital1, Spin2, Spin1)); // Hermite conjugate
    #else
    L->addTerm(Lattice::Term::Presets::Hopping(Label2, Label1, t, Orbital2, Orbital1, Spin2, Spin1));
    #endif
}

void LatticePresets::addHopping ( Lattice *L, const std::string& Label1, const std::string& Label2, MelemType t, unsigned short Orbital1, unsigned short Orbital2, unsigned short Spin)
{
    addHopping(L,Label1, Label2, t, Orbital1, Orbital2, Spin, Spin);
}

void LatticePresets::addHopping ( Lattice *L, const std::string& Label1, const std::string& Label2, MelemType t, unsigned short Orbital1, unsigned short Orbital2)
{
    if (L->Sites.find(Label1)==L->Sites.end()) { ERROR("No site " << Label1 << " found."); throw (Lattice::exWrongLabel()); };
    if (L->Sites.find(Label2)==L->Sites.end()) { ERROR("No site " << Label2 << " found."); throw (Lattice::exWrongLabel()); };
    if (Orbital1 >= L->Sites[Label1]->OrbitalSize || Orbital2 >= L->Sites[Label2]->OrbitalSize ){
        ERROR("Orbital or Spin index mismatch"); throw ( Lattice::Term::Presets::exWrongIndices() );
        };
    unsigned short Spins = L->Sites[Label1]->SpinSize;
    if ( Spins != L->Sites[Label1]->SpinSize ) { ERROR("Adjacent sites spin sizes do not match. Use Lattice::addTerm for specific interaction."); throw (Lattice::Term::Presets::exWrongIndices()); };
    for (int z=0; z<Spins; ++z) addHopping(L, Label1, Label2, t, Orbital1, Orbital2, z, z);
}

void LatticePresets::addHopping ( Lattice *L, const std::string& Label1, const std::string& Label2, MelemType t)
{
    if (L->Sites.find(Label1)==L->Sites.end()) { ERROR("No site" << Label1 << "found."); throw (Lattice::exWrongLabel()); };
    if (L->Sites.find(Label2)==L->Sites.end()) { ERROR("No site" << Label2 << "found."); throw (Lattice::exWrongLabel()); };
    unsigned short Orbitals = L->Sites[Label1]->OrbitalSize;
    unsigned short Spins = L->Sites[Label1]->SpinSize;
    if ( Orbitals != L->Sites[Label2]->OrbitalSize || Spins != L->Sites[Label1]->SpinSize ) {
        ERROR("Adjacent sites spin and orbital sizes do not match. Use Lattice::addTerm for specific interaction."); throw (Lattice::Term::Presets::exWrongIndices());
        };

    for (unsigned short z=0; z<Spins; ++z)
        for (unsigned short i=0; i<Orbitals; ++i)
            addHopping(L, Label1, Label2, t, i, i, z, z);

}

} // end of namespace Pomerol
