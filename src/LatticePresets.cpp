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

/** \file src/Lattice.h
** \brief A lattice handler. Reads and stores the lattice from a JSON file. Designed to be extended to other formats.
** 
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/

#include "LatticePresets.h"

namespace Pomerol {

//
// Lattice::Term::Presets
//

// 2-index presets

Lattice::Term* Lattice::Term::Presets::Hopping ( const std::string& Label1, const std::string& Label2, RealType Value, unsigned short orbital, unsigned short spin)
{
    Lattice::Term *T = new Term(2);
    bool Order[2] = { 1, 0 };
    std::string labels[2] = { Label1, Label2 };
    unsigned short Orbitals[2] = {orbital, orbital };
    unsigned short Spins[2] = { spin, spin };
    T->Order.assign(Order,Order+2);
    T->SiteLabels.assign(labels,labels+2); 
    T->Orbitals.assign(Orbitals,Orbitals+2);
    T->Spins.assign(Spins,Spins+2); 
    T->Value = Value;
    return T;
};

Lattice::Term* Lattice::Term::Presets::Level ( const std::string& Label, RealType Value, unsigned short orbital, unsigned short spin)
{
    return Hopping(Label, Label, Value, orbital, spin);
};

// 4-index presets

Lattice::Term* Lattice::Term::Presets::NupNdown ( const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2 )
{
    Lattice::Term *T = new Term(4); 
    bool order[4]              =      { 1,     0,     1,     0     };       
    unsigned short Spins[4]    =      { spin1, spin1, spin2, spin2  };
    std::string labels[4]      =      { Label, Label, Label, Label };
    unsigned short Orbitals[4] =      { orbital1, orbital1, orbital2, orbital2 }; 
    T->Order.assign(order,order+4);
    T->SiteLabels.assign(labels,labels+4); 
    T->Orbitals.assign(Orbitals,Orbitals+4);
    T->Spins.assign(Spins,Spins+4); 
    T->Value = Value;
    return T;
}

Lattice::Term* Lattice::Term::Presets::NupNdown ( const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2 )
{ 
    return NupNdown(Label, Value, orbital1, orbital2, up, down);
};

Lattice::Term* Lattice::Term::Presets::NupNdown ( const std::string& Label, RealType Value, unsigned short Orbital, unsigned short spin1, unsigned short spin2 )
{
    return NupNdown(Label, Value, Orbital, Orbital, spin1, spin2);
};

Lattice::Term* Lattice::Term::Presets::Spinflip ( const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2 )
{
    if (orbital1 == orbital2 || spin1 == spin2) { ERROR("Spinflips should have different spins and orbitals"); }
    Lattice::Term *T = new Term(4); 
    bool order[4]              =      { 1,     1,     0,     0 };       
    unsigned short Spins[4]    =      { spin1, spin2, spin1, spin2  };
    std::string labels[4]      =      { Label, Label, Label, Label };
    unsigned short Orbitals[4] =      { orbital1, orbital2, orbital2, orbital1 }; 
    T->Order.assign(order,order+4);
    T->Spins.assign(Spins,Spins+4); 
    T->SiteLabels.assign(labels,labels+4); 
    T->Orbitals.assign(Orbitals,Orbitals+4);
    T->Value = Value;
    return T;
}

Lattice::Term* Lattice::Term::Presets::PairHopping ( const std::string& Label, RealType Value, unsigned short orbital1, unsigned short orbital2, unsigned short spin1, unsigned short spin2 )
{
    if (orbital1 == orbital2 || spin1 == spin2) { ERROR("Pair hopping terms should have different spins and orbitals"); }
    Lattice::Term *T = new Term(4); 
    bool order[4]              =      { 1,     1,     0,     0 };       
    unsigned short Spins[4]    =      { spin1, spin2, spin1, spin2  };
    std::string labels[4]      =      { Label, Label, Label, Label };
    unsigned short Orbitals[4] =      { orbital1, orbital1, orbital2, orbital2 }; 
    T->Order.assign(order,order+4);
    T->Spins.assign(Spins,Spins+4); 
    T->SiteLabels.assign(labels,labels+4); 
    T->Orbitals.assign(Orbitals,Orbitals+4);
    T->Value = Value;
    return T;
}


//
// Lattice::Presets
//

void Lattice::Presets::addSSite(Lattice *L, const std::string& label, RealType U, RealType Level, unsigned short Orbitals, unsigned short Spins)
{
    Lattice::Site* current = new Lattice::Site(label, Orbitals, Spins);
    current->label=label;
    L->Sites[label]=current;
    Log.setDebugging(true);
    for (unsigned short i=0; i<Orbitals; ++i) 
        for (unsigned short z1=0; z1<Spins; ++z1){
            L->Terms->addTerm(Lattice::Term::Presets::Level(label, Level, i, z1 ));
            for (unsigned short z2=0; z2<z1; ++z2) {
                L->Terms->addTerm(Lattice::Term::Presets::NupNdown(label, U, i, i, z1, z2));
                /*if (Magnetization) { 
                    Term *T = Lattice::Term::Presets::Level(label,  Magnetization, i, z1); 
                    DEBUG("Adding term " << *T);
                    L->Terms->addTerm(T);
                    delete T;
                    T = Lattice::Term::Presets::Level(label, -Magnetization, i, z2); 
                    DEBUG("Adding term " << *T);
                    L->Terms->addTerm(T);
                    delete T;
                    };
                */
            };
       };
};

void Lattice::Presets::addPSite(Lattice *L, const std::string& label, RealType U, RealType U_p, RealType J, RealType Level, unsigned short Orbitals, unsigned short Spins)
{
    Lattice::Site* current = new Lattice::Site(label, Orbitals, Spins);
    L->Sites[label]=current;
    Log.setDebugging(true);
    for (unsigned short i=0; i<Orbitals; ++i){
        for (unsigned short z1=0; z1<Spins; ++z1){
            L->Terms->addTerm(Lattice::Term::Presets::Level(label, Level, i, z1 ));
            for (unsigned short j=0; j<Orbitals; ++j) if (i!=j) L->Terms->addTerm(Lattice::Term::Presets::NupNdown(label, (U_p-J)/2., i, j, z1, z1));
            for (unsigned short z2=0; z2<z1; ++z2){
                L->Terms->addTerm(Lattice::Term::Presets::NupNdown(label, U, i, i, z1, z2));
                for (unsigned short j=0; j<Orbitals; ++j){
                        if (i!=j) {
                            L->Terms->addTerm(Lattice::Term::Presets::NupNdown(label, U_p, i, j, z1, z2));
                            L->Terms->addTerm(Lattice::Term::Presets::Spinflip(label, -J, i, j, z1, z2)); 
                            L->Terms->addTerm(Lattice::Term::Presets::PairHopping(label, -J, i, j, z1, z2)); 
                            };
                        };
                    };
            };
        };
};

void Lattice::Presets::addPSite(Lattice *L, const std::string& label, RealType U, RealType J, RealType Level, unsigned short Orbitals, unsigned short Spins)
{
    addPSite(L, label, U, U-2.0*J, J, Level, Orbitals, Spins);
}

void Lattice::Presets::addPSite(Lattice *L, const std::string& label, RealType U, RealType J, RealType Level, unsigned Orbitals)
{
    addPSite(L, label, U, U-2.0*J, J, Level, Orbitals, 2);
}

void Lattice::Presets::addMagnetization(Lattice *L, const std::string& label, RealType Magnetization, unsigned short Orbitals, unsigned short Spins)
{
}


} // end of namespace Pomerol
