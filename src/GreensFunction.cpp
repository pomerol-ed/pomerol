/** \file src/GreensFunction.cpp
** \brief Thermal Green's function.
**
** \author Igor Krivenko (igor@shg.ru)
** \author Andrey Antipov (antipov@ct-qmc.org)
*/
#include "GreensFunction.h"

GreensFunction::GreensFunction(StatesClassification& S, Hamiltonian& H, 
                               AnnihilationOperator& C, CreationOperator& CX, DensityMatrix& DM,
                               output_handle &OUT) : S(S), H(H), C(C), CX(CX), DM(DM), parts(0)
{
    green_path = output_handle(OUT.path() + "/GreensFunction");   
    vanish = true;
}

GreensFunction::~GreensFunction()
{
    for(std::list<GreensFunctionPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
        delete *iter;
}

void GreensFunction::prepare(void)
{
    // Find out non-trivial blocks of C and CX.
    std::list<BlockMapping> CNontrivialBlocks = C.getNonTrivialIndices();
    std::list<BlockMapping> CXNontrivialBlocks = CX.getNonTrivialIndices();

    std::list<BlockMapping>::const_iterator Citer = CNontrivialBlocks.begin();
    std::list<BlockMapping>::const_iterator CXiter = CXNontrivialBlocks.begin();

    while(Citer != CNontrivialBlocks.end() && CXiter != CXNontrivialBlocks.end()){
        // <Cleft|C|Cright><CXleft|CX|CXright>
        BlockNumber Cleft = Citer->first;
        BlockNumber Cright = Citer->second;
        BlockNumber CXleft = CXiter->first;
        BlockNumber CXright = CXiter->second;

        // Select a relevant 'world stripe' (sequence of blocks).
        if(Cleft == CXright && Cright == CXleft){
              parts.push_back(new GreensFunctionPart(
                              (AnnihilationOperatorPart&)C.getPartFromLeftIndex(Cleft),
                              (CreationOperatorPart&)CX.getPartFromRightIndex(CXright),
                              H.part(Cright), H.part(Cleft),
                              DM.part(Cright), DM.part(Cleft)));
        }

        unsigned long CleftInt = Cleft;
        unsigned long CXrightInt = CXright;

        if(CleftInt <= CXrightInt) Citer++;
        if(CleftInt >= CXrightInt) CXiter++;
    }
if (parts.size() > 0) vanish = false;
}

void GreensFunction::compute(void)
{
      for(std::list<GreensFunctionPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
          (*iter)->compute();
}

ComplexType GreensFunction::operator()(long MatsubaraNum)
{
      ComplexType Value = 0;
      for(std::list<GreensFunctionPart*>::const_iterator iter = parts.begin(); iter != parts.end(); iter++)
          Value += (**iter)(MatsubaraNum);
      return Value;
}

unsigned short GreensFunction::getBit(size_t Position) const
{
    switch(Position){
        case 0: return C.getBit();
        case 1: return CX.getBit();
        default: assert(0);
    }
}

string GreensFunction::getPath()
{
    return green_path.fullpath();
}

bool GreensFunction::vanishes()
{
    return vanish;
}

//other functions

void GreensFunction::dumpMatsubara(unsigned short points)
{
    std::stringstream filename;
    unsigned short i=C.getBit();
    unsigned short j=CX.getBit();
    // TODO: Factor out this fragment to output_handle class.
    filename << green_path.path() << "/Gw" << i << j << ".dat";
    std::ofstream output;
    output.open(filename.str().c_str());

    for (int k=0; k<points; k++)
    {
        ComplexType G = (*this)(k);
        output << std::setprecision(9) << real(G) << "    " << imag(G) << "    " << endl;
    }

    output.close();
}

#warning"Dead code. Subject to remove."
//std::list<GreensFunctionPart::GreensTerm> GreensFunction::getTerms(void) const
// {
//     std::list<GreensFunctionPart::GreensTerm> allTerms;
//     for(std::list<GreensFunctionPart*>::const_iterator part = parts.begin(); part != parts.end(); part++){
//         std::list<GreensFunctionPart::GreensTerm> partTerms = (*part)->getTerms();
//         for(std::list<GreensFunctionPart::GreensTerm>::const_iterator term = partTerms.begin(); term != partTerms.end(); term++)
//             allTerms.push_back(*term);
//     }
//     return allTerms;
// }
