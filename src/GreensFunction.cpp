#include "GreensFunction.h"

extern IniConfig* pIni;

GreensFunction::GreensFunction(StatesClassification& S, Hamiltonian& H, 
                               AnnihilationOperator& C, CreationOperator& CX, DensityMatrix& DM,
                               output_handle &OUT) : parts(0), S(S), H(H), C(C), CX(CX), DM(DM)
{
    green_path = output_handle(OUT.path() + "/GreensFunction");   
}

GreensFunction::~GreensFunction()
{
      for(std::list<GreensFunctionPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
          delete *iter;
}

void GreensFunction::prepare(void)
{
    std::list<BlockMapping> CNontrivialBlocks = C.getNonTrivialIndices();
    std::list<BlockMapping> CXNontrivialBlocks = CX.getNonTrivialIndices();
    
    std::list<BlockMapping>::const_iterator Citer = CNontrivialBlocks.begin();
    std::list<BlockMapping>::const_iterator CXiter = CXNontrivialBlocks.begin();
    
    RealType GroundEnergy = H.getGroundEnergy();
    
    while(Citer != CNontrivialBlocks.end() && CXiter != CXNontrivialBlocks.end()){
        BlockNumber Cleft = Citer->first;
        BlockNumber Cright = Citer->second;
        BlockNumber CXleft = CXiter->first;
        BlockNumber CXright = CXiter->second;
  
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
}

void GreensFunction::compute(void)
{
      for(std::list<GreensFunctionPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
          (*iter)->compute();
}

ComplexType GreensFunction::operator()(ComplexType Frequency)
{     
      ComplexType Value = 0;
      for(std::list<GreensFunctionPart*>::iterator iter = parts.begin(); iter != parts.end(); iter++)
          Value += (**iter)(Frequency);
      return Value;
}

string GreensFunction::getPath()
{
    return green_path.fullpath();
}


//other functions

void GreensFunction::dumpMatsubara(unsigned short points)
{
	std::stringstream filename;
	unsigned short i=C.getBit();
	unsigned short j=CX.getBit();
	filename << green_path.path() << "/Gw" << i << j << ".dat";
	std::ofstream output;
	output.open(filename.str().c_str());
	
	for (int k=0; k<points; k++)
	{	
		ComplexType iw = M_PI/DM.getBeta()*(2*k+1.)*I;
        ComplexType G = (*this)(iw);
		output << k << "\t" << std::setprecision(9) << real(G) << "\t" << imag(G) << endl;
	}
		
	output.close();
}

