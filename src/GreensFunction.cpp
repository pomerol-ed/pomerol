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
    
    while(Citer != CNontrivialBlocks.end() && CXiter != CXNontrivialBlocks.end()){
        BlockNumber Cleft = Citer->first;
        BlockNumber Cright = Citer->second;
        BlockNumber CXleft = CXiter->first;
        BlockNumber CXright = CXiter->second;
  
        if(Cleft == CXright && Cright == CXleft){         
              parts.push_back(new GreensFunctionPart(
                              C.getPartFromLeftIndex(Cleft),
                              CX.getPartFromRightIndex(CXright),
                              H.part(Cright), H.part(Cleft), DM.part(Cright), DM.part(Cleft)));
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
/*
void GreensFunction::dump()
{
	bool matsubara = (*pIni)["Green Function:matsubara"];		//w = iw or w = w' + iw"
	int points = 0;

    RealType step;
    
	if (matsubara)
	{
		 points = (*pIni)["Green Function:points"];			//number of points Green Function 
	}
	else
	{
		int w_min = (*pIni)["Green Function:w_min"];		//minimum value of w'
		int w_max = (*pIni)["Green Function:w_max"];		//maximum value of w'
		step = (*pIni)["Green Function:step"];      		//step delta w'

		points = int((w_max-w_min)/step  + 1);		    	//number of points Green Function
	}

    int i = (*pIni)["Green Function:i"];
    int j = (*pIni)["Green Function:j"];
    
	std::stringstream filename;
	filename << green_path.path() << "/G(w)" << i << j << ".dat";
    ofstream output(filename.str());
	
	for (int k=0; k<points; k++)
		outHpart << k << "\t" << std::setprecision(9) << real(G_ij[k]) << "\t" << imag(G_ij[k]) << endl;
		
	outHpart.close();
}
*/
