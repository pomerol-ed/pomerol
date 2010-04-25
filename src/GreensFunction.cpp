#include "GreensFunction.h"

extern IniConfig* pIni;

GreensFunction::GreensFunction(StatesClassification& S, Hamiltonian& H, 
                               AnnihilationOperator& C, CreationOperator& CX, DensityMatrix& DM,
                               output_handle &OUT) : S(S)
{
    green_path = output_handle(OUT.path() + "/Green_func");   
    NumOfBlocks = S.NumberOfBlocks();
    
    parts = new GreensFunctionPart* [NumOfBlocks];
    for (BlockNumber current_block=0;current_block<NumOfBlocks;current_block++)
    {
      parts[current_block] = new GreensFunctionPart((AnnihilationOperatorPart&)C.getPartFromRightIndex(current_block),
                                                    (CreationOperatorPart&)CX.getPartFromRightIndex(current_block),
                                                    H.part(current_block), DM.part(current_block)); 
#warning - check parts
    }
}

GreensFunction::~GreensFunction()
{
    delete[] parts;
}

ComplexType GreensFunction::operator()(ComplexType Frequency)
{     
      ComplexType Value = 0;
      for (BlockNumber current_block=0;current_block<NumOfBlocks;current_block++)
          Value += (*parts[current_block])(Frequency);
      return Value;
}

string GreensFunction::getPath()
{
    return green_path.fullpath();
}

/*using std::stringstream;


void green::building()							//building Green Function
{
	
	IniConfig Ini("system.ini");					//parameters of building Green Function

	int beta = Ini["Green Function:beta"];			//inverse temperature	
	int matsubara = Ini["Green Function:matsubara"];		//w = iw or w = w' + iw"
	
	//numeric constant:

	//begining of calculating statistic summa

	RealType Stat_sum = 0;
	for (QuantumState m =0; m<S.N_st(); m++)
		Stat_sum += exp(-beta*H.eigenval(m));
	
	cout << "Partition function = " << Stat_sum << endl;
	
	//finishing of calculating
	
	int size1 = operatorsCiCXj.reVecC().size();			//size of C
	int size2 = operatorsCiCXj.reVecCX().size();			//size of CX

	if (matsubara == 1)
	{
		int points = Ini["Green Function:points"];		//number of points Green Function
		
		G_ij = new ComplexType [points];			//creation of massive GF
		
		for(int k=0; k<points; k++)
		{
		progressbar(int(k*1.0/points*100));
			G_ij[k] = 0.;
		
			if(i==j)
			{	
				for (int p=0; p<size2; p++)
				{	
					
					QuantumState n = operatorsCiCXj.reVecCX()[p].n;
					QuantumState m = operatorsCiCXj.reVecCX()[p].m;
					RealType CX_nm = operatorsCiCXj.reVecCX()[p].C;

					RealType Up = (1./Stat_sum)*CX_nm*CX_nm*( exp(-beta*H.eigenval(n)) + exp(-beta*H.eigenval(m)) );
						
					ComplexType Down = I*(Pi*(2*k + 1)/beta) + ( H.eigenval(n)-H.eigenval(m) );
					
					G_ij[k] += Up/Down;
				}
			}
			
			else
			{
				for (unsigned int p1=0; p1<operatorsCiCXj.reVecC().size(); p1++)
				{
					for (unsigned int p2=0; p2<operatorsCiCXj.reVecCX().size(); p2++)
					{
			
						QuantumState n = operatorsCiCXj.reVecCX()[p2].n;
						QuantumState m = operatorsCiCXj.reVecCX()[p2].m;

						int flag=0;
			
						if ( (operatorsCiCXj.reVecC()[p1].n == m) && (operatorsCiCXj.reVecC()[p1].m == n)  )
						{
			
							RealType CX_nm = operatorsCiCXj.reVecCX()[p2].C;
							RealType C_mn = operatorsCiCXj.reVecC()[p1].C;
	
							RealType Up = (1./Stat_sum)*CX_nm*C_mn*( exp(-beta*H.eigenval(n)) + exp(-beta*H.eigenval(m)) );
							
							ComplexType Down = I*(Pi*(2*k + 1)/beta) + ( H.eigenval(n)-H.eigenval(m) );
					
							G_ij[k] += Up/Down;
							
				//			operatorsCiCXj.reVecC().erase(operatorsCiCXj.reVecC().begin()+ p1);
				//			operatorsCiCXj.reVecCX().erase(operatorsCiCXj.reVecCX().begin() + p2);
											
							flag=1;
						}

						if (flag) break;
					}
				}
			}
		}
	}

	if (matsubara == 0)
	{
		RealType Imw = Ini["Green Function:Imw"];			//(matsubara=0) value of w"
		int w_min = Ini["Green Function:w_min"];		//minimum value of w'
		int w_max = Ini["Green Function:w_max"];		//maximum value of w'
		RealType step = Ini["Green Function:step"];		//step delta w'

		int Npoints = int((w_max-w_min)/step  + 1);			//number of points Green Function

		G_ij = new ComplexType [Npoints];			//creation of massive GF
		int old_percent = 0;
		
		for(int k=0; k<Npoints; k++)
		{
			G_ij[k] = 0.;
			int current_percent = (int) (k*100.0/Npoints);
			if (current_percent != old_percent) {progressbar(current_percent); old_percent = current_percent;};

			if(i==j)
			{	
				for (int p=0; p<size2; p++)
				{	
					
					QuantumState n = operatorsCiCXj.reVecCX()[p].n;
					QuantumState m = operatorsCiCXj.reVecCX()[p].m;
					RealType CX_nm = operatorsCiCXj.reVecCX()[p].C;

					RealType Up = (1./Stat_sum)*CX_nm*CX_nm*( exp(-beta*H.eigenval(n)) + exp(-beta*H.eigenval(m)) );
						
					ComplexType Down = I*Imw + ( (w_min + k*step) + (H.eigenval(n) - H.eigenval(m)) );
					
					G_ij[k] += Up/Down;
				}
			}
			
			else
			{
				for (int p1=0; p1<size1; p1++)
				{
					for (int p2=0; p2<size2; p2++)
					{
			
						QuantumState n = operatorsCiCXj.reVecCX()[p2].n;
						QuantumState m = operatorsCiCXj.reVecCX()[p2].m;
						
						int flag = 0;

						if ( (operatorsCiCXj.reVecC()[p1].n == m) && (operatorsCiCXj.reVecC()[p1].m == n)  )
						{
			
							RealType CX_nm = operatorsCiCXj.reVecCX()[p2].C;
							RealType C_mn = operatorsCiCXj.reVecC()[p1].C;

							RealType Up = (1./Stat_sum)*CX_nm*C_mn*( exp(-beta*H.eigenval(n)) + exp(-beta*H.eigenval(m)) );
						
							ComplexType Down = I*Imw + ( (w_min + k*step) + (H.eigenval(n)-H.eigenval(m)) );
					
							G_ij[k] += Up/Down;

							flag=1;
						}

						if (flag) break;
					}
				}
			}
		}
	}
}

*/

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
