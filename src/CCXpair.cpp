#include "CCXpair.h"
#include <sstream>
#include <fstream>

using std::stringstream;


//struct valC								//values rotated C or CX

valC::valC(int line, int column, RealType C_nm)				//inicialization valC
{	
	n=line;
	m=column;
	C=C_nm;
	
}
valC &valC::operator+=(const valC &rhs) 
{
	if (rhs.n == (*this).n && rhs.m == (*this).m)    
		(*this).C += rhs.C;
	return (*this);
}


//class matrixs								//rotates matrixes C and CX 

vector<valC>& matrixs::reVecC()						//return UXCU
{
	return UXCU;
}

vector<valC>& matrixs::reVecCX()					//return UXCXU(m,n)
{
	return UXCXU;
}

void matrixs::inimatrixs(int I, int J){
	i=I;
	j=J;
	
	matrixC_path = output_handle(OUT.path()+"/matrixC");
	matrixCX_path = output_handle(OUT.path()+"/matrixCX");
}

void matrixs::putMatrXC()
{
	//creation massive uxcu
	
	uxcu = new vector<valC> * [S.NumberOfBlocks()];
	for (int n=0; n<S.NumberOfBlocks(); n++)
	{
		uxcu[n] = new vector<valC> [S.NumberOfBlocks()];
		
		for (int m=0; m<S.NumberOfBlocks(); m++)
			uxcu[n][m] = UXCU;
	}
	
	//finish of creation

	//rotation of matrix C:

	for (QuantumState L=0; L<S.N_st(); L++)
	{
//		if (!(L%20)) cout << "Done " << int((L*1.0)/(1.0*S.N_st())*100.) << "%" << endl;
		progressbar(int((L*1.0)/(1.0*S.N_st())*100.));
		
		if (S.n_i(L,i)==0)
		{
			QuantumState K = retKforC(L);

			if( (mFuncC(L,K,i)!= 0) )
			{			
								
				int l=0, k=0;				// l,k in part of Hamilt			
				
				for ( unsigned int n=0; n<S.clstates(S.getStateInfo(L)).size(); n++ )
				{
				int flag=0;

					if (S.clstates(S.getStateInfo(L))[n] == L)
					{
					   l=n;				//get l
					   flag=1;
					}
					
					if(flag) break;
				}

				for ( unsigned int n=0; n<S.clstates(S.getStateInfo(K)).size(); n++ )
				{
					int flag=0;
					
					if (S.clstates(S.getStateInfo(K))[n] == K)
					{
						k=n;			//get k
						flag=1;
					}

					if(flag) break;
				}
				
				for ( unsigned int n=0; n<S.clstates(S.getStateInfo(L)).size(); n++)
				{
					if(H.block(S.getStateInfo(L)).reH(l,n)!=0)
					{
						for (unsigned int m=0; m<S.clstates(S.getStateInfo(K)).size(); m++)
						{
							int N = S.clstates(S.getStateInfo(L))[n];
							int M = S.clstates(S.getStateInfo(K))[m];
							RealType C_nm = H.block(S.getStateInfo(L)).reH(l,n)*mFuncC(L,K,i)*H.block(S.getStateInfo(K)).reH(k,m);
						
							int flag=0;
									
							if (C_nm!=0)
							{
								// indexes of uxcu:

								BlockNumber block_from = S.getBlockNumber(S.getStateInfo(N));
								BlockNumber block_to = S.getBlockNumber(S.getStateInfo(M));

								if (block_from.isCorrect() && block_to.isCorrect())
								{

								  for (unsigned int p=0; p<uxcu[block_from][block_to].size(); p++) 
								  {
								 	  if( (N==uxcu[block_from][block_to][p].n) && (M==uxcu[block_from][block_to][p].m))
									  {
									 	  uxcu[block_from][block_to][p].C+=C_nm; flag=1;
									  }
									  if (flag) break;
								  }
								  if (!flag)
								  {	
									  valC nCm(N,M,C_nm);
									  uxcu[block_from][block_to].push_back(nCm);
								  }
								}
							}
						}
					}
				}	
			}
		}	
	}
	
	for(int n=0; n<S.NumberOfBlocks(); n++)
	{
		for(int m=0; m<S.NumberOfBlocks(); m++)
		{
			for(unsigned int j=0; j<uxcu[n][m].size(); j++)
				UXCU.push_back(uxcu[n][m][j]);
		}
	}
	
}

void matrixs::putMatrXCX()
{
	//creation massive uxcxu
	
	uxcxu = new vector<valC> * [S.NumberOfBlocks()];
	for (int n=0; n<S.NumberOfBlocks(); n++)
	{
		uxcxu[n] = new vector<valC> [S.NumberOfBlocks()];
		
		for (int m=0; m<S.NumberOfBlocks(); m++)
			uxcxu[n][m] = UXCXU;
	}

	RealSparseMatrixType uxcxu2(S.N_st(),S.N_st());
	//finish of creation
		
	//rotation of matrix CX:
	//

	int old_percent=0;
	for (QuantumState L=0; L<S.N_st(); L++)
	{
//		if (!(L%20)) cout << "Done " << int((L*1.0)/(1.0*S.N_st())*100.) << "%" << endl;

		int current_percent = int((L*1.0)/(1.0*S.N_st())*100.);
		if (current_percent!=old_percent) { progressbar(current_percent); old_percent = current_percent; };
		

		if (S.n_i(L,j)==1)
		{


			QuantumState K = retKforCX(L);


			if( (mFuncCX(L,K,j)!= 0) )
			{	
				int l=0, k=0;				// l,k in part of Hamilt			


				for ( unsigned int n=0; n<S.clstates(S.getStateInfo(L)).size(); n++ )
				{

					int flag=0;
					
					if (S.clstates(S.getStateInfo(L))[n] == L)
					{
						l=n;			//get l
						flag=1;
					}

					if(flag) break;
				}

				for ( unsigned int n=0; n<S.clstates(S.getStateInfo(K)).size(); n++ )
				{
					int flag=0;
					
					if (S.clstates(S.getStateInfo(K))[n] == K)
					{
						k=n;			//get k
						flag=1;
					}

					if(flag) break;
				}

				for ( unsigned int n=0; n<S.clstates(S.getStateInfo(L)).size(); n++)
				{
					if(H.block(S.getStateInfo(L)).reH(l,n)!=0)
					{
						for (unsigned int m=0; m<S.clstates(S.getStateInfo(K)).size(); m++)
						{

							int N = S.clstates(S.getStateInfo(L))[n];
							int M = S.clstates(S.getStateInfo(K))[m];
							RealType CX_nm = H.block(S.getStateInfo(L)).reH(l,n)*mFuncCX(L,K,j)*H.block(S.getStateInfo(K)).reH(k,m);

							if (CX_nm!=0)
							{

								// indexes of uxcxu:

								BlockNumber block_from = S.getBlockNumber(S.getStateInfo(N));
								BlockNumber block_to = S.getBlockNumber(S.getStateInfo(M));

							   	int flag =0;
								if (block_from.isCorrect() && block_to.isCorrect())
								{
								//  valC nCXm(N,M,CX_nm);
								//  cout << N << M << endl;
								  uxcxu2.coeffRef(N,M)+=CX_nm;
								/*						  
								  for (unsigned int p=0; p<uxcxu[block_from][block_to].size(); p++) 
								  {
								  	  if( (N==uxcxu[block_from][block_to][p].n) && (M==uxcxu[block_from][block_to][p].m))
									  {
								  		  uxcxu[block_from][block_to][p].C+=CX_nm; flag=1;
									  }
									  if (flag) break;
								  }	
								  if (!flag)
								  {	
									  valC nCXm(N,M,CX_nm);
									  uxcxu[block_from][block_to].push_back(nCXm);
								  }*/
								}  
							}
						}

					}
				}	
			}
		}
	}

	for (int N=0; N<uxcxu2.outerSize(); ++N)
	 for (RealSparseMatrixType::InnerIterator it(uxcxu2,N); it; ++it)
		 {
			        if (fabs(it.value())>MATRIX_ELEMENT_TOLERANCE)
				UXCXU.push_back(valC(it.row(),it.col(),it.value()));

		 };

/*
	for(int n=0; n<S.NumberOfBlocks(); n++)
	{
		for(int m=0; m<S.NumberOfBlocks(); m++)
		{
			for(unsigned int d=0; d<uxcxu[n][m].size(); d++)
				UXCXU.push_back(uxcxu[n][m][d]);
		}
	}*/

}

int matrixs::retKforC(int L)							//return K for C
{	

	return( L + (1<<i) );

}

int matrixs::retKforCX(int L)							//return K for CX
{	

	return( L - (1<<j) );

}

int matrixs::mFuncC(long int state1, long int state2, int i)
{
	int flag=1, p=0;
	for (int m=0; m<S.N_b(); m++)
	{
		if( m == i )
		{
			if ((S.n_i(state2,i)==0) || (S.n_i(state1,i)==1) ) {flag=0;break;} else flag=1;
		}
		else
		{
			if (S.n_i(state1,m)==S.n_i(state2,m) ) flag=1; else {flag=0;break;}
		}
	}
	if (flag==1)
	{
		for (int m=0;m<i;m++) p+=S.n_i(state2,m);
	}
	return (flag*(1-2*(p%2)));
}

int matrixs::mFuncCX(long int state1, long int state2, int i)
{
	int flag=1, p=0;
	for (int m=0; m<S.N_b(); m++)
	{
		if( m == i )
		{
			if ((S.n_i(state2,i)==1) || (S.n_i(state1,i)==0) ) {flag=0;break;} else flag=1;
		}
		else
		{
			if (S.n_i(state1,m)==S.n_i(state2,m) ) flag=1; else {flag=0;break;}
		}
	}
	if (flag==1)
	{
		for (int m=0;m<i;m++) p+=S.n_i(state2,m);
	}
	return (flag*(1-2*(p%2)));
}

// other functions

void matrixs::print_to_screen()						//print to screen C and CX
{
	for (unsigned int i=0; i<UXCU.size(); i++)
		cout << UXCU[i].C << endl;
		
	for (unsigned int i=0; i<UXCXU.size(); i++)
		cout << UXCXU[i].C << endl;
		
}	

void matrixs::dump()							//writing matrixs C[M_sigma] and CX[M_sigma] in output file
{
	stringstream filename;
	filename << (*this).matrixC_path.path() << "//M_sig" << i << ".dat";
  	ofstream outHpart;
	outHpart.open(filename.str().c_str());
		
	for (unsigned int i=0; i<UXCU.size(); i++)
		outHpart << UXCU[i].n << "  " << UXCU[i].m << "  " << UXCU[i].C << endl;

	outHpart.close();
	
//	filename << std::flush;

	stringstream filename1;
	filename1 << (*this).matrixCX_path.path() << "//M_sig" << j << ".dat";
	outHpart.open(filename1.str().c_str());
		
	for (unsigned int i=0; i<UXCXU.size(); i++)
		outHpart << UXCXU[i].n << "  " << UXCXU[i].m << "  " << UXCXU[i].C << endl;
	
	outHpart.close();
}

string matrixs::path()
{
  return (*this).matrixCX_path.fullpath()+" ; "+(*this).matrixC_path.fullpath();
}

