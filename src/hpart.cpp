#include "hpart.h"
#include "getStates.h"
#include <fstream>
#include <sstream>

#include <json/json.h>

using std::stringstream;

// class getHpart

RealType getHpart::reH(int m, int n)								//return  H(m,n)
{
	return H(m,n);
}
RealType getHpart::reV(int m)									//return V(m)
{
	return V(m);
}

QuantumNumbers getHpart::id()
{
	return hpart_id;
}

void getHpart::inigetHpart( RealType J_c, double U_c, double Us_c, double mu_c, double mus_c, 
RealType t_c, double ts_c, const string &ev_path_, const string &ef_path_ )		//initialization getHpart
{ 		
	F_0 = U_c - 4*J_c/3;
	F_2 = 25*J_c/3;
	U=U_c;
	J=J_c;

	Us=Us_c;
	mu=mu_c;
	mus=mus_c;
	t=t_c;
	ts=ts_c;

	N_state_m = S.clstates(hpart_id).size();

	ev_path = ev_path_;
	ef_path = ef_path_;

	putmatrix();
	putHamilt();
}

void getHpart::putmatrix()										// inicialization {Wn}
{
	W1 = new int * [S.N_b_m()/2];
	for (int i=0;i<S.N_b_m()/2;i++)
	{
		W1[i]= new int [S.N_b_m()/2];
		for (int j=0;j<S.N_b_m()/2;j++) W1[i][j]=0;
	}

	W2 = new int * [S.N_b_m()/2];
	for (int k=0;k<S.N_b_m()/2;k++)
	{
		W2[k]= new int [S.N_b_m()/2];
		for (int l=0;l<S.N_b_m()/2;l++) W2[k][l]=0;
	}

	W3 = new int * [S.N_b_m()/2];
	for (int m=0;m<S.N_b_m()/2;m++)
	{
		W3[m]= new int [S.N_b_m()/2];
		for (int n=0;n<S.N_b_m()/2;n++) W3[m][n]=0;
	}
	if (S.N_b_m()==2)
	{
		W1[0][0]=1; 

		W2[0][0]=0; 

		W3[0][0]=0; 
	}

	if (S.N_b_m()==6)
	{
		W1[0][0]=1; W1[0][1]=-2; W1[0][2]=1;
		W1[1][0]=-2; W1[1][1]=4; W1[1][2]=-2;
		W1[2][0]=1; W1[2][1]=-2; W1[2][2]=1;

		W2[0][0]=0; W2[0][1]=3; W2[0][2]=6;
		W2[1][0]=3; W2[1][1]=0; W2[1][2]=3;
		W2[2][0]=6; W2[2][1]=3; W2[2][2]=0;

		W3[0][0]=0; W3[0][1]=-3; W3[0][2]=0;
		W3[1][0]=-3; W3[1][1]=0; W3[1][2]=-3;
		W3[2][0]=0; W3[2][1]=-3; W3[2][2]=0;
	}
	if (S.N_b_m()==10)
	{
		W1[0][0]=0;  W1[0][1]=0;  W1[0][2]=0;  W1[0][3]=0;  W1[0][4]=0;
		W1[1][0]=0;  W1[1][1]=0;  W1[1][2]=0;  W1[1][3]=0;  W1[1][4]=0;
		W1[2][0]=0;  W1[2][1]=0;  W1[2][2]=0;  W1[2][3]=0;  W1[2][4]=0;
		W1[3][0]=0;  W1[3][1]=0;  W1[3][2]=0;  W1[3][3]=0;  W1[3][4]=0;
		W1[4][0]=0;  W1[4][1]=0;  W1[4][2]=0;  W1[4][3]=0;  W1[4][4]=0;

		W2[0][0]=0;  W2[0][1]=0;  W2[0][2]=0;  W2[0][3]=0;  W2[0][4]=0;
		W2[1][0]=0;  W2[1][1]=0;  W2[1][2]=0;  W2[1][3]=0;  W2[1][4]=0;
		W2[2][0]=0;  W2[2][1]=0;  W2[2][2]=0;  W2[2][3]=0;  W2[2][4]=0;
		W2[3][0]=0;  W2[3][1]=0;  W2[3][2]=0;  W2[3][3]=0;  W2[3][4]=0;
		W2[4][0]=0;  W2[4][1]=0;  W2[4][2]=0;  W2[4][3]=0;  W2[4][4]=0;
		
		W3[0][0]=0;  W3[0][1]=0;  W3[0][2]=0;  W3[0][3]=0;  W3[0][4]=0;
		W3[1][0]=0;  W3[1][1]=0;  W3[1][2]=0;  W3[1][3]=0;  W3[1][4]=0;
		W3[2][0]=0;  W3[2][1]=0;  W3[2][2]=0;  W3[2][3]=0;  W3[2][4]=0;
		W3[3][0]=0;  W3[3][1]=0;  W3[3][2]=0;  W3[3][3]=0;  W3[3][4]=0;
		W3[4][0]=0;  W3[4][1]=0;  W3[4][2]=0;  W3[4][3]=0;  W3[4][4]=0;
		
	}
		
}

void getHpart::putHamilt()
{	

	H.resize(N_state_m,N_state_m);				//creation of H[i][j]=0 
	H.setZero();
		
	for (int st=0; st<N_state_m; st++)
	{
		add_diag(st,F_0,F_2);				//interactions on multiorbital(diagonal elements)
		add_U(st,Us);					//interactions on s-orbital
		add_mu(st, U*(1.5+S.L())-(5*S.L())*J);		//half-filling
		add_mu(st, mu);					//chem. potential on multiorbital
		add_mus(st, mus);				//chem. potential on s-orbitals
		add_mus(st, Us/2.);				//half-filling on s-orbitals
	}

	
	(*this).add_hopping(Formula.getHoppingMatrix());
	
  	for ( int st1=0; st1<N_state_m; st1++)
  	{
    		for ( int st2=0; st2<st1; st2++)
		{
			add_nondiag(st1,st2,F_2);		//interactions on multiorbital(nondiagonal elements)
		}
	}

	H.part<Eigen::UpperTriangular>() =  H.marked<Eigen::LowerTriangular>().transpose();  // Symmetric matrix

}

// multiorbital functions:

void getHpart::add_diag(int st, RealType F_0, double F_2)				//diagonal elements of Hamilt
{
	int state = S.cst(hpart_id,st);					//real state	

	for (int i=1; i<S.N_b_m()/2; i++)					//rule of filling multiorbital
	{
		for (int j=0; j<i; j++)
		{
			if ( (i>=S.N_b()/2) && (j<S.N_b()/2) )			//electrons with spin "down" and "up"
    			{			
      				H(st,st)+=F_0*S.n_i(state,i)*S.n_i(state,j);	
     				H(st,st)+=(F_2/25)*W1[(i-S.N_b()/2)%(S.N_b_m()/2)][j%(S.N_b_m()/2)]*S.n_i(state,i)*S.n_i(state,j);
      			
				if ( (i+j) == (S.N_b()/2+S.N_b_m()/2-1) ) 	//electrons with m1=-m2
       				
					H(st,st)+=(F_2/25)*W3[j%(S.N_b_m()/2)][j%(S.N_b_m()/2)]*S.n_i(state,i)*S.n_i(state,j);
 			     	
				else
        				H(st,st)+=0;
			}  

			if (fabs(i-j)!=S.N_b()/2)				//electrons with m1!=m2
			{
				if ( (i>=S.N_b()/2) && (j<S.N_b()/2) )
				        H(st,st)+=0;
			      	else
			        	H(st,st)+=(F_0 - F_2/5)*S.n_i(state,i)*S.n_i(state,j);
			}
			else
			{
				if ( (i>=S.N_b()/2) && (j<S.N_b()/2) )
			        	H(st,st)+=(F_2/25)*W2[j%(S.N_b_m()/2)][j%(S.N_b_m()/2)]*S.n_i(state,i)*S.n_i(state,j);
			      	else
				     	H(st,st)+=0;
			}
		}
	}

	for (int i=S.N_b()/2; i<(S.N_b()/2+S.N_b_m()/2); i++)				//rule of filling multiorbital
	{
		for (int j=0; j<i; j++)
		{
			if ( (j<S.N_b_m()/2) || (j>=S.N_b()/2) )
			{
				if ( (i>=S.N_b()/2) && (j<S.N_b()/2) )			//electrons with spin "down" and "up"
    				{			
      					H(st,st)+=F_0*S.n_i(state,i)*S.n_i(state,j);	
     					H(st,st)+=(F_2/25)*W1[(i-S.N_b()/2)%(S.N_b_m()/2)][j%(S.N_b_m()/2)]*S.n_i(state,i)*S.n_i(state,j);
      			
					if ( (i+j) == (S.N_b()/2+S.N_b_m()/2-1) ) 	//electrons with m1=-m2
       				
						H(st,st)+=(F_2/25)*W3[j%(S.N_b_m()/2)][j%(S.N_b_m()/2)]*S.n_i(state,i)*S.n_i(state,j);
 			     	
					else
        					H(st,st)+=0;
				}  

				if (fabs(i-j)!=S.N_b()/2)				//electrons with m1!=m2
				{
					if ( (i>=S.N_b()/2) && (j<S.N_b()/2) )
					        H(st,st)+=0;
			      		else
			        		H(st,st)+=(F_0 - F_2/5)*S.n_i(state,i)*S.n_i(state,j);
				}
				else
				{
					if ( (i>=S.N_b()/2) && (j<S.N_b()/2) )
				        	H(st,st)+=(F_2/25)*W2[j%(S.N_b_m()/2)][j%(S.N_b_m()/2)]*S.n_i(state,i)*S.n_i(state,j);
			      		else
				     		H(st,st)+=0;
				}
			}
		}
	}
}

int getHpart::measurefunc(long int state1, long int state2, int i, int j, int k, int l)		//auxiliary function
{
	int flag=1, p=0;
	for (int m=0; m<S.N_b();m++)								//checking of "hopping" between state1 and state2
	{
		if ( (m==i) || (m==j) || (m==k) || (m==l) )
		{
			if ( (S.n_i(state2,i)==0) || (S.n_i(state1,i)==1) ) {flag=0;break;} else flag=1;
			if ( (S.n_i(state2,j)==0) || (S.n_i(state1,j)==1) ) {flag=0;break;} else flag=1;
			if ( (S.n_i(state2,k)==1) || (S.n_i(state1,k)==0) ) {flag=0;break;} else flag=1;
			if ( (S.n_i(state2,l)==1) || (S.n_i(state1,l)==0) ) {flag=0;break;} else flag=1;
		}
		else
		{
			if (S.n_i(state1,m)==S.n_i(state2,m) ) flag=1; else {flag=0;break;}
		}
	}
	if (flag==1)										//choose of true sign
	{
		for (int m=0; m<j; m++) p+=S.n_i(state2,m);
		for (int s=0; s<i; s++) p+=S.n_i(state2,s);
		for (int t=0; t<k; t++) p+=S.n_i(state1,t);
		for (int u=0; u<l; u++) p+=S.n_i(state1,u);
	}
	return (flag*(1-2*(p%2)));
}

int getHpart::inhopfuncW_2(long int state1, long int state2,int i, int j)		//type of overturn
{
	int result = measurefunc(state1, state2, i, j, j+S.N_b()/2, i-S.N_b()/2);
	return result;
}

int getHpart::inhopfuncW_3(long int state1, long int state2,int i, int j, int * p)	//type of overturn
{
	int result=0;
	for (int m=1; m<S.N_b_m()/2; m++)
	{
		if ( ( (i+m)<(S.N_b()/2+S.N_b_m()/2) ) && ((j-m)>=0) )			//state2<state1
		{
			result = measurefunc(state1,state2,i,j,i+m,j-m);
			if (result!=0) { *p=-m; break; }
		}
	}
	return result;
}

void getHpart::add_nondiag(int st1, int st2, RealType F_2)				//nondiagonal elements of Hamilt 
{
	int state1 = S.cst(hpart_id,st1);						//real state1
	int state2 = S.cst(hpart_id,st2);						//real state2

	for (int i=1; i<S.N_b_m()/2; i++)						//rule of filling multiorbital
	{
		for (int j=0; j<i; j++)
		{
			if ( (i>=S.N_b()/2) && (j<S.N_b()/2) )
      			{
        			H(st1,st2)+=(F_2/25)*W2[(i-S.N_b()/2)%(S.N_b_m()/2)][j%(S.N_b_m()/2)]*inhopfuncW_2(state1,state2,i,j);
       				
				if ( (i+j)==(S.N_b()/2+S.N_b_m()/2-1) )
				
					H(st1,st2)+=inhopfuncW_3(state1,state2,i,j,&hop)*(F_2/25)*W3[j%(S.N_b_m()/2)][(j+hop)%(S.N_b_m()/2)];
				else
       			  		H(st1,st2)+=0;
			}
		}
	}
	
	for (int i=S.N_b()/2; i<(S.N_b()/2+S.N_b_m()/2); i++)				//rule of filling multiorbital
	{
		for (int j=0; j<i; j++)
		{
			if ( (j<S.N_b_m()/2) || (j>=S.N_b()/2) )
			{
				if ( (i>=S.N_b()/2) && (j<S.N_b()/2) )
      		 		{
        				H(st1,st2)+=(F_2/25)*W2[(i-S.N_b()/2)%(S.N_b_m()/2)][j%(S.N_b_m()/2)]*inhopfuncW_2(state1,state2,i,j);
       					
					if ( (i+j)==(S.N_b()/2+S.N_b_m()/2-1) )
				
						H(st1,st2)+=inhopfuncW_3(state1,state2,i,j,&hop)*(F_2/25)*W3[j%(S.N_b_m()/2)][(j+hop)%(S.N_b_m()/2)];
					else
       			  			H(st1,st2)+=0;
				}
			}
		}
	}
}
 
// s-orbital functions

void getHpart::add_U(int st, RealType Us)				//interactions on s-orbital
{
	int state = S.cst(hpart_id,st);			//real state

	for (int i=S.N_b_m()/2; i<S.N_b()/2; i++)		//rule of filling s-orbital

		H(st,st)+=Us*S.n_i(state,i)*S.n_i(state,i+S.N_b()/2);
}

// chem. potentials

void getHpart::add_mu(int st, RealType mu)			//adds chem. potential on multiorbital
{
	int state = S.cst(hpart_id,st);			//real state
	
	for (int j=0; j<S.N_b_m()/2; j++)
		H(st,st)-=mu*S.n_i(state,j);
	
	for (int j=S.N_b()/2; j<(S.N_b()/2+S.N_b_m()/2); j++)
		H(st,st)-=mu*S.n_i(state,j);
	
}

void getHpart::add_mus(int st, RealType mus)			//adds chem. potential on s-orbitals
{
	int state = S.cst(hpart_id,st);			//real state
	
	for (int j=S.N_b_m()/2; j<S.N_b()/2; j++)
		H(st,st)-=mus*S.n_i(state,j);

	for (int j=(S.N_b()/2+S.N_b_m()/2); j<S.N_b(); j++)
		H(st,st)-=mus*S.n_i(state,j);

}

// hopping functions

int getHpart::checkhop(long int state1, long int state2, int i, int j)			//analog measurefunc
{
	int flag=1, p=0;
	for (int m=0; m<S.N_b(); m++)
	{
		if( (m==i) || (m==j) )
		{
			if ((S.n_i(state2,i)==0) || (S.n_i(state1,i)==1) ) {flag=0;break;} else flag=1;
			if ((S.n_i(state2,j)==1) || (S.n_i(state1,j)==0) ) {flag=0;break;} else flag=1;
		}
		else
		{
			if (S.n_i(state1,m)==S.n_i(state2,m) ) flag=1; else {flag=0;break;}
		}
	}
	if (flag==1)
	{
		for (int m=0;m<i;m++) p+=S.n_i(state2,m);
		for (int u=0;u<j;u++) p+=S.n_i(state1,u);
	}
	return (flag*(1-2*(p%2)));
}

int getHpart::hoppingfunc(long int state1, long int state2, int i)				//checks hopping between  state1 and state2
{												//set rules of hopping between s-orbitals and Lz=0
	int result=0;
	int orbital;

	if (S.N_b_m()!=0)
		orbital=(S.N_b_m()/2-1)/2;
	else
	        orbital=-1;

	if ( (i+1)<S.N_b()/2)
	{
		if( (i<orbital) || ((i>orbital)&&(i<S.N_b_m()/2)) )
			result=0;
		else
		{
			if (i==orbital)
				result=checkhop(state1,state2,i,S.N_b_m()/2);			//state2<state1
			else
				result=checkhop(state1,state2,i,i+1);			
		}
	}
	if( (i>=S.N_b()/2)&&((i+1)<S.N_b()) )
	{
		if( (i<(S.N_b()/2+orbital)) || ((i>(S.N_b()/2+orbital))&&(i<(S.N_b()/2+S.N_b_m()/2))) )
			result=0;
		else
		{
			if (i==S.N_b()/2+orbital)
				result=checkhop(state1,state2,i,S.N_b()/2+S.N_b_m()/2);		//state2<state1
			else
				result=checkhop(state1,state2,i,i+1);
		}
	}
	return result;
}

void getHpart::add_hopping(RealMatrixType& HoppingMatrix)
{
  for (int i=0;i<HoppingMatrix.rows();i++)
	  for (int j=0;j<HoppingMatrix.cols();j++)
		  if (i!=j) add_hopping(i,j,HoppingMatrix(i,j));
}

void getHpart::add_hopping(int i, int j, RealType t)
{

  for ( int st1=0; st1<N_state_m; st1++)
    {
	QuantumState state1 = S.cst(hpart_id,st1);				//real state1
	QuantumState difference = (i>j)?(1<<i)-(1<<j):(1<<j)-(1<<i);
	if (!((difference > S.N_st()-state1 && i>j) || (i<j && difference > state1 )))
	{
		QuantumState state2=(i>j)?state1+difference:state1-difference;
		if (S.getStateInfo(state1) == S.getStateInfo(state2))
	  	{
	    		int st2 = S.inner_state(state2);
	    		if (st2>=0) H(st1,st2) = t*checkhop(state1,state2,i,j);
	  	}
    	}
    }
}
void getHpart::add_hopping_everywhere(int st1, int st2, RealType t, double ts) 	//adds hopping elements in Hamilt
{
	int state1 = S.cst(hpart_id,st1);				//real state1
	int state2 = S.cst(hpart_id,st2);				//real state2

	for (int j=0; j<S.N_b(); j++)
	{
		
		if (S.N_b_m()!=0)
		{
			int orbital=(S.N_b_m()/2-1)/2;
		
			if ( (j == orbital) || (j == (S.N_b()/2 + orbital) ) )
				H(st1,st2)+=t*hoppingfunc(state1,state2,j);			//hopping between "Lz=0" and s-orbital
			else
				H(st1,st2)+=ts*hoppingfunc(state1,state2,j);			//hoppings between s-orbitals
		}
		else
			{
				if (hoppingfunc(state1,state2,j)!=0) cout << "Hopping : " << state1 << "->" << state2 << endl;
				H(st1,st2)+=ts*hoppingfunc(state1,state2,j);				//hoppings between s-orbitals 
			}
	}
}

//other functions

void getHpart::diagonalization()					//method of diagonalization classificated part of Hamiltonian
{
	if (N_state_m == 1)
	{	
		V = H;
	 	H(0,0) = 1;
	}
	if (N_state_m > 1)
	{
		Eigen::SelfAdjointEigenSolver<RealMatrixType> Solver(H, true);

		H = Solver.eigenvectors();				// eigenfunctions are ready
		V = Solver.eigenvalues();				// eigenvectors are ready
	}
}

void getHpart::print_to_screen()					//ptint part of Hamiltonian to screen
{
	cout << H << endl;
	cout << endl;
}

void getHpart::dump()							//writing Eigen Values in output file
{
	int Lz_max=0;
	for(int L=(S.N_b_m()/2-1)/2;L>0;L--)
		Lz_max+=2*L;

	if(N_state_m!=0)
	{
		stringstream filename;
		filename << (*this).ef_path << "//ef" << hpart_id << ".dat";
  		ofstream outHpart;
		outHpart.open(filename.str().c_str());
		outHpart << H << endl;
		outHpart << endl;
  		outHpart.close();
	}

	if(N_state_m!=0)
	{
		stringstream filename;
		filename << (*this).ev_path <<"//ev" << hpart_id << ".dat";
		ofstream outHpart;
  		outHpart.open(filename.str().c_str());
		outHpart << V << endl;
  		outHpart.close();
	}
}


