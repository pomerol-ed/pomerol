#include "HamiltonianPart.h"
#include "StatesClassification.h"
#include <fstream>
#include <sstream>

#include <json/json.h>

using std::stringstream;

// class HamiltonianPart

RealType HamiltonianPart::reH(int m, int n)								//return  H(m,n)
{
	return H(m,n);
}
RealType HamiltonianPart::reV(int m)									//return V(m)
{
	return V(m);
}

InnerQuantumState HamiltonianPart::size(void)
{
    return N_state_m;
}

QuantumNumbers HamiltonianPart::id()
{
	return hpart_id;
}

void HamiltonianPart::enter( double mu_c, double mus_c )		//initialization HamiltonianPart
{ 		
	mu=mu_c;
	mus=mus_c;

	N_state_m = S.clstates(hpart_id).size();

	H.resize(N_state_m,N_state_m);				//creation of H[i][j]=0 
	H.setZero();
		
	for (InnerQuantumState st=0; st<N_state_m; st++)
	{
//		add_mu(st, U*(1.5+S.L())-(5*S.L())*J);		//half-filling
//		add_mu(st, mu);					//chem. potential on multiorbital
//		add_mus(st, Us/2.);				//half-filling on s-orbitals
		
		// loop over terms
/*		std::list<Term*>::iterator it1;
		for ( it1=Formula.getTermsList().getTerms(2).begin() ; it1 != Formula.getTermsList().getTerms(2).end(); it1++ )
		{
			if (( *it1)->type == "n") add_nTerm(st,(nTerm*) *it1);
		};
*/		
		std::list<Term*>::iterator it2;
		for ( it2=Formula.getTermsList().getTerms(4).begin() ; it2 != Formula.getTermsList().getTerms(4).end(); it2++ )
		{
			if ( (*it2)->type == "nn") add_nnTerm(st,(nnTerm*) *it2);
			else 
			if ( (*it2)->type == "spinflip") add_spinflipTerm(st,(spinflipTerm*) *it2);
		};
	}

	
	(*this).add_hopping(Formula.getHoppingMatrix());
	
	H.part<Eigen::UpperTriangular>() =  H.marked<Eigen::LowerTriangular>().transpose();  // Symmetric matrix
//	H.part<Eigen::SelfAdjoint>() = H.eval();
}

void HamiltonianPart::add_nTerm(InnerQuantumState inner_state,nTerm *N)
{
	QuantumState state = S.cst(hpart_id,inner_state);					//real state	
      	H(inner_state,inner_state)+=N->Value*S.n_i(state,N->bit[0]);	
};
void HamiltonianPart::add_nnTerm(InnerQuantumState inner_state, nnTerm *T)
{
	QuantumState state = S.cst(hpart_id,inner_state);					//real state	
      	H(inner_state,inner_state)+=T->Value*S.n_i(state,T->bit[0])*S.n_i(state,T->bit[2]);	
};

void HamiltonianPart::add_spinflipTerm(InnerQuantumState inner_state, spinflipTerm *T)
{
	QuantumState in = S.cst(hpart_id,inner_state); //real state	
	if ( S.n_i(in,T->bit[0]) || S.n_i(in,T->bit[1]) || !S.n_i(in,T->bit[2]) || !S.n_i(in,T->bit[3])) { return; }
	QuantumState diff1 = (1<<T->bit[0]) + (1<<T->bit[1]);
	QuantumState diff2 = (1<<T->bit[2]) + (1<<T->bit[3]);
	if ((diff2 > in + diff1) || ((diff2 <= in) && (S.N_st() - diff1 <= in - diff2))) return;
	QuantumState out = in + diff1 - diff2;	
	if (out<in ) return;
//	#warning The fact that hamiltonian is hermitian isn't taken into account. There may be slight perfomance issues.
	QuantumNumbers out_info = S.getStateInfo(out);
	if (out_info==(QuantumNumbers) hpart_id) 
	{
		InnerQuantumState out_inner_state = S.getInnerState(out);
		H(out_inner_state,inner_state)+=T->Value*measurefunc(in,out,T->bit[0],T->bit[1],T->bit[2],T->bit[3]);
//		DEBUG(in << "->" << out << " " << inner_state << "->" << out_inner_state << " | Term : " << (Term&) *T << " | " << T->Value*measurefunc(in,out,T->bit[0],T->bit[1],T->bit[2],T->bit[3]) << " = " << H(out_inner_state,inner_state));
	}
	

};

int HamiltonianPart::measurefunc(QuantumState state1, QuantumState state2, int i, int j, int k, int l)		//auxiliary function
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

// s-orbital functions

//void HamiltonianPart::add_U(int st, unsigned short bit_up, unsigned short bit_down, RealType Us)				//interactions on s-orbital
// chem. potentials

void HamiltonianPart::add_mu(int st, RealType mu)			//adds chem. potential on multiorbital
{
	int state = S.cst(hpart_id,st);			//real state
	
	for (int j=0; j<S.N_b_m()/2; j++)
		H(st,st)-=mu*S.n_i(state,j);
	
	for (int j=S.N_b()/2; j<(S.N_b()/2+S.N_b_m()/2); j++)
		H(st,st)-=mu*S.n_i(state,j);
	
}

void HamiltonianPart::add_mus(int st, RealType mus)			//adds chem. potential on s-orbitals
{
	int state = S.cst(hpart_id,st);			//real state
	
	for (int j=S.N_b_m()/2; j<S.N_b()/2; j++)
		H(st,st)-=mus*S.n_i(state,j);

	for (int j=(S.N_b()/2+S.N_b_m()/2); j<S.N_b(); j++)
		H(st,st)-=mus*S.n_i(state,j);

}

// hopping functions

int HamiltonianPart::checkhop(long int state1, long int state2, int i, int j)			//analog measurefunc
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


void HamiltonianPart::add_hopping(RealMatrixType& HoppingMatrix)
{
  for (int i=0;i<HoppingMatrix.rows();i++)
	  for (int j=0;j<HoppingMatrix.cols();j++)
		  if (i!=j) add_hopping(i,j,HoppingMatrix(i,j));
}

void HamiltonianPart::add_hopping(int i, int j, RealType t)
{

  for ( InnerQuantumState st1=0; st1<N_state_m; st1++)
    {
	QuantumState state1 = S.cst(hpart_id,st1);				//real state1
	QuantumState difference = (i>j)?(1<<i)-(1<<j):(1<<j)-(1<<i);
	if (!((difference > S.N_st()-state1 && i>j) || (i<j && difference > state1 )))
	{
		QuantumState state2=(i>j)?state1+difference:state1-difference;
		if (S.getStateInfo(state1) == S.getStateInfo(state2))
	  	{
	    		int st2 = S.getInnerState(state2);
	    		if (st2>=0) H(st1,st2)+= t*checkhop(state1,state2,i,j);
	  	}
    	}
    }
}

//other functions

void HamiltonianPart::diagonalization()					//method of diagonalization classificated part of Hamiltonian
{
	if (N_state_m == 1)
	{	
		V = H;
	 	H(0,0) = 1;
	}
	if (N_state_m > 1)
	{
/*		cout << hpart_id << endl;
		cout << "___________" << endl;
		cout << H << " " << endl;
		cout << "___________" << endl;
	*/
		Eigen::SelfAdjointEigenSolver<RealMatrixType> Solver(H, true);
		H = Solver.eigenvectors();
		V = Solver.eigenvalues();				// eigenvectors are ready
	}
}

void HamiltonianPart::print_to_screen()					//ptint part of Hamiltonian to screen
{
	cout << H << endl;
	cout << endl;
}

void HamiltonianPart::dump()							//writing Eigen Values in output file
{
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

