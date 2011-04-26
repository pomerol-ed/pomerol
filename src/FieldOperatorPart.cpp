#include <fstream>

#include "FieldOperatorPart.h"

using std::stringstream;

FieldOperatorPart::FieldOperatorPart(
        int i, StatesClassification &S, HamiltonianPart &h_from,  HamiltonianPart &h_to) : 
        ComputableObject(), i(i), S(S), h_from(h_from), h_to(h_to)
{};


ColMajorMatrixType& FieldOperatorPart::getColMajorValue()
{
    return elementsColMajor;
}

RowMajorMatrixType& FieldOperatorPart::getRowMajorValue()
{
    return elementsRowMajor;
}

void FieldOperatorPart::print_to_screen()  //print to screen C and CX
{
    QuantumNumbers to   = h_to.id();
    QuantumNumbers from = h_from.id();
    for (int P=0; P<elementsColMajor.outerSize(); ++P)
        for (ColMajorMatrixType::InnerIterator it(elementsColMajor,P); it; ++it)
        {
                QuantumState N = S.clstates(to)[it.row()];
                QuantumState M = S.clstates(from)[it.col()];
                std::cout << N <<" " << M << " : " << it.value() << std::endl;
        };
}

void FieldOperatorPart::dump() //writing FieldOperatorPart C[M_sigma] and CX[M_sigma] in output file
{
DEBUG("No dump method is implemented for field operators");
assert(0);
/*    std::stringstream filename;
    filename << (*this).OUT.fullpath() << "/" << "part" << i << "_" << h_from.id() << "->" << h_to.id() << ".dat";
    ofstream outCpart;
    outCpart.open(filename.str().c_str());
     for (int P=0; P<elementsColMajor.outerSize(); ++P)
        for (ColMajorMatrixType::InnerIterator it(elementsColMajor,P); it; ++it)
        {
                QuantumState N = it.row();//S.clstates(to)[it.row()];
                QuantumState M = it.col();//S.clstates(from)[it.col()];
                outCpart << S.clstates(h_to.id())[N] <<" " << S.clstates(h_from.id())[M] << "  " << it.value() << endl;
        };


 //   outCpart << std::setprecision(DUMP_FLOATING_POINT_NUMBERS) << elements.toDense() << endl;

    outCpart.close();
    cout << "The part of field operator " << h_from.id() << "->" << h_to.id() << " is dumped to " << filename.str() << "." << endl;
*/
}

void FieldOperatorPart::compute()
{
    QuantumNumbers to = h_to.id();
    QuantumNumbers from = h_from.id();

    const std::vector<QuantumState>& toStates = S.clstates(to);
    const std::vector<QuantumState>& fromStates = S.clstates(from);
    
    DynamicSparseMatrixType tempElements(toStates.size(),fromStates.size());
    
    for (std::vector<QuantumState>::const_iterator current_state = toStates.begin();
                                                   current_state < toStates.end(); current_state++)
    {     
        QuantumState L=*current_state;
        if (checkL(L))
        {
            int K = retK(L);
            if( (mFunc(L,K,i)!= 0) )
            {                                                 
                int l=S.getInnerState(L), k=S.getInnerState(K);                             // l,k in part of Hamilt                        
               
                for ( unsigned int n=0; n<toStates.size(); n++)
                {
                    if(h_to.reH(l,n)!=0)
                    {
                        for (unsigned int m=0; m<fromStates.size(); m++)
                        {
                            RealType C_nm = h_to.reH(l,n)*mFunc(L,K,i)*h_from.reH(k,m);
                            if (fabs(C_nm)>MatrixElementTolerance)
                            {       
                                tempElements.coeffRef(n,m) += C_nm;
                            }
                        }
                    }
                }
            }
        }       
    }
    tempElements.prune(MatrixElementTolerance);
    elementsRowMajor = tempElements;
    elementsColMajor = tempElements;
}

// Functions of specialized classes

AnnihilationOperatorPart::AnnihilationOperatorPart(int i, StatesClassification &S, 
                                                   HamiltonianPart &h_from, HamiltonianPart &h_to) :
FieldOperatorPart(i,S,h_from,h_to)
{}

CreationOperatorPart::CreationOperatorPart(int i, StatesClassification &S, 
                                                  HamiltonianPart &h_from, HamiltonianPart &h_to) :
FieldOperatorPart(i,S,h_from,h_to)
{}

QuantumState AnnihilationOperatorPart::retK(QuantumState L)							//return K for C
{	
	return( L + (1<<i) );
}

QuantumState CreationOperatorPart::retK(QuantumState L)							//return K for CX
{	

	return( L - (1<<i) );

}

int AnnihilationOperatorPart::mFunc(QuantumState state1, QuantumState state2, int i)
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

int CreationOperatorPart::mFunc(QuantumState state1, QuantumState state2, int i)
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


bool AnnihilationOperatorPart::checkL(QuantumState L)
{
	return (S.n_i(L,i)==0);
}

bool CreationOperatorPart::checkL(QuantumState L)
{
	return (S.n_i(L,i)==1);
}

BlockNumber FieldOperatorPart::getLeftIndex()
{
	return h_to.getId();
};

BlockNumber FieldOperatorPart::getRightIndex()
{
	return h_from.getId();
};

CreationOperatorPart& AnnihilationOperatorPart::transpose()
{
	CreationOperatorPart *CX = new CreationOperatorPart(i, S, h_to, h_from); // swapped h_to and h_from
	CX->elementsRowMajor = elementsRowMajor.transpose();
	CX->elementsColMajor = elementsColMajor.transpose();
	return *CX;
};

AnnihilationOperatorPart& CreationOperatorPart::transpose()
{
	AnnihilationOperatorPart *C = new AnnihilationOperatorPart(i, S, h_to, h_from); // swapped h_to and h_from
	C->elementsRowMajor = elementsRowMajor.transpose();
	C->elementsColMajor = elementsColMajor.transpose();
	return *C;
};
