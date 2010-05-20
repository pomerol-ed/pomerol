//class FieldOperatorPart                               //rotates matrixes C and CX 
#include <sstream>
#include <fstream>
#include <iomanip>

template<int StorageOrder> 
FieldOperatorPart<StorageOrder>::FieldOperatorPart(
        int i, StatesClassification &S, HamiltonianPart &h_from,  HamiltonianPart &h_to, output_handle OUT) : 
        i(i), S(S), h_from(h_from), h_to(h_to), OUT(OUT)
{};

template<int StorageOrder>
Eigen::SparseMatrix<RealType,StorageOrder>& FieldOperatorPart<StorageOrder>::value()
{
    return elements;
}

// other functions

template<int StorageOrder> 
void FieldOperatorPart<StorageOrder>::print_to_screen()  //print to screen C and CX
{
    QuantumNumbers to   = h_to.id();
    QuantumNumbers from = h_from.id();
    for (int P=0; P<elements.outerSize(); ++P)
        for (typename Eigen::SparseMatrix<RealType,StorageOrder>::InnerIterator it(elements,P); it; ++it)
        {
                QuantumState N = S.clstates(to)[it.row()];
                QuantumState M = S.clstates(from)[it.col()];
                cout << N <<" " << M << " : " << it.value() << endl;
        };
}

template<int StorageOrder>
void FieldOperatorPart<StorageOrder>::dump() //writing FieldOperatorPart C[M_sigma] and CX[M_sigma] in output file
{
    std::stringstream filename;
    filename << (*this).OUT.fullpath() << "/" << "C" << i << "_" << h_from.id() << "->" << h_to.id() << ".dat";
    ofstream outCpart;
    outCpart.open(filename.str().c_str());
        
    outCpart << std::setprecision(DUMP_FLOATING_POINT_NUMBERS) << elements.toDense() << endl;

    outCpart.close();
    cout << "The part of field operator " << h_from.id() << "->" << h_to.id() << " is dumped to " << filename.str() << "." << endl;
}

template<int StorageOrder>
const string &FieldOperatorPart<StorageOrder>::path()
{
    static string str=(*this).OUT.fullpath(); return str;
}

template<int StorageOrder> 
void FieldOperatorPart<StorageOrder>::compute()
{
    QuantumNumbers to = h_to.id();
    QuantumNumbers from = h_from.id();

    const vector<QuantumState>& toStates = S.clstates(to);
    const vector<QuantumState>& fromStates = S.clstates(from);
    
    Eigen::DynamicSparseMatrix<RealType,StorageOrder> tempElements(toStates.size(),fromStates.size());
    
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
                            if (fabs(C_nm)>MATRIX_ELEMENT_TOLERANCE)
                            {       
                                tempElements.coeffRef(n,m) += C_nm;
                            }
                        }
                    }
                }
            }
        }       
    }
    
    tempElements.prune(MATRIX_ELEMENT_TOLERANCE);
    elements = tempElements;
}
