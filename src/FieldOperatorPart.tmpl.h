//class FieldOperatorPart                               //rotates matrixes C and CX 
#include <sstream>
#include <fstream>
#include <iomanip>

template<class StorageType> 
FieldOperatorPart<StorageType>::FieldOperatorPart(
        int i, StatesClassification &S, HamiltonianPart &h_from,  HamiltonianPart &h_to, output_handle OUT) : 
        i(i), S(S), h_from(h_from), h_to(h_to), OUT(OUT)
{};

template<class StorageType>
StorageType& FieldOperatorPart<StorageType>::value()
{
    return elements;
}

// other functions

template<class StorageType> 
void FieldOperatorPart<StorageType>::print_to_screen()  //print to screen C and CX
{
    QuantumNumbers to   = h_to.id();
    QuantumNumbers from = h_from.id();
    for (int P=0; P<elements.outerSize(); ++P)
        for (typename StorageType::InnerIterator it(elements,P); it; ++it)
        {
                QuantumState N = S.clstates(to)[it.row()];
                QuantumState M = S.clstates(from)[it.col()];
                cout << N <<" " << M << " : " << it.value() << endl;
        };
}

template<class StorageType>
void FieldOperatorPart<StorageType>::dump() //writing FieldOperatorPart C[M_sigma] and CX[M_sigma] in output file
{
    std::stringstream filename;
    filename << (*this).OUT.fullpath() << "/" << "C" << i << "_" << h_from.id() << "->" << h_to.id() << ".dat";
    ofstream outCpart;
    outCpart.open(filename.str().c_str());
        
    outCpart << std::setprecision(DUMP_FLOATING_POINT_NUMBERS) << elements.toDense() << endl;

    outCpart.close();
    cout << "The part of field operator " << h_from.id() << "->" << h_to.id() << " is dumped to " << filename.str() << "." << endl;
}

template<class StorageType>
const string &FieldOperatorPart<StorageType>::path()
{
    static string str=(*this).OUT.fullpath(); return str;
}

template<class StorageType> 
RealType FieldOperatorPart<StorageType>::computeElement(const QuantumState row, const QuantumState col,
                                                        const QuantumNumbers &from, const QuantumNumbers &to)
{
    RealType value = 0.0;
    for (std::vector<QuantumState>::const_iterator current_state = S.clstates(to).begin(); 
                                                   current_state < S.clstates(to).end(); current_state++){
        QuantumState L = *current_state;
        if (checkL(L))
        {
            int K = retK(L);
            if(mFunc(L,K,i) != 0){
                InnerQuantumState l = S.getInnerState(L), k = S.getInnerState(K);             // l,k in part of Hamiltonian
                if(h_to.reH(l,row) != 0){
                    value += h_to.reH(l,row)*mFunc(L,K,i)*h_from.reH(k,col);
                }
            }
        }
    }
    
    return value;
}

template<> 
void FieldOperatorPart<RowMajorMatrixType>::compute()
{
    QuantumNumbers to   = h_to.id();
    QuantumNumbers from = h_from.id();
    elements.resize(S.clstates(to).size(),S.clstates(from).size());
    
    elements.startFill();
    for (QuantumState n=0; n<S.clstates(to).size(); n++)
    for (QuantumState m=0; m<S.clstates(from).size(); m++){
        RealType matrixElement = computeElement(n,m,from,to);
        if(fabs(matrixElement)>MATRIX_ELEMENT_TOLERANCE) elements.fill(n,m) = matrixElement;
    }
    elements.endFill();
}

template<> 
void FieldOperatorPart<ColMajorMatrixType>::compute()
{
    QuantumNumbers to   = h_to.id();
    QuantumNumbers from = h_from.id();
    elements.resize(S.clstates(to).size(),S.clstates(from).size());
    
    elements.startFill();
    for (QuantumState m=0; m<S.clstates(from).size(); m++)
    for (QuantumState n=0; n<S.clstates(to).size(); n++){
        RealType matrixElement = computeElement(n,m,from,to);
        if(fabs(matrixElement)>MATRIX_ELEMENT_TOLERANCE) elements.fill(n,m) = matrixElement;
    }
    elements.endFill();
}
