#include <iostream>
#include<fstream>
#include<complex>

typedef std::complex<double> ComplexType;
typedef double RealType;

using namespace std;
int main(int argc, char *argv[])
{
  RealType acc=1e-10;
  ComplexType G;
  int Z1, Z2, W1, W1_, W2, W2_, N1, N1_, N2, N2_;
  std::ifstream g4_str(argv[1],ios::in | ios::binary);
  for(;;){
    g4_str.read(reinterpret_cast<char *>(&G),sizeof(std::complex<double>));
    if(abs(G)<0.1*acc*acc) break;
    g4_str.read(reinterpret_cast<char *>(&Z1),sizeof(int));
    g4_str.read(reinterpret_cast<char *>(&Z2),sizeof(int));

    g4_str.read(reinterpret_cast<char *>(&W1),sizeof(int));
    g4_str.read(reinterpret_cast<char *>(&W1_),sizeof(int));
    g4_str.read(reinterpret_cast<char *>(&W2),sizeof(int));
    g4_str.read(reinterpret_cast<char *>(&W2_),sizeof(int));

    g4_str.read(reinterpret_cast<char *>(&N1),sizeof(int));
    g4_str.read(reinterpret_cast<char *>(&N1_),sizeof(int));
    g4_str.read(reinterpret_cast<char *>(&N2),sizeof(int));
    g4_str.read(reinterpret_cast<char *>(&N2_),sizeof(int));
    std::cout<<"read:1 "<<Z1<<" "<<Z2<<" "<<W1<<" "<<W1_<<" "<<W2<<" "<<W2_<<" "<<N1<<" "<<N1_<<" "<<N2<<" "<<N2_<<" "<<G<<std::flush<<std::endl;
    }
}
