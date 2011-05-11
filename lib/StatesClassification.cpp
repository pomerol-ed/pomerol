//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


#include "StatesClassification.h"

bool BlockNumber::operator<(const BlockNumber& rhs) const {return number<rhs.number;}
bool BlockNumber::operator==(const BlockNumber& rhs) const {return number==rhs.number;}

QuantumNumbers::QuantumNumbers(int LZ, int N_UP, int N_DOWN) : std::vector<short>(3)
{
    (*this)[0]=LZ;
    (*this)[1]=N_UP;
    (*this)[2]=N_DOWN;
}

QuantumNumbers::QuantumNumbers() : std::vector<short>()
{}

bool QuantumNumbers::operator<(const QuantumNumbers &rhs) const
{
    if((*this)[0]!= rhs[0]) return (*this)[0] < rhs[0];
    if((*this)[1] != rhs[1]) return (*this)[1] < rhs[1];
    return (*this)[2] < rhs[2];
}

bool QuantumNumbers::operator==(const QuantumNumbers &rhs) const 
{
    return (rhs[0] == (*this)[0] && rhs[1] == (*this)[1] && rhs[2] == (*this)[2] );
}

std::ostream& operator<<(std::ostream& output,const QuantumNumbers& out)
{
  const short& Lz = out[0];
  const short& N_up = out[1];
  const short& N_down = out[2];

  output << "(" << Lz << "," << N_up << "," << N_down << ")";
  return output;
}

//class StatesClassification

const std::vector<QuantumState>& StatesClassification::getQuantumStates( QuantumNumbers in )            //return st[in.Lz][in.N_up][in.N_down]
{
    return st[in[0]][in[1]][in[2]];
}


const QuantumState StatesClassification::getQuantumState( QuantumNumbers in, int m )            //return st[in.Lz][in.N_up][in.N_down][m]
{
    return st[in[0]][in[1]][in[2]][m];
}

const InnerQuantumState StatesClassification::getInnerState(QuantumState state)
{
  int ST=-1;                // "state" in part of Hamilt            
  for (unsigned int n=0; n<(*this).getQuantumStates((*this).getStateInfo(state)).size(); n++ )
    {    
        if ( (*this).getQuantumState((*this).getStateInfo(state),n)== state) ST=n;    //get ST
    }
 return ST;
}

bool StatesClassification::checkQuantumNumbers(QuantumNumbers in)
{
    return (QuantumToBlock.count(in)>0);
};


BlockNumber StatesClassification::getBlockNumber(QuantumNumbers in)
{
 /* if ((in.Lz*(N_bit/2 + 1)*(N_bit/2 + 1) + in.N_up*(N_bit/2 + 1) + in.N_down) <= ((int) BLOCKNUMBERLIMIT) 
    && 
      (in.Lz*(N_bit/2 + 1)*(N_bit/2 + 1) + in.N_up*(N_bit/2 + 1) + in.N_down) >= 0) 
  {  
     return num_bl[in.Lz*(N_bit/2 + 1)*(N_bit/2 + 1) + in.N_up*(N_bit/2 + 1) + in.N_down];    //number from formula
  }
  else return ERROR_BLOCK_NUMBER;
  */
    return (QuantumToBlock.count(in))?QuantumToBlock[in]:ERROR_BLOCK_NUMBER;
}

QuantumNumbers StatesClassification::getBlockInfo(BlockNumber in)
{
      return (BlockToQuantum.count(in))?BlockToQuantum[in]:ERROR_QUANTUM_NUMBERS;
}

void StatesClassification::getSiteInfo(int bit, int& lz, int& spin)
{
  spin = (bit>=N_bit/2)?-1:1; 
  bit%=(N_bit/2);
  lz = (bit>=N_bit_m/2)?0:bit-(N_bit_m/2-1)/2;
}
BlockNumber StatesClassification::NumberOfBlocks()
{
    return maximumBlockNumber_+1; 
}
const QuantumState StatesClassification::getNumberOfStates()
{
    return N_state;
}

void StatesClassification::compute()             //initalize StatesClassification class by sorting all quantum states in system
{
    N_bit = IndexInfo.getIndexSize(); 
    N_state = ( 1 << N_bit ); 
    N_bit_m=0;
    #warning : bad N_bit_m definition. Actually existence of N_bit_m is definetely bad.
      for (std::vector<SingleIndex*>::iterator it=IndexInfo.getSingleIndexList().begin(); it != IndexInfo.getSingleIndexList().end(); it++ ) 
    {
        if ((*it)->type == p) { N_bit_m = 6; break;};
        if ((*it)->type == d) { N_bit_m = 10; break; };
        if ((*it)->type == f) { N_bit_m = 14; break; };
    }
    
    std::vector<QuantumState> null;                        //creation null vectors
    int Lz_max=0;
    for(int L=(N_bit_m/2-1)/2;L>0;L--)
        Lz_max+=2*L;
    st = new std::vector<QuantumState> ** [2*Lz_max+1];
    for (int i=0; i<(2*Lz_max+1); i++)
    {
        st[i] = new std::vector<QuantumState> * [N_bit/2+1];
        for (int j=0; j<(N_bit/2 + 1); j++)
        {
            st[i][j] = new std::vector<QuantumState> [N_bit/2+1];
            for (int k=0; k<(N_bit/2 +1); k++)
            {
                st[i][j][k] = null;
            }
        }
    }
    for(QuantumState state=0;state<N_state;state++)            //get "classificated" vectors
    {
        int Lz_get=0, N_up_get=0, N_down_get=0;
        for (int k=0;k<N_bit;k++)
        {
            if (N_bit_m!=0)
            {
                if ( (k<N_bit_m/2) )
                    Lz_get+=n_i(state,k)*( k%(N_bit_m/2) - (N_bit_m/2-1)/2 );
                if ( (k>=N_bit/2) && (k<N_bit/2+N_bit_m/2) )
                    Lz_get+=n_i(state,k)*( (k-N_bit/2)%(N_bit_m/2) - (N_bit_m/2-1)/2 );                
            }    
            if (k<N_bit/2)
                N_up_get+=n_i(state,k);
            else
                N_down_get+=n_i(state,k);
        }
        Lz_get+= Lz_max;
        st[Lz_get][N_up_get][N_down_get].push_back(state);
        
    }


    BLOCKNUMBERLIMIT = 2*Lz_max*(N_bit/2 + 1)*(N_bit/2 + 1) + (N_bit/2)*(N_bit/2 + 1) + (N_bit/2);
    maximumBlockNumber_ = BLOCKNUMBERLIMIT;
    
    int iter = 0;

    for (int Lz=0; Lz<(2*Lz_max + 1); Lz++)
    {
        for (int N_up=0; N_up<(N_bit/2 + 1); N_up++)
        {
            for (int N_down=0; N_down<(N_bit/2 + 1); N_down++)
            {
                

                if(st[Lz][N_up][N_down].size()==0)
                {
                    iter++;
                    maximumBlockNumber_=maximumBlockNumber_-1;
                }
                else
                {
                    int num_f = Lz*(N_bit/2 + 1)*(N_bit/2 + 1) + N_up*(N_bit/2 + 1) + N_down;    //Assign a unique number for a given QuantumNumbers combination
                    QuantumNumbers tmp=QuantumNumbers(Lz,N_up,N_down);
                    BlockNumber tmpBlockNumber = num_f - iter;

                    BlockToQuantum[tmpBlockNumber] = tmp;
                    QuantumToBlock[tmp] = tmpBlockNumber;
                }
            }
        }
    }
}

QuantumNumbers StatesClassification::getStateInfo(QuantumState in)                        //returns Lz,N_up,N_down for number in
{
    
    if (in >= N_state) return ERROR_QUANTUM_NUMBERS;
    int Lz_max=0;                                //begining of calculating Lz,N_up,N_down
    
    for(int L=(N_bit_m/2-1)/2;L>0;L--)
        Lz_max+=2*L;

    int Lz_st=0, N_up_st=0, N_down_st=0;

    for (int p=0;p<N_bit;p++)
    {    
        if (N_bit_m!=0)
        {
            if ( (p<N_bit_m/2) )
                Lz_st+=(*this).n_i(in,p)*( p%(N_bit_m/2) - (N_bit_m/2-1)/2 );
            
            if ( (p>=N_bit/2) && (p<N_bit/2+N_bit_m/2) )
                Lz_st+=(*this).n_i(in,p)*( (p-N_bit/2)%(N_bit_m/2) - (N_bit_m/2-1)/2 );    
        }            
        
        if (p<N_bit/2)
            N_up_st+=(*this).n_i(in,p);
    
        else
            N_down_st+=(*this).n_i(in,p);
    }
    
    Lz_st+= Lz_max;                                //finish of calculating
    
    QuantumNumbers N(Lz_st,N_up_st,N_down_st);

    return N;

}

BlockNumber StatesClassification::getBlockNumber(QuantumState in)
{
    return (*this).getBlockNumber((*this).getStateInfo(in));
}


