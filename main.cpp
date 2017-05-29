#include "FactorBase.h"
#include "Polynomial.h"
#include "AlgebraicFactorBase.h"
#include "NumberFieldSieve.h"

using namespace NTL;

int main()
{
    ZZ p = NextPrime(ZZ(1000000));
    std::cout<<p<<std::endl;
    ZZ_p::init(p);
    mat_ZZ* sieveres;
    std::vector<std::pair<ZZ, ZZ>> pairs;
    sieveres = sieve(pairs);


    /*long* schir;
    schir = schirokauer_map(ZZ(5), ZZ(7), ZZ(13), 3);
    for(int i=0; i<3; ++i)
    {
        std::cout<<schir[i]<<" ";
    }*/
    return 0;
}