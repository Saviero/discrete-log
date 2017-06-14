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


    /*mat_ZZ* sieveres;
    std::vector<std::pair<ZZ, ZZ>> pairs;
    sieveres = sieve(pairs);
    std::cout<<(*sieveres)<<"\n";*/



    vec_ZZ schir;
    Polynomial f;
    ZZX poly;
    SetCoeff(poly, 0, ZZ(27));
    SetCoeff(poly, 1, ZZ(1));
    SetCoeff(poly, 2, ZZ(1));
    f.f = poly;
    f.d = ZZ(2);
    f.m = ZZ(31);
    schir = schirokauer_map(ZZ(9), ZZ(25), ZZ(509), f);
    for(int i=0; i<2; ++i)
    {
        std::cout<<schir[i]<<" ";
    }
    return 0;
}