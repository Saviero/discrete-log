#include "FactorBase.h"
#include "Polynomial.h"
#include "AlgebraicFactorBase.h"
#include "NumberFieldSieve.h"

using namespace NTL;

int main()
{
    ZZ p = NextPrime(conv<ZZ>("1001003"));
    std::cout<<p<<std::endl;
    ZZ_p::init(p);

    std::cout<<log(ZZ(74375), ZZ(556480));
    return 0;
}