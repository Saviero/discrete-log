#include "FactorBase.h"
#include "Polynomial.h"
#include "AlgebraicFactorBase.h"
#include "NumberFieldSieve.h"

using namespace NTL;

int main()
{
    ZZ p = NextPrime(ZZ(11026));
    std::cout<<p<<std::endl;
    ZZ_p::init(p);

    std::cout<<log(ZZ(2), ZZ(4096));
    return 0;
}