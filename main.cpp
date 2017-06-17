#include "FactorBase.h"
#include "Polynomial.h"
#include "AlgebraicFactorBase.h"
#include "NumberFieldSieve.h"

using namespace NTL;

int main()
{
    ZZ p = NextPrime(ZZ(1019));
    std::cout<<p<<std::endl;
    ZZ_p::init(p);
    ZZ a = ZZ(-17);
    ZZ b = ZZ(1);
    ZZ mod = ZZ(2);
    std::cout<<(a % mod)<<" "<<(-b*b % mod)<<std::endl;
    std::cout<<log(ZZ(33), ZZ(625));
    return 0;
}