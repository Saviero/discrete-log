#include "FactorBase.h"
#include "Polynomial.h"
#include "AlgebraicFactorBase.h"

using namespace NTL;

int main()
{
    ZZ p = NextPrime(ZZ(1000000000));
    std::cout<<p;
    ZZ_p::init(p);
    Polynomial f;
    long q = NextPrime(1000000);
    AlgebraicFactorBase base(q, f);
    return 0;
}