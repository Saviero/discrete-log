#include"FactorBaseGenerator.h"
#include "PolynomialGenerator.h"

using namespace NTL;

int main()
{
    ZZ p;
    std::cin>>p;
    ZZ_p::init(p);
    PolynomialGenerator fgen;
    ZZ_pX poly = fgen.generate();
    FactorBaseGenerator basegen;
    basegen.generate(ZZ(10000), poly, conv<ZZ_p>(fgen.getM()));
    
}