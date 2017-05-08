//
// Created by Илья Моржаретто on 28.04.17.
//

#include "PolynomialGenerator.h"

PolynomialGenerator::PolynomialGenerator()
{
    RR p = conv<RR>(ZZ_p::modulus());
    d = RoundToZZ(pow(log(p)/(log(log(p))), RR(1./3)));
    RR dinv = RR(1)/conv<RR>(d);
    m = TruncToZZ(pow(p, dinv));
}

ZZ_pX PolynomialGenerator::generate() {
    ZZ_pX poly;
    ZZ p = ZZ_p::modulus();
    for(int i=0; i<=d; ++i)
    {
        SetCoeff(poly, i, conv<ZZ_p>(p%m));
        p = p/m;
    }
    return poly;
}

ZZ PolynomialGenerator::getM()
{
    return m;
}