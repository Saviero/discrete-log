//
// Created by Илья Моржаретто on 28.04.17.
//

#include "Polynomial.h"

Polynomial::Polynomial()
{
    RR p = conv<RR>(ZZ_p::modulus());
    d = RoundToZZ(pow(log(p)/(log(log(p))), RR(1./3)));
    RR dinv = RR(1)/conv<RR>(d);
    m = TruncToZZ(pow(p, dinv));

    f.SetLength(conv<long>(d));

    ZZ mod = ZZ_p::modulus();
    for(int i=0; i<=d; ++i)
    {
        SetCoeff(f, i, mod%m);
        mod = mod/m;
    }
}