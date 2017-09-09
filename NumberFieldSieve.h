//
// Created by Илья Моржаретто on 17.05.17.
//

#ifndef DISCRETE_LOG_NUMBERFIELDSIEVE_H
#define DISCRETE_LOG_NUMBERFIELDSIEVE_H

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include "AlgebraicFactorBase.h"

//#define DEBUG_LOG
//#define DEBUG_LANCZOS
//#define DEBUG
//#define DEBUG_SCHIR

using namespace NTL;

class Logarithm {
private:
    Polynomial f;
    FactorBase rfb;
    AlgebraicFactorBase alfb;
    std::vector<std::pair<ZZ, ZZ>> pairs;
    mat_ZZ *sieveres;
    unsigned long dim;

    vec_ZZ_p gausselim(mat_ZZ_p& a);
    vec_ZZ schirokauer_map(const ZZ &a, const ZZ &b, const ZZ &l, Polynomial f);
    ZZ chineserem(std::vector<std::pair<ZZ, ZZ>> s);
    mat_ZZ* sieve(std::vector<std::pair<ZZ, ZZ>>& s,
                  const Polynomial& f,
                  const FactorBase rfb,
                  const AlgebraicFactorBase& alfb,
                  const ZZ& v);
public:
    Logarithm() : f(), rfb(), alfb(f)
    {
        ZZ p = ZZ_p::modulus();
        RR logp = NTL::log(conv<RR>(p));
        ZZ v = TruncToZZ(exp(pow(logp, RR(1./3))*pow(NTL::log(logp), RR(2./3))));
        rfb.setBound(v);
        alfb.generate(rfb);
        dim = rfb.getSize() + alfb.getSize()+conv<long>(f.d);
        sieveres = sieve(pairs, f, rfb, alfb, v);
    };

    ZZ log(ZZ t, ZZ g);

};

#endif //DISCRETE_LOG_NUMBERFIELDSIEVE_H
