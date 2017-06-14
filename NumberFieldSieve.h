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

using namespace NTL;

vec_ZZ schirokauer_map(const ZZ& a, const ZZ& b, const ZZ& l, Polynomial f);
ZZ chineserem(std::vector<std::pair<ZZ, ZZ>> s);
mat_ZZ* sieve(std::vector<std::pair<ZZ, ZZ>>& s);


#endif //DISCRETE_LOG_NUMBERFIELDSIEVE_H
