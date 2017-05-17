//
// Created by Илья Моржаретто on 14.05.17.
//

#ifndef DISCRETE_LOG_ALGEBRAICFACTORBASE_H
#define DISCRETE_LOG_ALGEBRAICFACTORBASE_H
#include"FactorBase.h"
#include "Polynomial.h"


class AlgebraicFactorBase {
private:
    Polynomial poly;
    FactorBase fb = FactorBase();
    std::vector<ZZ> a;
    long num;
public:
    AlgebraicFactorBase(long bound, const Polynomial& f);
    bool factor(long* f, const ZZ& a, ZZ& b, ZZ& rem);
};


#endif //DISCRETE_LOG_ALGEBRAICFACTORBASE_H
