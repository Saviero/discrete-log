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
    std::vector<std::vector<ZZ>> a;
    long num;
public:
    FactorBase fb = FactorBase();
    AlgebraicFactorBase(ZZ bound, const Polynomial& f);
    bool factor(std::vector<long>& f, const ZZ& a, const ZZ& b)const;
    long getTotalSize()const;
};


#endif //DISCRETE_LOG_ALGEBRAICFACTORBASE_H
