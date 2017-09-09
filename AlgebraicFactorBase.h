//
// Created by Илья Моржаретто on 14.05.17.
//

#ifndef DISCRETE_LOG_ALGEBRAICFACTORBASE_H
#define DISCRETE_LOG_ALGEBRAICFACTORBASE_H
#include"FactorBase.h"
#include "Polynomial.h"
//#define FACTOR_DEBUG
//#define GEN_DEBUG

class AlgebraicFactorBase {
private:
    Polynomial poly;
    std::vector<std::vector<ZZ>> a;
    unsigned long num;
public:
    AlgebraicFactorBase(const Polynomial& f);
    void generate(const FactorBase fb);
    bool factor(std::vector<long>& f, const ZZ& a, const ZZ& b)const;
    unsigned long getSize()const;
};


#endif //DISCRETE_LOG_ALGEBRAICFACTORBASE_H
