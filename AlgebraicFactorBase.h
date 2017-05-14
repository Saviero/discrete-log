//
// Created by Илья Моржаретто on 14.05.17.
//

#ifndef DISCRETE_LOG_ALGEBRAICFACTORBASE_H
#define DISCRETE_LOG_ALGEBRAICFACTORBASE_H
#include"FactorBase.h"


class AlgebraicFactorBase {
private:
    FactorBase fb;
    std::vector<ZZ> a;
public:
    AlgebraicFactorBase(ZZ bound);


};


#endif //DISCRETE_LOG_ALGEBRAICFACTORBASE_H
