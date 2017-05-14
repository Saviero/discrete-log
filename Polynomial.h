//
// Created by Илья Моржаретто on 28.04.17.
//

#ifndef DISCRETE_LOG_POLYNOMIALGENERATOR_H
#define DISCRETE_LOG_POLYNOMIALGENERATOR_H
#include<NTL/ZZ_p.h>
#include<NTL/ZZX.h>
#include<NTL/ZZ.h>
#include<NTL/RR.h>

using namespace NTL;

class Polynomial {
public:
    ZZ m;
    ZZ d;
    ZZX f;
    Polynomial();
};


#endif //DISCRETE_LOG_POLYNOMIALGENERATOR_H
