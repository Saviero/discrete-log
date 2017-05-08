//
// Created by Илья Моржаретто on 28.04.17.
//

#ifndef DISCRETE_LOG_POLYNOMIALGENERATOR_H
#define DISCRETE_LOG_POLYNOMIALGENERATOR_H
#include<NTL/ZZ_p.h>
#include<NTL/ZZ_pX.h>
#include<NTL/ZZ.h>
#include<NTL/RR.h>



using namespace NTL;


class PolynomialGenerator {
private:
    ZZ m;
    ZZ d;
public:
    PolynomialGenerator();
    ZZ getM();
    ZZ_pX generate();
};


#endif //DISCRETE_LOG_POLYNOMIALGENERATOR_H
