
#ifndef DISCRETE_LOG_SMOOTHGENERATOR_H
#define DISCRETE_LOG_SMOOTHGENERATOR_H

#include"FactorBaseGenerator.h"
#include<utility>
#include<armadillo>

class SmoothGenerator {
public:
    void generate(FactorBaseGenerator base, ZZ bound);
private:
    std::vector<std::pair<ZZ, ZZ> > pairs;
};


#endif //DISCRETE_LOG_SMOOTHGENERATOR_H
