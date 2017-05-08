#ifndef DISCRETE_LOG_FACTORBASEGENERATOR_H
#define DISCRETE_LOG_FACTORBASEGENERATOR_H
#include<NTL/ZZ.h>
#include<NTL/ZZ_pX.h>
#include<vector>
#include<utility>
#include<fstream>

using namespace NTL;

class FactorBaseGenerator {
private:
    ZZ bound;
    long countA;
    long countR;
public:
    std::vector<std::pair<ZZ, ZZ> > R;
    std::vector<std::pair<ZZ, ZZ> > A;
    FactorBaseGenerator();
    void generate(ZZ bound, ZZ_pX f, ZZ_p m);

};


#endif //DISCRETE_LOG_FACTORBASEGENERATOR_H
