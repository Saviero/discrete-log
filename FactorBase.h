#ifndef DISCRETE_LOG_FACTORBASEGENERATOR_H
#define DISCRETE_LOG_FACTORBASEGENERATOR_H
#include<NTL/ZZ.h>
#include<NTL/ZZ_pX.h>
#include<vector>
#include<utility>
#include<fstream>
//#define DEBUG
using namespace NTL;

class FactorBase {
private:
    long count;
    long smallInd;
public:
    std::vector<ZZ> r;
    FactorBase(ZZ bound);
    FactorBase();
    void setBound(ZZ bound);
    bool factor(std::vector<long>& f, const ZZ& n)const;
    unsigned long getSize()const;
};

inline ZZ rema(ZZ& q, ZZ a, ZZ b);
#endif //DISCRETE_LOG_FACTORBASEGENERATOR_H
