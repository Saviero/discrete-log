#include <NTL/RR.h>

#include "FactorBase.h"

//#define DEBUG
#ifdef DEBUG
#include<iostream>
#endif
// TODO: test new class

FactorBase::FactorBase(ZZ bound)
{
    setBound(bound);
}

inline ZZ rema(ZZ& q, ZZ a, ZZ b)
{
#ifdef DEBUG
    std::cerr<<"remainder params: "<<q<<" "<<a<<" "<<b<<"\n";
#endif
    q = a/b;
    ZZ r = a%b;
    while(r < 0)
    {
        r += b;
    }
    return r;
}

bool FactorBase::factor(std::vector<long>& f, const ZZ& _n)
{
    f.resize(r.size());
    for(int i = 0; i < r.size(); ++i)
    {
        f[i] = 0;
    }
    ZZ n(_n);
    ZZ q;
    for(long i=0; i<smallInd; ++i)
    {
        if(n % r[i] == 0)
        {
            n = n / r[i];
            ++f[i];
        }
        while(n % r[i]  == 0)
        {
            n = n / r[i];
            ++f[i];
        }

        if(n == 1)
        {
            return true;
        }
    }

    if(n <= r[r.size()-1])
    {
        for(long i=0; i<r.size(); ++i)
        {
            if (n == r[i])
            {
                ++f[i];
                return true;
            }
        }
    }

    for(long i=smallInd; i<r.size(); ++i)
    {
        if(n % r[i] == 0)
        {
            n = n / r[i];
            ++f[i];
        }
        while(n % r[i]  == 0)
        {
            n = n / r[i];
            ++f[i];
        }
        if(n <= r[r.size()-1])
        {
            if (n == 1)
            {
                return true;
            }
            for(long j=i+1; j<r.size(); ++j)
            {
                if (n == r[j])
                {
                    ++f[j];
                    return true;
                }
            }
        }
    }
    return false;
}

FactorBase::FactorBase() {
    count = 0;
    smallInd = -1;
    r.reserve(0);
}

void FactorBase::setBound(ZZ bound) {
    ZZ q = ZZ(2);
    ZZ elem;
    smallInd = 0;
    while(q<bound)
    {
        elem = q;
        r.push_back(elem);
        count++;
        q = NextPrime(q+1);
#ifdef DEBUG
            std::cerr<<q<<"\n";
#endif
    }
    while(r[smallInd] < TruncToZZ(sqrt(conv<RR>(r[r.size()-1]))))
    {
        ++smallInd;
    }
}
