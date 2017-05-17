#include <NTL/RR.h>

#include "FactorBase.h"

#define DEBUG
#ifdef DEBUG
#include<iostream>
#endif
// TODO: test new class

FactorBase::FactorBase(long bound)
{
    setBound(bound);
}

inline ZZ rema(ZZ& q, ZZ a, ZZ b)
{
    q = a/b;
    ZZ r = a%b;
    while(r < 0)
    {
        r += b;
    }
    return r;
}

bool FactorBase::factor(long* f, const ZZ& _n, ZZ& rem)
{
    ZZ n(_n);
    if (f)
    {
        memset(f, 0, sizeof(long)*r.size());
    }
    rem = 0;
    ZZ q;
    for(long i=0; i<smallInd; ++i)
    {
        if(rema(q, n, r[i]) == 0)
        {
            n = q;
            if (f) ++f[i];
        }
        while(rema(q, n, r[i])  == 0)
        {
            n = q;
            if (f) ++f[i];
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
                if (f) ++f[i];
                return true;
            }
        }
    }

    for(long i=smallInd; i<r.size(); ++i)
    {
        if(rema(q, n, r[i]) == 0)
        {
            n = q;
            if (f) ++f[i];
        }
        while(rema(q, n, r[i])  == 0)
        {
            n = q;
            if (f) ++f[i];
        }
        if(n <= r[r.size()-1])
        {
            if (n == 1)
            {
                return true;
            }
            for(long j=i+1; j<r.size(); ++i)
            {
                if (n == r[j])
                {
                    if (f) ++f[j];
                    return true;
                }
            }
        }
    }
    rem = n;
    return false;
}

FactorBase::FactorBase() {
    count = 0;
    smallInd = -1;
    r.reserve(0);
}

void FactorBase::setBound(long bound) {
    long q = 2;
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
