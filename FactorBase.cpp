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

bool FactorBase::factor(std::vector<long>& f, const ZZ& _n)const
{
#ifdef DEBUG
    std::cerr<<"Value of factored num is: "<<_n<<std::endl;
#endif
    if(_n == 0)
    {
        return false;
    }
    f.resize(r.size());
    for(int i = 0; i < r.size(); ++i)
    {
        f[i] = 0;
    }
    ZZ n(_n);
    if(n < 0)
    {
        f[0] = 1;
        n = -n;
    }
    else
    {
        f[0] = 0;
    }
    ZZ q;
    for(long i=1; i<smallInd; ++i)
    {
        while(n % r[i]  == 0)
        {
#ifdef DEBUG
            std::cerr<<"N="<<n<<" is divisible by r="<<r[i]<<std::endl;
#endif
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
    ZZ q = ZZ(-1);
    ZZ elem;
    smallInd = 0;
#ifdef DEBUG
    std::cerr<<"Rational factor base elements are:\n";
#endif
    while(q<bound)
    {
#ifdef DEBUG
        std::cerr<<q<<" ";
#endif
        elem = q;
        r.push_back(elem);
        count++;
        q = NextPrime(q+1);
    }
    while(r[smallInd] < TruncToZZ(sqrt(conv<RR>(r[r.size()-1]))))
    {
        ++smallInd;
    }
#ifdef DEBUG
    std::cerr<<"\n";
#endif
}
