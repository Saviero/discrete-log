//
// Created by Илья Моржаретто on 14.05.17.
//

#include "AlgebraicFactorBase.h"
//#define DEBUG

inline ZZ norm(const ZZ& a, const ZZ& b, const Polynomial& poly)
{
    ZZ* powersb = new ZZ;
    ZZ powerb = ZZ(1);
    ZZ powera = ZZ(1);
    ZZ answer = ZZ(0);

    memset(powersb, 0, (deg(poly.f)+1)*sizeof(ZZ));
    for(long i=0; i<=deg(poly.f); ++i)
    {
        powersb[i] = powerb;
        powerb *= (-b);
    }

    for(long i=0; i<=deg(poly.f); ++i)
    {
        answer += powera*poly.f[i]*powersb[deg(poly.f) - i];
        powera *= a;
    }

    return answer;
}

bool AlgebraicFactorBase::factor(long *f, const ZZ& a, ZZ& b, ZZ &rem)
{
    long* zfactor = new long;
    ZZ n = a+b*poly.m;
    ZZ nor = norm(a, b, poly);
    ZZ zrem;
    if(fb.factor(zfactor, n, zrem))
    {
        if (f)
        {
            memset(f, 0, sizeof(long)*(fb.r.size()+num));
        }
        for(long i = 0; i<fb.r.size(); ++i)
        {
            f[i] = zfactor[i];
        }
        ZZ q;
        for(long i=0; i<this->a.size(); ++i)
        {
            if (this->a[i] != -1)
            {
                if(rema(q, n, fb.r[i]) == 0)
                {
                    nor = q;
                    if (f) ++f[fb.r.size()+i];
                }
                while(rema(q, n, fb.r[i])  == 0)
                {
                    nor = q;
                    if (f) ++f[fb.r.size()+i];
                }

                if(nor == 1)
                {
                    return true;
                }
            }
        }
    }
    else
    {
        return false;
    }

}

inline ZZ_p evalp(const Polynomial& poly, const ZZ& x)
{
    ZZ_p answer = ZZ_p(0);
    ZZ_p xp = conv<ZZ_p>(x);
    for(long i=deg(poly.f); i>=0; ++i)
    {
        answer *= xp;
        answer += conv<ZZ_p>(poly.f[i]);
    }
    return answer;
}

AlgebraicFactorBase::AlgebraicFactorBase(long bound, const Polynomial& _f)
{
    poly = _f;
    fb.setBound(bound);
    a.reserve(fb.r.size());
    num = 0;
    for(long i=0; i<a.size(); ++i)
    {
        a[i] = -1;
    }
    for(long i=0; i<fb.r.size(); ++i)
    {
        ZZ q = fb.r[i];
        ZZ_pPush push;
        ZZ_p::init(q);
        for(ZZ r = ZZ(0); r < q; ++r)
        {
            ZZ_p x = evalp(poly, r);
            if(rep(x) == 0)
            {
                a[i] = r;
                ++num;
#ifdef DEBUG
                std::cerr<<r<<"\n";
#endif
                break;
            }
        }
    }
}