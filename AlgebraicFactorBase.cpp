//
// Created by Илья Моржаретто on 14.05.17.
//

#include "AlgebraicFactorBase.h"
//#define EVALP_DEBUG

inline ZZ norm(const ZZ& a, const ZZ& b, const Polynomial& poly)
{
    std::vector<ZZ> powersb;
    powersb.resize(deg(poly.f)+1);
    ZZ powerb = ZZ(1);
    ZZ powera = ZZ(1);
    ZZ answer = ZZ(0);

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
#ifdef FACTOR_DEBUG
    std::cerr<<"Norm for ("<<a<<", "<<b<<") is "<<answer<<"\n";
#endif
    return answer;
}

bool AlgebraicFactorBase::factor(std::vector<long>& f, const ZZ& a, const ZZ& b)
{
    std::vector<long> zfactor;
    ZZ n = a+b*poly.m;
    ZZ nor = abs(norm(a, b, poly));
    if(fb.factor(zfactor, n))
    {
        f.resize(fb.r.size()+num);
        for(long i = 0; i<fb.r.size(); ++i)
        {
            f[i] = zfactor[i];
        }
        long j = 0;
        for(long i=0; i<this->a.size(); ++i)
        {
            if (this->a[i] != -1)
            {
                f[fb.r.size()+j] = 0;
                if(nor % fb.r[i] == 0)
                {
                    nor = nor / fb.r[i];
                    ++f[fb.r.size()+j];
                }
                while(nor % fb.r[i] == 0)
                {
                    nor = nor / fb.r[i];
                    ++f[fb.r.size()+j];
                }
                ++j;
                if(nor == 1)
                {
#ifdef FACTOR_DEBUG
                    std::cerr<<"Pair ("<<a<<", "<<b<<") factors completely\n";
#endif
                    return true;
                }
            }
        }
    }
    return false;
}

inline ZZ_p evalp(const Polynomial& poly, const ZZ& x)
{
#ifdef EVALP_DEBUG
    std::cerr<<"Polynomial is "<<poly.f<<"\n"<<"P is "<<ZZ_p::modulus()<<"\n"<<"X is "<<x<<"\n";
#endif
    ZZ_p answer = ZZ_p(1);
    ZZ_p xp = conv<ZZ_p>(x);
    for(long i=deg(poly.f)-1; i>=0; --i)
    {
        answer *= xp;
        answer += conv<ZZ_p>(poly.f[i]);
    }
#ifdef EVALP_DEBUG
    std::cerr<<"Answer is "<<rep(answer)<<"\n";
#endif
    return answer;
}

AlgebraicFactorBase::AlgebraicFactorBase(ZZ bound, const Polynomial& _f)
{
    poly = _f;
    fb.setBound(bound);
    a.resize(fb.r.size());
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
#ifdef GEN_DEBUG
                std::cerr<<"R for q = "<<q<<" is "<<r<<"\n";
#endif
                break;
            }
        }
    }
}

long AlgebraicFactorBase::getTotalSize() {
    return num+fb.r.size();
}
