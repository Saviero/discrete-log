#include "AlgebraicFactorBase.h"

inline ZZ norm(const ZZ& a, const ZZ& b, const Polynomial& poly)
{
    /*
     * Calculating norm N(a) for ideals <a+b*alpha>.
     */
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

bool AlgebraicFactorBase::factor(std::vector<long>& f, const ZZ& a, const ZZ& b)const {
    /*
     * Factorization of an ideal <a+b*alpha> over algebraic factor base.
     */
    ZZ nor = abs(norm(a, b, poly));
#ifdef FACTOR_DEBUG
    std::cerr<<"Factoring ("<<a<<", "<<b<<")\n";
#endif
    f.resize(num);
    long j = 0;
    for(long i=0; i<this->a.size(); ++i)
    {
        for(unsigned long k=1; k<this->a[i].size(); ++k)
        {
            f[j] = 0;
#ifdef FACTOR_DEBUG
            std::cerr<<"For ideal q="<<this->a[i][0]<<" r="<<this->a[i][k]<<": a mod q = "<<(a % this->a[i][0])<<"; -b*r mod q = "<<((-b*this->a[i][k]) % this->a[i][0])<<std::endl;
#endif
            if((a % this->a[i][0]) == ((-b*this->a[i][k]) % this->a[i][0])) {
                while (nor % this->a[i][0] == 0) {
                    nor = nor / this->a[i][0];
                    ++f[j];
                }
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

AlgebraicFactorBase::AlgebraicFactorBase(const Polynomial &f)
{
    poly = f;
    a.resize(0);
    num = 0;
}

unsigned long AlgebraicFactorBase::getSize()const {
    return num;
}

void AlgebraicFactorBase::generate(const FactorBase fb)
{
    a.resize(fb.r.size());
    num = 0;
    for(long i=1; i<fb.r.size(); ++i)
    {
        ZZ q = fb.r[i];
        ZZ_pPush push;
        ZZ_p::init(q);
        for(ZZ r = ZZ(0); r < q; ++r)
        {
            ZZ_p x = evalp(poly, r);
            if(rep(x) == 0)
            {
                if(a[i].size() == 0)
                {
                    a[i].push_back(q);
                }
                a[i].push_back(r);
                ++num;
#ifdef GEN_DEBUG
                std::cerr<<"R for q = "<<q<<" is "<<r<<"\n";
#endif
            }
        }
    }
}