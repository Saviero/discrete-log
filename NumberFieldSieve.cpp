#include "NumberFieldSieve.h"
#include "AlgebraicFactorBase.h"
#define DEBUG

long* schirokauer_map(const ZZ& a, const ZZ& b, const ZZ& l, long d)
{
    ZZ mod = l*l;
    ZZ_pPush push;
    ZZ_p::init(mod);
    mat_ZZ_p matr;
    matr.SetDims(d, d);
    for(long i=0; i<d; ++i)
    {
        matr[i][i] = conv<ZZ_p>(l);
    }
    vec_ZZ_p v;
    v.SetLength(d);
    v[0] = conv<ZZ_p>(a);
    v[1] = conv<ZZ_p>(b);
    ZZ_p discr;
    vec_ZZ_p x;
    x.SetLength(d);
    solve(discr, x, matr, v);
    long* answer = new long;
    memset(answer, 0, sizeof(long)*d);
    for(int i=0; i<d; ++i)
    {
        answer[i] = conv<long>(rep(x[i]));
    }
    return answer;
}

mat_ZZ* sieve(std::vector<std::pair<ZZ, ZZ>>& s)
{
    ZZ p = ZZ_p::modulus();
    RR logp = log(conv<RR>(p));
    ZZ v = TruncToZZ(exp(pow(logp, RR(1./3))*pow(log(logp), RR(2./3))));
#ifdef DEBUG
    std::cerr<<"Bound is: "<<v<<"\n";
#endif
    Polynomial f;
#ifdef DEBUG
    std::cerr<<"Generating factor bases\n";
#endif
    AlgebraicFactorBase base(v, f);
#ifdef DEBUG
    std::cerr<<"FB done\n";
#endif
    mat_ZZ* res = new mat_ZZ();
    res->SetDims(base.getTotalSize()+conv<long>(f.d), base.getTotalSize()+conv<long>(f.d));
    long numOfRows = 1;

    for(ZZ b = ZZ(1); b < v; ++b)
    {
        for(ZZ a = -v; a < v; ++a)
        {
            if (GCD(a, b) == ZZ(1))
            {
                std::vector<long> expvec;
#ifdef DEBUG
                std::cerr<<"Trying to factor: ("<<a<<"; "<<b<<")\n";
#endif
                if(base.factor(expvec, a, b))
                {
#ifdef DEBUG
                    std::cerr<<"Factored\n";
#endif
                    std::pair<ZZ, ZZ> ab;
                    ab.first = a;
                    ab.second = b;
                    s.push_back(ab);
                    for(long i=0; i<base.getTotalSize(); ++i)
                    {
                        (*res)[numOfRows][i] = expvec[i]; // TODO: Add checking for numOfRows
                    }
                    ++numOfRows;
                }
                else
                {
#ifdef DEBUG
                    std::cerr<<"Not smooth\n";
#endif
                }
            }
        }
    }
    return res;
}

ZZ chineserem(std::vector<std::pair<ZZ, ZZ>> s)
{
    ZZ m = ZZ(1);
    for(int i=0; i<s.size(); ++i)
    {
        m *= s[i].second;
    }
    ZZ q, t, res, d;
    res = 0;
    for(int i=0; i<s.size(); ++i)
    {
        XGCD(d, q, t, s[i].second, m/s[i].second);
        res += (m/s[i].second)*s[i].first*t;
    }
    return res % m;
}