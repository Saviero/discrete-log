#include "NumberFieldSieve.h"
//#define DEBUG

vec_ZZ schirokauer_map(const ZZ& a, const ZZ& b, const ZZ& l, Polynomial f) // TODO: make this work!
{
    // Moving to the field mod l
    ZZ_pPush push;
    ZZ_p::init(l);
    ZZ_pX poly = conv<ZZ_pX>(f.f);
#ifdef DEBUG
    std::cerr<<"Polynomial is: "<<poly<<std::endl;
#endif
    // Determining epsilon
    vec_pair_ZZ_pX_long factorf = SquareFreeDecomp(poly);
    ZZ epsilon = l-1;
    for(auto i=factorf.begin(); i<factorf.end(); ++i)
    {
        epsilon = epsilon*(l-i->b)/GCD(epsilon, (l-i->b));
    }
#ifdef DEBUG
    std::cerr<<"Epsilon is: "<<epsilon<<std::endl;
#endif;

    // Determining matrix M
    ZZ_p::init(l*l);
    mat_ZZ_p m;
    m.SetDims(conv<long>(f.d), conv<long>(f.d));
    ZZX abpoly, mulpoly;
    SetCoeff(abpoly, 0, a);
    SetCoeff(abpoly, 1, b);

#ifdef DEBUG
    std::cerr<<"ABpoly is: "<<abpoly<<std::endl;
#endif;

    for(unsigned long i=0; i<f.d; ++i)
    {
        ZZX mult;
        SetCoeff(mult, i, 1);

#ifdef DEBUG
        std::cerr<<"Mult is: "<<mult<<std::endl;
#endif;

        mulpoly = abpoly*mult % f.f;

#ifdef DEBUG
        std::cerr<<"Mulpoly is: "<<mulpoly<<std::endl;
#endif;
        for(unsigned long j = 0; j < f.d; ++j)
        {
            m[j][i] = conv<ZZ_p>(coeff(mulpoly, j));
        }
    }

#ifdef DEBUG
    std::cerr<<"Matrix M is:\n"<<m<<std::endl;
#endif

    // Calculating lambdas
    m = power(m, epsilon);

#ifdef DEBUG
    std::cerr<<"Matrix M power epsilon is:\n"<<m<<std::endl;
#endif

    vec_ZZ_p lambdas;
    vec_ZZ_p eye;
    eye.SetLength(conv<long>(f.d));
    eye[0] = ZZ_p(1);
    mul(lambdas, m, eye);

#ifdef DEBUG
    std::cerr<<"Lambdas is: "<<lambdas<<std::endl;
#endif

    lambdas = lambdas - eye;

#ifdef DEBUG
    std::cerr<<"Final lambdas is: "<<lambdas<<std::endl;
#endif
    vec_ZZ res = conv<vec_ZZ>(lambdas);
    for(auto i = res.begin(); i<res.end(); ++i)
    {
        (*i) /= l;
    }
    return res;
}

mat_ZZ* sieve(std::vector<std::pair<ZZ, ZZ>>& s, const Polynomial& f, const AlgebraicFactorBase& base, const ZZ& v)
{
#ifdef DEBUG
    std::cerr<<"Bound is: "<<v<<"\n";
#endif
    // Preparing output matrix
    mat_ZZ* res = new mat_ZZ();
    res->SetDims(base.getTotalSize()+conv<long>(f.d), base.getTotalSize()+conv<long>(f.d));
    long numOfRows = 1;

    for(ZZ b = ZZ(1); b < v; ++b)
    {
        if(numOfRows >= base.getTotalSize()+conv<long>(f.d))
        {
            break;
        }
        for(ZZ a = -v; a < v; ++a)
        {
            if(numOfRows >= base.getTotalSize()+conv<long>(f.d))
            {
                break;
            }
            if (GCD(a, b) == ZZ(1)) // we need only coprime a and b
            {
                std::vector<long> expvec;
#ifdef DEBUG
                std::cerr<<"Trying to factor: ("<<a<<"; "<<b<<")\n";
#endif
                if(base.factor(expvec, a, b)) // factoring a+bm and a+b*alpha
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
                        (*res)[numOfRows][i] = expvec[i];
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

ZZ log(ZZ t, ZZ g)
{
    Polynomial f;
    mat_ZZ* sieveres;

    ZZ p = ZZ_p::modulus();
    RR logp = log(conv<RR>(p));
    ZZ v = TruncToZZ(exp(pow(logp, RR(1./3))*pow(log(logp), RR(2./3))));

    AlgebraicFactorBase base(v, f);

    std::vector<std::pair<ZZ, ZZ>> pairs;
    sieveres = sieve(pairs, f, base, v);
    std::vector<long> factort;
    if(base.fb.factor(factort, t))
    {
        for(unsigned long i=0; i < base.fb.r.size(); ++i)
        {
            (*sieveres)[0][i] = factort.at(i);
        }
        for(unsigned long i=base.fb.r.size(); i < base.getTotalSize()+f.d; ++i)
        {
            (*sieveres)[0][i] = 0;
        }
    }

    ZZ q = p-1;
    ZZ bound = TruncToZZ(sqrt(conv<RR>(q)));
    std::vector<std::pair<ZZ, ZZ>> factor;
    std::pair<ZZ, ZZ> primediv;
    for(ZZ l = ZZ(2); l < bound; l = NextPrime(l+1))
    {
        primediv.first = l;
        primediv.second = 0;
        while(q % l == 0)
        {
            ++primediv.second;
            q /= l;
        }
        if(primediv.second > 0)
        {
            factor.push_back(primediv);
        }
    }

    std::vector<long> gfactor;
    vec_ZZ gfactorZZ;
    if(base.fb.factor(gfactor, g))
    {
        for(long i=0; i<gfactor.size(); ++i) {
            gfactorZZ.append(conv<ZZ>(gfactor[i]));
        }
    }
    for(long i=0; i<factor.size(); ++i)
    {
        q = factor[i].first;
        for(long j=0; j<pairs.size(); ++j)
        {
            vec_ZZ schirmap = schirokauer_map(pairs[j].first, pairs[j].second, q, f);
            for(long k=0; k < f.d; ++k)
            {
                (*sieveres)[j+1][k+base.getTotalSize()] = schirmap[k];
            }
        }
        ZZ det;
        vec_ZZ x;
        solve(det, x, *sieveres, gfactorZZ);
        // TODO: finish the method
    }
}