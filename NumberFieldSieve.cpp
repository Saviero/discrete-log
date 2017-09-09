#include "NumberFieldSieve.h"

vec_ZZ_p Logarithm::gausselim(mat_ZZ_p& a)
{
    vec_ZZ_p ans;
    ans.SetLength(a.NumRows());
    long m = a.NumCols()-1;
    long n = a.NumRows();
    std::vector<long> where(m, -1);
    for (long col = 0, row = 0; col < m && row < n; ++col) {
        long sel = row;
        for (long i = row; i < n; ++i) {
            if (rep(a[i][col]) > rep(a[sel][col])) {
                swap (a[i], a[row]);
                break;
            }
        }

        if (IsZero(a[row][col])) {
            continue;
        }
        where[col] = row;

        for (long i = 0; i < n; ++i) {
            if (i != row) {
                ZZ_p c = a[i][col] / a[row][col];
                for (long j=col; j<=m; ++j)
                    a[i][j] -= a[row][j] * c;
            }
        }
        ++row;
    }
    for (int i=0; i<m; ++i) {
        if (where[i] != -1) {
            ans[i] = a[where[i]][m] / a[where[i]][i];
        }
    }
    return ans;
}

vec_ZZ Logarithm::schirokauer_map(const ZZ& a, const ZZ& b, const ZZ& l, Polynomial f)
{
    // Moving to the field mod l
    ZZ_pPush push;
    ZZ_p::init(l);
    ZZ_pX poly = conv<ZZ_pX>(f.f);
#ifdef DEBUG_SCHIR
    std::cerr<<"Polynomial is: "<<poly<<", L is "<<l<<std::endl;
#endif
    // Determining epsilon
    vec_pair_ZZ_pX_long factorf = SquareFreeDecomp(poly);
#ifdef DEBUG_SCHIR
    std::cerr<<"Factorisation of F mod L is "<<factorf<<std::endl;
#endif
    ZZ epsilon = l-1;
    for(auto i=factorf.begin(); i<factorf.end(); ++i)
    {
        epsilon = epsilon*(l-i->b)/GCD(epsilon, (l-i->b));
    }
#ifdef DEBUG_SCHIR
    std::cerr<<"Epsilon is: "<<epsilon<<std::endl;
#endif;

    // Determining matrix M
    ZZ_p::init(l*l);
    mat_ZZ_p m;
    m.SetDims(conv<long>(f.d), conv<long>(f.d));
    ZZX abpoly, mulpoly;
    SetCoeff(abpoly, 0, a);
    SetCoeff(abpoly, 1, b);

#ifdef DEBUG_SCHIR
    std::cerr<<"ABpoly is: "<<abpoly<<std::endl;
#endif;

    for(unsigned long i=0; i<f.d; ++i)
    {
        ZZX mult;
        SetCoeff(mult, i, 1);

#ifdef DEBUG_SCHIR
        std::cerr<<"Mult is: "<<mult<<std::endl;
#endif;

        mulpoly = abpoly*mult % f.f;

#ifdef DEBUG_SCHIR
        std::cerr<<"Mulpoly is: "<<mulpoly<<std::endl;
#endif;
        for(unsigned long j = 0; j < f.d; ++j)
        {
            m[j][i] = conv<ZZ_p>(coeff(mulpoly, j));
        }
    }

#ifdef DEBUG_SCHIR
    std::cerr<<"Matrix M is:\n"<<m<<std::endl;
#endif

    // Calculating lambdas
    m = power(m, epsilon);

#ifdef DEBUG_SCHIR
    std::cerr<<"Matrix M power epsilon is:\n"<<m<<std::endl;
#endif

    vec_ZZ_p lambdas;
    vec_ZZ_p eye;
    eye.SetLength(conv<long>(f.d));
    eye[0] = ZZ_p(1);
    mul(lambdas, m, eye);

#ifdef DEBUG_SCHIR
    std::cerr<<"Lambdas is: "<<lambdas<<std::endl;
#endif

    lambdas = lambdas - eye;

#ifdef DEBUG_SCHIR
    std::cerr<<"Final lambdas is: "<<lambdas<<std::endl;
#endif
    vec_ZZ res = conv<vec_ZZ>(lambdas);
    for(auto i = res.begin(); i<res.end(); ++i)
    {
        (*i) /= l;
    }
    return res;
}

mat_ZZ* Logarithm::sieve(std::vector<std::pair<ZZ, ZZ>>& s,
                         const Polynomial& f,
                         const FactorBase rfb,
                         const AlgebraicFactorBase& alfb,
                         const ZZ& v)
{
#ifdef DEBUG
    std::cerr<<"Bound is: "<<v<<"\n";
#endif
    // Preparing output matrix
    mat_ZZ* res = new mat_ZZ();
    res->SetDims(rfb.getSize() + alfb.getSize()+conv<long>(f.d), rfb.getSize() + alfb.getSize()+conv<long>(f.d));
    long numOfRows = 1;
    ZZ newv = v;
    ZZ oldv = ZZ(1);
    while(numOfRows < dim) {
        for (ZZ a = oldv; a < newv; ++a) {
            if (numOfRows >= dim) {
                break;
            }
            for (ZZ b = -newv; b < newv; ++b) {
                if (b == 0) {
                    continue;
                }
                if (numOfRows >= dim) {
                    break;
                }
                if (GCD(a, b) == ZZ(1)) // we need only coprime a and b
                {
                    std::vector<long> expvec;
                    std::vector<long> rexpvec;
                    std::vector<long> alexpvec;

#ifdef DEBUG
                    std::cerr<<"Trying to factor: ("<<a<<"; "<<b<<")\n";
#endif
                    if (rfb.factor(rexpvec, a+b*f.m)) // factoring a+bm
                    {
#ifdef DEBUG
                        std::cerr<<"Rexpvec for a="<<a<<", b="<<b<<" is: ";
                        for(unsigned long i=0; i<rexpvec.size(); ++i)
                        {
                            std::cerr<<rexpvec[i]<<" ";
                        }
                        std::cerr<<"\n";
#endif
                        for(unsigned long i=0; i<rexpvec.size(); ++i)
                        {
                            expvec.push_back(rexpvec[i]);
                        }
                        if (alfb.factor(alexpvec, a, b)) // factoring <a+b*alpha>
                        {
#ifdef DEBUG
                            std::cerr<<"Alexpvec for a="<<a<<", b="<<b<<" is: ";
                            for(unsigned long i=0; i<alexpvec.size(); ++i)
                            {
                                std::cerr<<alexpvec[i]<<" ";
                            }
                            std::cerr<<"\n";
#endif
#ifdef DEBUG
                            std::cerr<<"Factored\n";
#endif
                            for(unsigned long i=0; i<alexpvec.size(); ++i)
                            {
                                expvec.push_back(alexpvec[i]);
                            }
#ifdef DEBUG
                            std::cerr<<"Expvec for a="<<a<<", b="<<b<<" is: ";
                            for(unsigned long i=0; i<expvec.size(); ++i)
                            {
                                std::cerr<<expvec[i]<<" ";
                            }
                            std::cerr<<"\n";
#endif
                            std::pair<ZZ, ZZ> ab;
                            ab.first = a;
                            ab.second = b;
                            s.push_back(ab);
                            for (long i = 0; i < dim - conv<long>(f.d); ++i) {
                                (*res)[i][numOfRows] = expvec[i];
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

        }
        ZZ temp;
        temp = newv;
        newv = 2*newv;
        oldv = temp;
    }
    return res;
}

ZZ Logarithm::chineserem(std::vector<std::pair<ZZ, ZZ>> s)
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


ZZ Logarithm::log(ZZ t, ZZ g)
{
    std::vector<long> factort;
    if(rfb.factor(factort, t))
    {
#ifdef DEBUG_LOG
    std::cerr<<"tfactor's results are: ";
    for(int i=0; i<factort.size(); ++i)
    {
        std::cerr<<factort[i]<<" ";
    }
#endif
        for(unsigned long i=0; i < rfb.getSize(); ++i)
        {
            (*sieveres)[i][0] = factort.at(i);
        }
        for(unsigned long i=rfb.getSize(); i < rfb.getSize()  +alfb.getSize() + f.d; ++i)
        {
            (*sieveres)[i][0] = 0;
        }
    }
    else
    {
        std::cerr<<"T is not smooth!\n";
        return ZZ(-1);
    }

    ZZ q = ZZ_p::modulus()-1;
    ZZ bound = TruncToZZ(sqrt(conv<RR>(q)));
    std::vector<std::pair<ZZ, ZZ>> factor;
    std::pair<ZZ, ZZ> primediv;
    for(ZZ l = ZZ(2); l <= bound; l = NextPrime(l+1))
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
    if(q > ZZ(1))
    {
        primediv.first = q;
        primediv.second = 1;
        factor.push_back(primediv);
    }

    std::vector<long> gfactor;
    vec_ZZ gfactorZZ;
    gfactorZZ.SetLength(dim);
#ifdef DEBUG_LOG
    std::cerr<<"G is: "<<g<<std::endl;
#endif
    if(rfb.factor(gfactor, g))
    {
        for(long i=0; i<gfactor.size(); ++i) {
            gfactorZZ[i] = conv<ZZ>(gfactor[i]);
        }
    }
    else
    {
        std::cerr<<"G is not smooth!\n";
        return ZZ(-1);
    }
#ifdef DEBUG_LOG
    std::cerr<<"gfactor's results are: ";
    for(int i=0; i<gfactor.size(); ++i)
    {
        std::cerr<<gfactor[i]<<" ";
    }
    std::cerr<<"\ngfactorZZ: "<<gfactorZZ<<std::endl;
#endif
    std::vector<std::pair<ZZ, ZZ>> chineseinput;
    std::pair<ZZ, ZZ> sol;
    for(long i=0; i<factor.size(); ++i)
    {
        q = factor[i].first;
        if(q == ZZ(2))
        {
            continue;
        }
        ZZ_pPush push;
        ZZ_p::init(q);
        mat_ZZ_p matrixl;
        matrixl = conv<mat_ZZ_p>((*sieveres));
        for(long j=0; j<pairs.size(); ++j)
        {
            vec_ZZ schirmap = schirokauer_map(pairs[j].first, pairs[j].second, q, f);
            for(long k=0; k < f.d; ++k)
            {
                matrixl[k+dim - conv<long>(f.d)][j+1] = conv<ZZ_p>(schirmap[k]);
            }
        }
        ZZ_p det;
        vec_ZZ_p x;
        vec_ZZ_p gfactorZZP = -conv<vec_ZZ_p>(gfactorZZ);
#ifdef DEBUG_LOG
        std::cerr<<"Matrix is:\n"<<matrixl<<"\ngfactorZZP is: "<<gfactorZZP<<std::endl;
#endif
        mat_ZZ_p matrixlb;
        matrixlb.SetDims(matrixl.NumRows(), matrixl.NumCols()+1);
        for(unsigned long j=0; j<matrixl.NumRows(); ++j)
        {
            for(unsigned long k=0; k<matrixl.NumCols(); ++k)
            {
                matrixlb[j][k] = matrixl[j][k];
            }
        }
        for(unsigned long j=0; j<gfactorZZP.length(); ++j)
        {
            matrixlb[j][matrixlb.NumCols()-1] = gfactorZZP[j];
        }

        x = gausselim(matrixlb);
#ifdef DEBUG_LOG
        std::cerr<<"Solution is "<<x<<std::endl;
#endif
        sol.first = conv<ZZ>(x[0]);
        sol.second = q;
        chineseinput.push_back(sol);
    }
    sol.first = 0;
    sol.second = 2;
    chineseinput.push_back(sol);
    ZZ ans = (-chineserem(chineseinput)) % (ZZ_p::modulus()-1);
#ifdef DEBUG_LOG
    std::cerr<<"Chineserem answer for 0 mod 2 is: "<<ans<<std::endl;
#endif
    ZZ_p tp = conv<ZZ_p>(t), gp = conv<ZZ_p>(g);
    ZZ_p check = power(tp, ans);
#ifdef DEBUG_LOG
    std::cerr<<"Check is: "<<check<<std::endl;
#endif
    if(check != gp)
    {
        chineseinput.pop_back();
        sol.first = 1;
        sol.second = 2;
        chineseinput.push_back(sol);
        ans = (-chineserem(chineseinput)) % (ZZ_p::modulus()-1);
#ifdef DEBUG_LOG
        std::cerr<<"Chineserem answer for 1 mod 2 is: "<<ans<<std::endl;
#endif
        check = power(tp, ans);
#ifdef DEBUG_LOG
        std::cerr<<"Check is: "<<check<<std::endl;
#endif
        if (check != gp)
        {
            std::cerr<<"Wrong answer!"<<std::endl;
            return ZZ(-1);
        }
    }
    return ans;
}