#include "NumberFieldSieve.h"
//#define DEBUG_LOG
//#define DEBUG_LANCZOS
//#define DEBUG
//#define DEBUG_SCHIR

ZZ_p scalmul(vec_ZZ_p x, vec_ZZ_p y)
{
    if(x.length() != y.length())
    {
        std::cerr<<"Length of vectors are not the same!\n";
        return ZZ_p(0);
    }
    ZZ_p result = ZZ_p(0);
    for(unsigned long i = 0; i<x.length(); ++i)
    {
        result += x[i]*y[i];
    }
#ifdef DEBUG_LANCZOS
    std::cerr<<"Scalar mul for x="<<x<<", y="<<y<<"is "<<result<<std::endl;
#endif
    return result;
}

ZZ_p scalmul_mod(vec_ZZ_p x, vec_ZZ_p y, mat_ZZ_p a)
{
    if(x.length() != a.NumCols())
    {
        std::cerr<<"Invalid matrix's dims!\n";
        return ZZ_p(0);
    }
    vec_ZZ_p xnew;
    mul(xnew, a, x);
    std::cerr<<"Scalar mul_mod for x="<<x<<", y="<<y<<"; new x is "<<xnew<<std::endl;
    return scalmul(xnew, y);
}

vec_ZZ_p lanczos(mat_ZZ_p A, vec_ZZ_p b)
{
    if(IsZero(b))
    {
#ifdef DEBUG_LANCZOS
        std::cerr<<"B is zero vector; getting the last column of A as a new B"<<std::endl;
#endif
        b = -A[A.NumCols()-1];
#ifdef DEBUG_LANCZOS
        std::cerr<<"B now is"<<b<<std::endl;
#endif
        mat_ZZ_p Atemp = A;
        A.SetDims(A.NumRows(), A.NumCols()-1);
        for(unsigned long i=0; i<A.NumRows(); ++i)
        {
            for(unsigned long j=0; j<A.NumCols(); ++j)
            {
                A[i][j] = Atemp[i][j];
            }
        }
#ifdef DEBUG_LANCZOS
        std::cerr<<"A now is"<<A<<std::endl;
#endif
    }
    vec_ZZ_p answer;
    mat_ZZ_p D;
    D.SetDims(A.NumRows(), A.NumRows());
    for(unsigned long i=0; i<D.NumCols(); ++i)
    {
        ZZ_p rand = random_ZZ_p();
        while(rand == 0) {rand = random_ZZ_p();}
        D[i][i] = rand;
    }
#ifdef DEBUG_LANCZOS
    std::cerr<<"Matrix D is\n"<<D<<std::endl;
#endif
    mat_ZZ_p Atemp, Anew;
    mul(Atemp, transpose(A), power(D, 2));
    vec_ZZ_p bnew;
    mul(bnew, b, Atemp);
    mul(Anew, Atemp, A);
#ifdef DEBUG_LANCZOS
    std::cerr<<"New A is\n"<<Anew<<"\nNew B is "<<bnew<<std::endl;
#endif
    std::vector<vec_ZZ_p> omega;
    vec_ZZ_p curromega = bnew;
    omega.push_back(curromega);
#ifdef DEBUG_LANCZOS
    std::cerr<<"Omega[0] is "<<curromega<<std::endl;
#endif
    mul(curromega, Anew, omega[0]);
    ZZ_p coeff1 = scalmul_mod(omega[0], curromega, Anew)/scalmul_mod(omega[0], omega[0], Anew);
    ZZ_p coeff2;
    vec_ZZ_p aomeg1;
    vec_ZZ_p aomeg2;
    curromega -= coeff1*omega[0];
    omega.push_back(curromega);
#ifdef DEBUG_LANCZOS
    std::cerr<<"Omega[1] is "<<curromega<<std::endl;
#endif
    unsigned long i = 2;
    while(!IsZero(curromega))
    {
        if(scalmul_mod(curromega, curromega, Anew) == 0)
        {
            std::cerr<<"Cannot solve system!"<<std::endl;
            return answer;
        }
        mul(curromega, Anew, omega[i-1]);
        mul(aomeg1, Anew, omega[i-1]);
        mul(aomeg2, Anew, omega[i-2]);
        coeff1 = -(scalmul_mod(omega[i-1], aomeg1, Anew)/scalmul_mod(omega[i-1], omega[i-1], Anew));
        coeff2 = -(scalmul_mod(omega[i-1], aomeg2, Anew)/scalmul_mod(omega[i-2], aomeg2, Anew));
        curromega += coeff1*omega[i-1];
        curromega += coeff2*omega[i-2];
        omega.push_back(curromega);
#ifdef DEBUG_LANCZOS
        std::cerr<<"Omega["<<i<<"] is "<<curromega<<std::endl;
#endif
        ++i;
    }
    vec_ZZ_p result;
    result.SetLength(A.NumRows());
    for(unsigned long j=0; j < omega.size()-1; ++j)
    {
        result += (scalmul(omega[j], bnew)/scalmul(omega[j], omega[j]))*omega[j];
    }
    vec_ZZ_p check = A*result;
    if(check != b)
    {
#ifdef DEBUG_LANCZOS
        std::cerr<<"Wrong result! Result is: "<<result<<std::endl;
        return lanczos(A, b);
#endif
    }
    return result;
}

vec_ZZ schirokauer_map(const ZZ& a, const ZZ& b, const ZZ& l, Polynomial f)
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

mat_ZZ* sieve(std::vector<std::pair<ZZ, ZZ>>& s, const Polynomial& f, const AlgebraicFactorBase& base, const ZZ& v)
{
#ifdef DEBUG
    std::cerr<<"Bound is: "<<v<<"\n";
#endif
    // Preparing output matrix
    mat_ZZ* res = new mat_ZZ();
    res->SetDims(base.getTotalSize()+conv<long>(f.d), base.getTotalSize()+conv<long>(f.d));
    long numOfRows = 1;
    ZZ newv = v;
    ZZ oldv = ZZ(1);
    while(numOfRows < base.getTotalSize()+conv<long>(f.d)) {
        for (ZZ a = oldv; a < newv; ++a) {
            if (numOfRows >= base.getTotalSize() + conv<long>(f.d)) {
                break;
            }
            for (ZZ b = -newv; b < newv; ++b) {
                if (b == 0) {
                    continue;
                }
                if (numOfRows >= base.getTotalSize() + conv<long>(f.d)) {
                    break;
                }
                if (GCD(a, b) == ZZ(1)) // we need only coprime a and b
                {
                    std::vector<long> expvec;
#ifdef DEBUG
                    std::cerr<<"Trying to factor: ("<<a<<"; "<<b<<")\n";
#endif
                    if (base.factor(expvec, a, b)) // factoring a+bm and a+b*alpha
                    {
#ifdef DEBUG
                        std::cerr<<"Factored\n";
#endif
                        std::pair<ZZ, ZZ> ab;
                        ab.first = a;
                        ab.second = b;
                        s.push_back(ab);
                        for (long i = 0; i < base.getTotalSize(); ++i) {
                            (*res)[i][numOfRows] = expvec[i];
                        }
                        ++numOfRows;
                    } else {
#ifdef DEBUG
                        std::cerr<<"Not smooth\n";
#endif
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
#ifdef DEBUG_LOG
    std::cerr<<"Polynomial is "<<f.f<<"; M is "<<f.m<<std::endl;
#endif
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
            (*sieveres)[i][0] = factort.at(i);
        }
        for(unsigned long i=base.fb.r.size(); i < base.getTotalSize()+f.d; ++i)
        {
            (*sieveres)[i][0] = 0;
        }
    }
    else
    {
        std::cerr<<"T is not smooth!\n";
        return ZZ(-1);
    }

    ZZ q = p-1;
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
    gfactorZZ.SetLength(base.getTotalSize()+conv<long>(f.d));
#ifdef DEBUG_LOG
    std::cerr<<"G is: "<<g<<std::endl;
#endif
    if(base.fb.factor(gfactor, g))
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
                matrixl[k+base.getTotalSize()][j+1] = conv<ZZ_p>(schirmap[k]);
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
        long rank = gauss(matrixlb);
#ifdef DEBUG_LOG
        std::cerr<<"Echeloned matrix with B is: \n"<<matrixlb<<"Rank is: "<<rank<<std::endl;
#endif
        matrixl.SetDims(rank, rank);
        gfactorZZP.SetLength(rank);
        for(unsigned long j=0; j<rank; ++j)
        {
            for(unsigned long k=0; k<rank; ++k)
            {
                matrixl[j][k] = matrixlb[j][k];
            }
        }
        for(unsigned long j=0; j<rank; ++j)
        {
            gfactorZZP[j] = matrixlb[j][matrixlb.NumCols()-1];
        }
        for(unsigned long j=rank; j<matrixlb.NumRows(); ++j)
        {
            if(matrixlb[j][matrixlb.NumCols()-1] != 0)
            {
                std::cerr<<"Error! No solution to matrix"<<std::endl;
                return ZZ(-1);
            }
        }
#ifdef DEBUG_LOG
        std::cerr<<"Matrix is:\n"<<matrixl<<"gfactorZZP is: "<<gfactorZZP<<std::endl;
#endif
        matrixl = transpose(matrixl);
        solve(det, x, matrixl, gfactorZZP);
        std::pair<ZZ, ZZ> sol;
#ifdef DEBUG_LOG
        std::cerr<<"Determinant is "<<det<<std::endl;
        std::cerr<<"Solution is "<<x<<std::endl;
#endif
        sol.first = conv<ZZ>(x[0]);
        sol.second = q;
        chineseinput.push_back(sol);
    }
    std::pair<ZZ, ZZ> sol;
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
        ans = -(-chineserem(chineseinput)) % (ZZ_p::modulus()-1);
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