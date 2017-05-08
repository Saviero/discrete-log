#include "FactorBaseGenerator.h"

FactorBaseGenerator::FactorBaseGenerator()
{
    countA = 0;
    countR = 0;
}

void FactorBaseGenerator::generate(ZZ bound, ZZ_pX f, ZZ_p m)
{
    std::ofstream outR("/fbR.txt", std::ofstream::out);
    std::ofstream outA("/fbA.txt", std::ofstream::out);
    ZZ q = ZZ(2);
    std::pair<ZZ, ZZ> elem;
    while(q<bound)
    {
        ZZ r = rep(m)%q;
        elem.first = q;
        elem.second = r;
        R.push_back(elem);
        outR<<elem.first<<" "<<elem.second<<"\n";
        countR++;
        std::cout<<elem.first<<" "<<elem.second<<"\n";
        for(ZZ rr = ZZ(0); rr<q ; ++rr)
        {
            ZZ_p x;
            ZZ_p rp = conv<ZZ_p>(rr);
            eval(x, f, rp);

            if(rep(x)%q == 0)
            {
                elem.second = rr;
                A.push_back(elem);
                outA<<elem.first<<" "<<elem.second<<"\n";
                std::cout<<elem.first<<" "<<elem.second<<"\n";
                countA++;
                break;
            }
        }
        std::cout<<"\n";
        q = NextPrime(q+1);
    }
    outA.close();
    outR.close();
}