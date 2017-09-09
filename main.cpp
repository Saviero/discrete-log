#include "NumberFieldSieve.h"
#include <string>
#include <cstring>
using namespace NTL;

int main()
{
    std::cout<<"Discrete logarithm using NFS.\nType 'exit' to exit.\n";
    std::cout<<"Enter field characteristic (integer prime number):\n";
    std::string pstr;
    std::cin>>pstr;
    if(pstr == "exit")
    {
        return 0;
    }
    ZZ p = NextPrime(conv<ZZ>(pstr.c_str()));
    ZZ_p::init(p);
    Logarithm logarithm;
    std::string t;
    std::string g;
    while(1) {
        std::cout<<"Enter t:\n";
        std::cin>>t;
        if(t == "exit")
        {
            break;
        }
        std::cout<<"Enter g:\n";
        std::cin>>g;
        if(g == "exit")
        {
            break;
        }
        std::cout<<"The answer is: "<<logarithm.log(conv<ZZ>(t.c_str()), conv<ZZ>(g.c_str()))<<"\n";
    }
    return 0;
}