#include "itensor/all.h"

using namespace itensor;

int main(){
    int N = 20;
auto sites = SpinHalf(N);
auto state = InitState(sites,"Up");
auto psi = randomMPS(state);
auto phi = randomMPS(state);

int i = 4;
int j = 10;
auto op_i = op(sites,"Sz",i);
auto op_j = op(sites,"Sz",j);

auto phidag = dag(phi);

auto M = psi(1)*phidag(1);
for(auto n : range1(2,N))
    {
    M *= psi(n);
    if(n == i)
        {
        M *= op_i*prime(phidag(i),"Site");
        }
    else if(n == j)
        {
        M *= op_j*prime(phidag(j),"Site");
        }
    else
        {
        M *= phidag(n);
        }
    }
auto result = elt(M);

print(result);

return 0;
}