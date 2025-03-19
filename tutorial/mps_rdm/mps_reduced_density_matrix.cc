#include "itensor/all.h"

using namespace itensor;

int main(){
    //Given an MPS called "psi",
//and assuming j > i
int N = 20;
auto sites = SpinHalf(N);
auto state = InitState(sites,"Up");
auto psi = randomMPS(state);

auto i = 4;
auto j = 10;

//'gauge' the MPS to site i
//any 'position' between i and j, inclusive, would work here
psi.position(i); 

auto psidag = dag(psi);
psidag.prime("Link");

//index linking i to i-1:
auto li_1 = leftLinkIndex(psi,i);

auto rho = prime(psi(i),li_1)*prime(psidag(i),"Site");
for(int k = i+1; k < j; ++k)
    {
    rho *= psi(k);
    rho *= psidag(k);
    }
//index linking j to j+1:
auto lj = rightLinkIndex(psi,j);

rho *= prime(psi(j),lj);
rho *= prime(psidag(j),"Site");

print(rho);

return 0;
}