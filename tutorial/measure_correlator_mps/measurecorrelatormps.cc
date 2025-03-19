#include "itensor/all.h"

using namespace itensor;

int
main(){
    //Number of sites
auto N = 100;

//Measure a correlation function
//at sites i and j
//Below we will assume j > i
auto i = 45;
auto j = 55;

auto sites = SpinHalf(N);

//Make a random MPS for testing
auto state = InitState(sites,"Up");
auto psi = randomMPS(state);

//Make the operators you want to measure
auto op_i = op(sites,"Sx",i);
auto op_j = op(sites,"Sz",j);

//'gauge' the MPS to site i
//any 'position' between i and j, inclusive, would work here
psi.position(i); 

//Create the bra/dual version of the MPS psi
auto psidag = dag(psi);

//Prime the link indices to make them distinct from
//the original ket links
psidag.prime("Link");

//index linking i-1 to i:
auto li_1 = leftLinkIndex(psi,i);

auto C = prime(psi(i),li_1)*op_i;
C *= prime(psidag(i),"Site");
for(int k = i+1; k < j; ++k)
    {
    C *= psi(k);
    C *= psidag(k);
    }
//index linking j to j+1:
auto lj = rightLinkIndex(psi,j);

C *= prime(psi(j),lj)*op_j;
C *= prime(psidag(j),"Site");

auto result = elt(C); //or eltC(C) if expecting complex

print(result);

return 0;

}
