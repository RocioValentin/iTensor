#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include <chrono>

using namespace itensor;
using std::vector;

int main()
{
    int N = 20; //number of sites
    Real tstep = 1; //0.02; //time step (smaller is generally more accurate)
    Real ttotal_max = 10; //1.0; //total time to evolve
    Real cutoff = 1E-8;
// int N = 50; //number of sites
// Real tstep = 0.02; //time step (smaller is generally more accurate)
// Real ttotal_max = 2.0; //total time to evolve
// Real cutoff = 1E-8; //truncation error cutoff when restoring MPS form

//Define a site set object "sites" which lets us
//easily obtain Site indices defining our Hilbert space
//and S=1/2 single-site operators
auto sites = SpinHalf(N);

//Make initial MPS psi to be in the Neel state
auto state = InitState(sites);
for(auto j : range1(N))
    {
    state.set(j,j%2==1?"Up":"Dn");
    }
auto psi = MPS(state);

//Create a std::vector (dynamically sizeable array)
//to hold the Trotter gates
auto gates = vector<BondGate>();

//Create the gates exp(-i*tstep/2*hterm)
//and add them to gates
for(int b = 1; b <= N-1; ++b)
    {
    auto hterm = op(sites,"Sz",b)*op(sites,"Sz",b+1);
    hterm += 0.5*op(sites,"S+",b)*op(sites,"S-",b+1);
    hterm += 0.5*op(sites,"S-",b)*op(sites,"S+",b+1);

    auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
    gates.push_back(g);
    }
//Create the gates exp(-i*tstep/2*hterm) in reverse
//order (to get a second order Trotter breakup which
//does a time step of "tstep") and add them to gates
for(int b = N-1; b >= 1; --b)
    {
    auto hterm = op(sites,"Sz",b)*op(sites,"Sz",b+1);
    hterm += 0.5*op(sites,"S+",b)*op(sites,"S-",b+1);
    hterm += 0.5*op(sites,"S-",b)*op(sites,"S+",b+1);

    auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
    gates.push_back(g);
    }

//Save initial state;
auto psi0 = psi;

//Time evolve, overwriting psi when done
// gateTEvol(gates,ttotal,tstep,psi,{"Cutoff=",cutoff,"Verbose=",true});

// printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));

    // Create an output file to store the data
    std::ofstream outfile("maxLinkDim_vs_ttotal3.csv");
    outfile << "ttotal,maxLinkDim,frequency,delta\n";

    // Time evolve and record maxLinkDim at each time step
    for (Real ttotal = 1; ttotal <= ttotal_max; ttotal += 1)
    {
        auto start = std::chrono::high_resolution_clock::now();
        // Perform time evolution
        gateTEvol(gates, ttotal, tstep, psi, {"Cutoff=", cutoff, "Verbose=", true});

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration = end - start;
        double delta = duration.count();
        double frequency = 1.0 / delta;
        
        // Get the max link dimension and write it to the file
        int max_dim = maxLinkDim(psi);
        outfile << ttotal << "," << max_dim << "," << frequency << "," << delta << "\n";
    }

//Print overlap of final state with initial state
//(Will be complex so using innerC which can return complex);
Print(innerC(psi,psi0));

return 0;
}