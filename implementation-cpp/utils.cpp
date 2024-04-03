#include <cassert>
#include <unordered_map>
#include <vector>
#include <cassert>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>

#include "ec.hpp"
#include "ecp.hpp"
#include "utils.hpp"


int _BSGS(ecp const &Q, ecp const &base, int base_order) {

    assert (Q.curve() == base.curve());

    size_t m = NTL::SqrRoot(base_order);

    std::unordered_map<ecp, size_t> baby_steps;

    ecp G = 0*base;

    baby_steps[G] = 0;


    for (size_t i = 1; i <= m+1; i++) {
        G += base;
        baby_steps[G] = i;
    }

    ecp min_G = -G;
    ecp gamma = Q;

    for (size_t i = 0; i < m; i++) {
        auto hit = baby_steps.find(gamma);

        if (hit != baby_steps.end())
            return hit->second + m*i;
        gamma += min_G;
    }

    return -1;
}

NTL::ZZ DLP(ecp const &Q, ecp const &base, int ell, int e) {
    ecp gamma = NTL::power(NTL::ZZ(ell), e-1)*base;

    NTL::ZZ ell_to_k(1);
    NTL::ZZ x(0);
    for (int k = 0; k < e; k++) {
        ecp H = NTL::power(NTL::ZZ(ell), e-1-k)*(Q - x*base);
        int d = _BSGS(H, gamma, ell);
        if (d == -1) {
            return NTL::ZZ(-1);
        }
        x += ell_to_k*d;
        ell_to_k *= ell;
    }
    return x;
}

std::unordered_map<int, int> factor(Integer const &N) {
    // Trial division for now
    int p = 1;

    std::unordered_map<int, int> factored;
    auto Ni = N;
    while (Ni > 1) {
        p = NTL::NextPrime(p+1);
        bool first = true;
        while ((Ni % p) == 0) {
            Ni /= p;
            if (first) {
                factored[p] = 1;
                first = false;
            } else {
                factored[p] += 1;
            }
        }
    }

    return factored;
}

size_t _myroots(NTL::ZZ_pE *roots, NTL::ZZ_pEX const &f)
{
    //XXX This is only marginally faster than NTL::CanZass().
    //XXX Is it worth the extra (potentially buggy) code?
    assert(NTL::IsOne(NTL::LeadCoeff(f)));
    NTL::ZZ_pEX x; NTL::SetX(x);
    auto sfd = NTL::SquareFreeDecomp(f);
    size_t idx = 0;
    for (auto const &[g,m]: sfd) {
        NTL::ZZ_pEXModulus F;
        NTL::build(F, g);
        auto h = NTL::PowerXMod(NTL::ZZ_pE::cardinality(), F);
        auto g1 = NTL::GCD(g, h - x);
        for (auto const &r: NTL::FindRoots(g1))
            for (ssize_t i = 0; i < m; ++i)
                roots[idx++] = r;
    }
    return idx;
}

std::optional<NTL::ZZ_pE> sqrt(NTL::ZZ_pE const &alpha) 
// Take EXTREME caution to only call this when ZZ_pE agrees with alpha 
// This is the same as with all arithmetic in field extensions...
{
    NTL::ZZ_pEX f;
    NTL::ZZ_pE roots[2];
    SetCoeff(f, 2);
    f[0] = -alpha;
    if (_myroots(roots, f) == 0) {
        return {};
    }
    return roots[0];
}