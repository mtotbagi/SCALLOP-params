#pragma once

#include <cassert>
#include <optional>
#include <vector>

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include "ec.hpp"
#include "ecp.hpp"

class isog
{
private:
    std::shared_ptr<const ec> domain;
    std::shared_ptr<const ec> codomain;
    int _degree;
    NTL::ZZ_pEX psi, psi_der, psi_der_der, psi_der_der_der;
    NTL::ZZ_pE t; // Useful for degree 2 isogenies

public:
    isog(std::shared_ptr<const ec> E0, std::shared_ptr<const ec> E1, int degree, NTL::ZZ_pEX const &psi, NTL::ZZ_pEX const &psi_der, NTL::ZZ_pEX const &psi_der_der, NTL::ZZ_pEX const &psi_der_der_der) :
        domain(E0),
        codomain(E1),
        _degree(degree),
        psi(psi),
        psi_der(psi_der),
        psi_der_der(psi_der_der),
        psi_der_der_der(psi_der_der_der) {};

    isog(std::shared_ptr<const ec> E0, std::shared_ptr<const ec> E1, int degree, NTL::ZZ_pEX const &psi, NTL::ZZ_pE const &t) :
        domain(E0),
        codomain(E1),
        _degree(degree),
        psi(psi),
        t(t) {};

    ec const &get_domain() const { return *domain; }
    std::shared_ptr<const ec> const &get_domain_ptr() const { return domain; }
    ec const &get_codomain() const { return *codomain; }
    std::shared_ptr<const ec> const &get_codomain_ptr() const { return codomain; }

    ecp operator()(ecp const &P) const;
    ecp _even_evaluation(ecp const &P) const;
};

isog isogeny(ecp const &K, int degree);

NTL::ZZ_pEX kernel_polynomial(ecp const &K, int degree);

class isog_2_chain
{
    private:
        std::vector<isog> isog_2_parts;
    
    public:
        isog_2_chain(ecp const &P, std::vector<int> const &strategy);
        ecp operator()(ecp const &P) const {
            ecp Pi = P;
            for (isog phi : isog_2_parts) {
                Pi = phi(Pi);
            }
            return Pi;
        };
        ec const &get_domain() const { return isog_2_parts.front().get_domain(); }
        ec const &get_codomain() const { return isog_2_parts.back().get_codomain(); }
};

class isomorphism
{
    private:
        NTL::ZZ_pE u, r, s, t;
        std::shared_ptr<const ec> domain;
        std::shared_ptr<const ec> codomain;

    public:
        isomorphism(std::shared_ptr<const ec> E, std::shared_ptr<const ec> F);
        ec const &get_domain() const { return *domain; }
        ec const &get_codomain() const { return *codomain; }
        
        ecp operator()(ecp const &P) const;
};