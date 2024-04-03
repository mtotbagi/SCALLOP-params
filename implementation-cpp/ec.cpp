#include <cassert>
#include <optional>

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ.h>
#include "ec.hpp"
#include "ecp.hpp"
#include "utils.hpp"

ec::ec(NTL::ZZ_pE const &a, NTL::ZZ_pE const &b) : _a{a}, _b{b}
{
    if ((4*this->_a*this->_a*this->_a + 27*this->_b*this->_b) == 0)
        throw std::logic_error("curve is singular");
}

ecp ec::operator()(NTL::ZZ_pE new_x, NTL::ZZ_pE new_y) const {
    return ecp(std::make_shared<const ec>(*this), new_x, new_y);
}

ecp ec::identity() const {
    return ecp(std::make_shared<const ec>(*this));
}

std::optional<ecp> ec::lift_x( NTL::ZZ_pE const &x) const
{
    NTL::ZZ_pE rhs = (x*x + this->_a)*x + this->_b;
    auto y = sqrt(rhs);
    if (y)
        return ecp(std::make_shared<const ec>(*this), x, *y);
    return {};
}

ecp ec::random_point() const
{
    while (true) {
        auto pt = this->lift_x(NTL::random_ZZ_pE());
        if (pt) {
            return *pt;
        }
    }
}

NTL::ZZ_pE ec::j_invariant() const
{
    auto &a = _a, &b = _b;
    auto a3 = NTL::ZZ_pE(4)*a*a*a;
    auto b2 = NTL::ZZ_pE(27)*b*b;
    return NTL::ZZ_pE(1728) * a3/(a3 + b2);
}

ec ec::from_j(NTL::ZZ_pE const &j)
{
    if (j == 0) { return ec(NTL::ZZ_pE(0), NTL::ZZ_pE(1));}
    NTL::ZZ_pE const f(1728);
    if (j == f) { return ec(NTL::ZZ_pE(1), NTL::ZZ_pE(0));}
    NTL::ZZ_pE two(2), three(3);
    auto j2 = j*j;
    auto j3 = j2*j;
    NTL::ZZ_pE a = three * (f*j - j2);
    NTL::ZZ_pE b = two * (-j3 + two*f*j2 - f*f*j);
    auto E = ec(a, b);
    assert(E.j_invariant() == j);
    return E;
}

ecp ec::random_point_of_order(NTL::ZZ cof, int ell, int e) const {
    ecp P = this->random_point();
    P = cof*P;
    while (!(NTL::power(NTL::ZZ(ell), e-1)*P)) {
        P = this->random_point();
        P = cof*P;
    }

    assert (NTL::power(NTL::ZZ(ell), e-1)*P);
    assert (!(NTL::power(NTL::ZZ(ell), e)*P));
    return P;
}

ecp ec::random_point_of_order_dividing(NTL::ZZ D) const {

    NTL::ZZ p = NTL::ZZ_p::modulus();
    NTL::ZZ cof = p + 1;
    assert (cof % D == 0);
    cof /= D;

    ecp P = this->random_point();
    P = cof*P;

    assert (!(D*P));
    return P;
}

std::pair<ecp, ecp> const ec::torsionBasis(int ell, int e) const
{
    // What if we are on the twist? 
    NTL::ZZ p = NTL::ZZ_p::modulus();
    NTL::ZZ cof = p + 1;

    NTL::ZZ ellcof = NTL::power(NTL::ZZ(ell), e-1);
    NTL::ZZ le = ellcof * ell;
    assert (cof % le == 0);
    cof /= le;

    ecp P = this->random_point_of_order(cof, ell, e);
    auto ellP = ellcof*P;

    while (true) {
        ecp Q = this->random_point_of_order(cof, ell, e);
        if (DLP(ellcof*Q, ellP, ell, 1) == NTL::ZZ(-1))
            return {P, Q};
    }
}

ecp const ec::completeBasis(int ell, int e, ecp const &P) const
{
    assert (ell == 2);
    NTL::ZZ p = NTL::ZZ_p::modulus();
    NTL::ZZ cof = p + 1;

    NTL::ZZ ellcof = NTL::power(NTL::ZZ(ell), e-1);
    NTL::ZZ le = ellcof * ell;
    assert (cof % le == 0);
    cof /= le;

    ecp ellP = ellcof*P;
    while (true) {
        ecp Q = this->random_point_of_order(cof, ell, e);
        if (DLP(ellcof*Q, ellP, ell, 1) == NTL::ZZ(-1))
            return Q;
    }
}