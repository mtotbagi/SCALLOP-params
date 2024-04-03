#include <cassert>
#include <optional>

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include "ec.hpp"
#include "ecp.hpp"
#include "isog.hpp"

//TODO rework, this is slow from Field extension code

NTL::ZZ_pEX _derivative(NTL::ZZ_pEX const &F)
{
    NTL::ZZ_pEX f;
    if (NTL::IsZero(F)) {
        f = 0;
        return f;
    }
    size_t deg = NTL::deg(F);
    if (NTL::deg(F) == 0) {
        f = 0;
        return f;
    }
    NTL::SetCoeff(f, deg-1);

    for (size_t i = 1; i <= deg; i++) {
        NTL::ZZ_pE coeff = NTL::coeff(F, i);
        f[i-1] = i*coeff;
    }

    return f;
}

isog isogeny(ecp const &K, int degree) {

    NTL::ZZ_pEX h_K;
    if (degree > 2) {
        h_K = kernel_polynomial(K, degree);
    } else {
        NTL::SetCoeff(h_K, 1);
        h_K[0] = -K.aff_x();
    }
    if (degree > 2) {
        NTL::ZZ_pEX psi_der = _derivative(h_K);
        NTL::ZZ_pEX psi_der_der = _derivative(psi_der);
        NTL::ZZ_pEX psi_der_der_der = _derivative(psi_der_der);
        int n = (degree - 1)/2;
        assert (n == NTL::deg(h_K));

        NTL::ZZ_pE s1 = -NTL::coeff(h_K, n-1), s2 = NTL::coeff(h_K, n-2), s3 = -NTL::coeff(h_K, n-3);
        NTL::ZZ_pE a = K.curve().a(), b = K.curve().b();

        auto s1_sq = s1*s1;
        NTL::ZZ_pE t = 6*(s1_sq - 2*s2) + n*2*a;
        NTL::ZZ_pE w = 10*(s1_sq*s1 - 3*s1*s2 + 3*s3) + 6*a*s1 + n*4*b;

        ec E1(a - 5*t, b - 7*w);

        return isog(K.curve_ptr(), std::make_shared<const ec>(E1), degree, h_K, psi_der, psi_der_der, psi_der_der_der);
    } else {
        NTL::ZZ_pE a = K.curve().a(), b = K.curve().b();

        NTL::ZZ_pE x0 = -coeff(h_K, 0);
        NTL::ZZ_pE t = 3*x0*x0 + a;
        NTL::ZZ_pE w = x0*t;

        ec E1(a - 5*t, b - 7*w);
        return isog(K.curve_ptr(), std::make_shared<const ec>(E1), degree, h_K, t);
    }
}

std::pair<NTL::ZZ_pE, NTL::ZZ_pE> _prod_derivative_pair(std::pair<NTL::ZZ_pE, NTL::ZZ_pE> const &f1, std::pair<NTL::ZZ_pE, NTL::ZZ_pE> const &f2) {
    std::pair<NTL::ZZ_pE, NTL::ZZ_pE> out(f1.first*f2.first, f1.first*f2.second + f1.second*f2.first);
    return out;
}

ecp isog::operator()(ecp const &P) const
{
    assert (P.curve() == this->get_domain());

    if (this->_degree == 2) {
        return this->_even_evaluation(P);
    }

    NTL::ZZ_pE x = P.aff_x();

    NTL::ZZ_pEX psi = this->psi;
    NTL::ZZ_pEX psi_der = this->psi_der;
    NTL::ZZ_pEX psi_der_der = this->psi_der_der;
    NTL::ZZ_pEX psi_der_der_der = this->psi_der_der_der;

    NTL::ZZ_pE psi_x = NTL::coeff(psi, 0);
    NTL::ZZ_pE psi_der_x = NTL::coeff(psi_der, 0);
    NTL::ZZ_pE psi_der_der_x = NTL::coeff(psi_der_der, 0);
    NTL::ZZ_pE psi_der_der_der_x = NTL::coeff(psi_der_der_der, 0);

    size_t n = NTL::deg(psi);
    NTL::ZZ_pE xi(1);

    for (size_t i = 1; i <= n; i++) { //NTL::coeff returns zero when i > degree
        xi *= x;
        psi_x += xi*NTL::coeff(psi, i);
        psi_der_x += xi*NTL::coeff(psi_der, i);
        psi_der_der_x += xi*NTL::coeff(psi_der_der, i);
        psi_der_der_der_x += xi*NTL::coeff(psi_der_der_der, i);
    }

    // Evaluating the derivative without expanding the polynomial Pogchamp
    std::pair<NTL::ZZ_pE, NTL::ZZ_pE> psi_x_pair(psi_x, psi_der_x);
    std::pair<NTL::ZZ_pE, NTL::ZZ_pE> psi_der_x_pair(psi_der_x, psi_der_der_x);
    std::pair<NTL::ZZ_pE, NTL::ZZ_pE> psi_der_der_x_pair(psi_der_der_x, psi_der_der_der_x);

    NTL::ZZ_pE b4 = this->get_domain().a()*2, b6 = this->get_domain().b()*4;
    std::pair<NTL::ZZ_pE, NTL::ZZ_pE> aux_1_pair((4*x*x*x + 2*b4*x + b6), (12*x*x + 2*b4));
    std::pair<NTL::ZZ_pE, NTL::ZZ_pE> aux_2_pair((6*x*x + b4), (12*x));
    std::pair<NTL::ZZ_pE, NTL::ZZ_pE> aux_3_pair((this->_degree*x + 2*NTL::coeff(psi, n-1)), this->_degree);

    //begin the evaluation
    auto psi_x_sq_pair = _prod_derivative_pair(psi_x_pair, psi_x_pair);
    auto psi_der_x_sq_pair = _prod_derivative_pair(psi_der_x_pair, psi_der_x_pair);
    auto psi_x_psi_der_der_x_pair = _prod_derivative_pair(psi_x_pair, psi_der_der_x_pair);

    std::pair<NTL::ZZ_pE, NTL::ZZ_pE> temp(psi_der_x_sq_pair.first - psi_x_psi_der_der_x_pair.first, psi_der_x_sq_pair.second - psi_x_psi_der_der_x_pair.second);
    auto summand_1 = _prod_derivative_pair(aux_1_pair, temp);

    auto psi_x_psi_der_x_pair = _prod_derivative_pair(psi_der_x_pair, psi_x_pair);
    auto summand_2 = _prod_derivative_pair(psi_x_psi_der_x_pair, aux_2_pair);

    auto summand_3 = _prod_derivative_pair(psi_x_sq_pair, aux_3_pair);
    
    std::pair<NTL::ZZ_pE, NTL::ZZ_pE> N_x_pair(summand_1.first - summand_2.first + summand_3.first, summand_1.second - summand_2.second + summand_3.second);

    NTL::ZZ_pE one(1);
    if (psi_x_sq_pair.first == 0) {
        return this->get_codomain().identity();
    }
    std::pair<NTL::ZZ_pE, NTL::ZZ_pE> psi_x_sq_inv_pair(one/psi_x_sq_pair.first, -psi_x_sq_pair.second/(psi_x_sq_pair.first*psi_x_sq_pair.first));
    auto eval_res = _prod_derivative_pair(N_x_pair, psi_x_sq_inv_pair);

    auto new_x = eval_res.first;
    auto new_y = P.aff_y()*eval_res.second;

    auto phi_P = this->get_codomain()(new_x, new_y);

    return phi_P;
}

ecp isog::_even_evaluation(ecp const &P) const
{
    NTL::ZZ_pE x0 = -NTL::coeff(this->psi, 0), x = P.aff_x(), y = P.aff_y();
    NTL::ZZ_pE xmx0 = x-x0;
    NTL::ZZ_pE new_x = x + t/xmx0;
    NTL::ZZ_pE new_y = y*(1 - t/(xmx0*xmx0));
 
    ecp phi_P(this->get_codomain_ptr(), new_x, new_y);
    return phi_P;
}

NTL::ZZ_pEX kernel_polynomial(ecp const&K, int degree) {
    /* 
    K point of order degree, 
    returns the kernel poly of <K>
    */

    NTL::ZZ_pEX Psi;

    NTL::ZZ_pE x = K.aff_x();
    if (degree == 3) {
        NTL::SetCoeff(Psi, 1);
        Psi[0] = -x;  
        return Psi; 
    }

    NTL::ZZ_pEX X;
    NTL::SetX(X);

    NTL::ZZ_pEX f0 = X - x;

    size_t m = (degree-1)/2;

    std::vector<NTL::ZZ_pEX> fs;
    fs.reserve(m);
    
    fs.push_back(f0);
    ecp Ki = K;

    for (size_t i = 1; i < m; i++) {
        Ki = K + Ki;
        fs.push_back(X - Ki.aff_x());
    }

    Psi = product_tree<NTL::ZZ_pEX>(fs);

    return Psi;
}

isog_2_chain::isog_2_chain(ecp const &Q, std::vector<int> const &strategy)
{
    size_t MAX = strategy.size();
    std::vector<ecp> pts;
    pts.reserve(MAX);
    std::vector<size_t> pts_index;
    pts_index.reserve(MAX);
	size_t npts = 0;
	size_t ii = 0;

    int ell = 2;

    ecp P = Q;
    size_t index = 0;

    for (size_t row = 1; row < MAX; row++) {
        while (index < (MAX - row)) {
            if (!(row == 1)) {
                pts[npts] = P;
                pts_index[npts] = index;
            } else {
                pts.push_back(P);
                pts_index.push_back(index);
            }
            npts += 1;
            int m = strategy[ii];
            ii += 1;
            P = NTL::power(NTL::ZZ(2), m)*P;
            index += m;
        }

        if (!(row == 1)) {
            isog phi = isogeny(P, ell);
            isog_2_parts.push_back(phi);
            for (size_t i = 0; i < npts; i++) {
                pts[i] = phi(pts[i]);
            }

        } else {
            //four isogenies for some reason
            isog phi = isogeny(2*P, ell);
            isog_2_parts.push_back(phi);
            for (size_t i = 0; i < npts; i++) {
                pts[i] = phi(pts[i]);
            }
            isog phi2 = isogeny(phi(P), ell);
            isog_2_parts.push_back(phi2);
            for (size_t i = 0; i < npts; i++) {
                pts[i] = phi2(pts[i]);
            }
        }

        P = pts[npts-1];
        index = pts_index[npts - 1];
        npts -= 1;
    }

    isog phi = isogeny(P, 2);
    isog_2_parts.push_back(phi);
}

isomorphism::isomorphism(std::shared_ptr<const ec> EE, std::shared_ptr<const ec> FF) {
    domain = EE;
    codomain = FF;

    ec E = *EE;
    ec F = *FF;
    NTL::ZZ_pE j = E.j_invariant();
    assert (F.j_invariant() == j);
    //Probably neeeever gonna happen
    assert (!(j == 0));
    assert (!(j == 1728));

    NTL::ZZ_pE u_sqr = (E.a()*F.b())/(F.a()*E.b());
    u = *sqrt(u_sqr);

    // Maybe later we will use montgomery form, but for now these are 0
    s = 0;
    r = 0;
    t = 0;
}

ecp isomorphism::operator()(ecp const &P) const {
    NTL::ZZ_pE u2 = u*u;
    NTL::ZZ_pE u3 = u2*u;
    return this->get_codomain()(u2*P.aff_x(), u3*P.aff_y());
}