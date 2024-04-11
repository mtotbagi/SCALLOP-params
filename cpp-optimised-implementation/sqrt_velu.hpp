#pragma once

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <vector>

#include "fp2.hpp"
#include "montgomery.hpp"

// Currently not working....
// Also, NTL resultant computation is hella slow??

std::vector<fp2_elem> Poly_mult(std::vector<fp2_elem> const &a, std::vector<fp2_elem> const &b) {
    std::vector<fp2_elem> ab;

    size_t m = a.size() - 1;
    size_t n = b.size() - 1;
    for (size_t i = 0 ; i <= m + n ; i++) {
        fp2_elem coeff = Fp2_zero();
        for (size_t j = 0 ; j <= std::min(i, m) ; j++) {
            if (i-j <= n) {
                coeff = Fp2_add(coeff, Fp2_mul(a[j], b[i-j]));
            }
        }
        ab.push_back(coeff);
    }

    return ab;
}

std::vector<fp2_elem> Poly_add(std::vector<fp2_elem> const &a, std::vector<fp2_elem> const &b) {
    std::vector<fp2_elem> a_sum_b;

    size_t m = a.size() - 1;
    size_t n = b.size() - 1;

    for (size_t i = 0 ; i <= std::min(m, n) ; i++) {
        a_sum_b.push_back(Fp2_add(a[i], b[i]));
    }
    if (m > n) {
        for (size_t i = n; i <= m ; i++) {
            a_sum_b.push_back(a[i]);
        }
    } else if (n > m) {
        for (size_t i = m; i <= n ; i++) {
            a_sum_b.push_back(b[i]);
        }
    }
    return a_sum_b;
}

std::vector<fp2_elem> Poly_scale(fp2_elem const &a, std::vector<fp2_elem> const &b) {
    std::vector<fp2_elem> ab;

    for (fp2_elem const &b_coeff : b) {
        ab.push_back(Fp2_mul(a, b_coeff));
    }

    return ab;
}

std::vector<fp2_elem> Xmin(fp2_elem const &a) {
    return {Fp2_negative(a), Fp2_one()};
}

std::vector<fp2_elem> F0(fp2_elem const &a) {

    std::vector<fp2_elem> f0 = Xmin(a);

    return Poly_mult(f0, f0);
}

std::vector<fp2_elem> F1(fp2_elem const &a, fp2_elem const &A) {

    std::vector<fp2_elem> t0{Fp2_zero(), a};
    std::vector<fp2_elem> t1{Fp2_one(), a};
    std::vector<fp2_elem> t2{a, Fp2_one()};
    t1 = Poly_mult(t1, t2);
    t0 = Poly_scale(Fp2_add(A, A), t0);

    std::vector<fp2_elem> f1 = Poly_add(t0, t1);

    fp2_elem min_two{NTL::ZZ_p(-2), NTL::ZZ_p(0)};
    return Poly_scale(min_two, f1);
}

std::vector<fp2_elem> F2(fp2_elem const &a) {
    std::vector<fp2_elem> f2{Fp2_negative(Fp2_one()), a};
    return Poly_mult(f2, f2);
}

// Product tree (from SCALLOP code)
std::vector<fp2_elem> product_tree(std::vector<std::vector<fp2_elem>> const &leaves) {
            if (leaves.empty())
                throw std::logic_error("no leaves");
            auto prev = leaves; //<- copies leaves to not modify the original list...
            while (prev.size() > 1) {
                std::vector<std::vector<fp2_elem>> next;
                {
                    for (size_t i = 0; i < prev.size()-1; i += 2)
                        next.push_back(Poly_mult(prev[i], prev[i+1]));
                    if (prev.size() % 2)
                        next.push_back(prev.back());
                }
                prev = next;
            }
            return prev[0];
        };

NTL::ZZ_pE ntl(fp2_elem const &a) {
    NTL::ZZ_pX a_poly;
    NTL::ZZ_pE a_ntl;


    NTL::SetCoeff(a_poly, 1);
    a_poly[0] = a.first;
    if (!IsZero(a.second)) {
        a_poly[1] = a.second;
    }

    NTL::conv(a_ntl, a_poly);

    return a_ntl;
}

fp2_elem fp2_conv(NTL::ZZ_pE const &a) {
    return {NTL::rep(a)[0], NTL::rep(a)[1]};
}

NTL::ZZ_pE Resultant(std::vector<fp2_elem> const &a, std::vector<fp2_elem> const &b) {
    NTL::vec_ZZ_pE a_vec, b_vec;
    NTL::ZZ_pEX a_ntl, b_ntl;
    for (fp2_elem const &a_coeff : a) {
        a_vec.append(ntl(a_coeff));
    }
    for (fp2_elem const &b_coeff : b) {
        b_vec.append(ntl(b_coeff));
    }
    NTL::conv(a_ntl, a_vec);
    NTL::conv(b_ntl, b_vec);

    return NTL::resultant(a_ntl, b_ntl);
}

ProjA SqrtVELU(xPoint &Q, ProjA &A, int ell, std::vector<xPoint> &evalPts) {
    
    NormalizeCoeff(A); 
    NormalizePoint(Q);
    NTL::ZZ_pE A_ntl = ntl(A.first);
    fp2_elem Del_IJ;

    std::vector<std::vector<fp2_elem>> h_I_list, D_J_list;   
    std::vector<std::vector<std::vector<fp2_elem>>> E_Js_list;

    std::vector<NTL::ZZ_pE> alphas;
    std::vector<fp2_elem> eval_xi;

    eval_xi.push_back(Fp2_one());
    eval_xi.push_back(Fp2_negative(Fp2_one()));


    for (xPoint P : evalPts) {
        NormalizePoint(P);
        eval_xi.push_back(P.first);        
        eval_xi.push_back(Fp2_inv(P.first));
    }

    for (fp2_elem const &xi : eval_xi) {
        alphas.push_back(ntl(xi));
    }


    std::vector<int> I, J, K;

    int m = NTL::SqrRoot(ell - 1)/2;
    int mm = (ell + 1)/(4*m);

    size_t I_start = 2*m;
    size_t I_step = 4*m;
    size_t I_end = 2*m*(2*mm-1)+1;

    size_t J_start = 1;
    size_t J_step = 2;
    size_t J_end = 2*m;

    size_t K_start = 4*m*mm+1;
    size_t K_step = 2;
    size_t K_end = ell;

    // loop through I
    xPoint R = xMUL(Q, NTL::ZZ(I_start), A);
    xPoint diff = xMUL(Q, NTL::ZZ((I_step - I_start) % ell), A);
    NormalizePoint(R);
    h_I_list.push_back(Xmin(R.first));

    xPoint sQ = xMUL(Q, NTL::ZZ(I_step), A);

    for (size_t i = I_start + I_step ; i < I_end ; i += I_step) {
        xPoint Rold = R;
        R = xADD(R, sQ, diff);
        NormalizePoint(R);

        h_I_list.push_back(Xmin(R.first));
        diff = Rold;
    }

    std::vector<fp2_elem> h_I = product_tree(h_I_list);

    size_t num_pts = alphas.size();

    // loop through J
    R = Q;
    fp2_elem xq = R.first;

    std::vector<fp2_elem> F0q = F0(xq);
    std::vector<fp2_elem> F1q = F1(xq, A.first);
    std::vector<fp2_elem> F2q = F2(xq);

    for (fp2_elem const &alpha : eval_xi) {
        std::vector<fp2_elem> f = Poly_scale(Fp2_sqr(alpha), F0q);
        f = Poly_add(f, Poly_scale(alpha, F1q));
        f = Poly_add(f, F2q);
        E_Js_list.push_back({f});
    }

    D_J_list.push_back(F0q);

    sQ = xDBL(Q, A);
    diff = Q;

    for (size_t j = J_start + J_step; j < J_end; j += J_step) {
        xPoint Rold = R;
        R = xADD(R, sQ, diff);
        NormalizePoint(R);

        assert (PointEqual(R, xMUL(Q, NTL::ZZ(j), A)));

        NTL::ZZ_pE xq = ntl(R.first);

        for (size_t i = 0 ; i < num_pts ; i++) {        
            std::vector<fp2_elem> f = Poly_scale(Fp2_sqr(eval_xi[i]), F0q);
            f = Poly_add(f, Poly_scale(eval_xi[i], F1q));
            f = Poly_add(f, F2q);
            E_Js_list[i].push_back(f);
        }
        D_J_list.push_back(F0q);
        diff = Rold;
    }

    std::vector<fp2_elem> D_J = product_tree(D_J_list);
    std::vector<std::vector<fp2_elem>> E_Js;

    for (std::vector<std::vector<fp2_elem>> E_J_list : E_Js_list) {
        E_Js.push_back(product_tree(E_J_list));
    }

    Del_IJ = fp2_conv(Resultant(h_I, D_J));

    std::vector<fp2_elem> h_S;

    //looping through K
    R = xMUL(Q, NTL::ZZ(K_start), A);
    NormalizePoint(R);
    for (fp2_elem alpha : eval_xi) {
        h_S.push_back(Fp2_sub(alpha, R.first));
    }

    diff = xMUL(Q, NTL::ZZ((K_start - K_step) % ell), A);
    for (size_t k = K_start + K_step; k < K_end ; k += K_step) {
        xPoint Rold = R;
        R = xADD(R, sQ, diff);
        NormalizePoint(R);

        assert (PointEqual(R, xMUL(Q, NTL::ZZ(k), A)));

        for (size_t i = 0 ; i < num_pts ; i++) {
            h_S[i] = Fp2_mul(h_S[i], Fp2_sub(eval_xi[i], R.first));
        }
        diff = Rold;
    }

    for (size_t i = 0 ; i < num_pts ; i++) {
        fp2_elem Ri = fp2_conv(Resultant(h_I, E_Js[i]));
        h_S[i] = Fp2_div(Fp2_mul(h_S[i], Ri), Del_IJ);
    }

    //Codomain curve
    fp2_elem two{NTL::ZZ_p(2), NTL::ZZ_p(0)};

    fp2_elem t0 = Fp2_sub(A.first, two);
    fp2_elem t1 = Fp2_add(A.first, two);

    fp2_elem d = Fp2_div(t0, t1);
    d = Fp2_pow(d, ell);

    t0 = Fp2_div(h_S[0], h_S[1]);
    t0 = Fp2_pow(t0, 8);
    d = Fp2_mul(d, t0);

    t0 = Fp2_add(d, Fp2_one());
    t0 = Fp2_add(t0, t0);

    t1 = Fp2_sub(Fp2_one(), d);

    ProjA A_new{t0, t1};

    //Points
    for (size_t i = 2; i < num_pts; i += 2) {
        fp2_elem X = Fp2_pow(eval_xi[i], ell);
        X = Fp2_mul(X, Fp2_sqr(h_S[i+1]));
        fp2_elem Z = Fp2_sqr(h_S[i]);

        evalPts[(i/2) - 1].first = X;
        evalPts[(i/2) - 1].second = Z;
    }

    NormalizeCoeff(A_new);

    return A_new;
}