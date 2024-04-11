#pragma once

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <vector>

#include "fp2.hpp"
#include "montgomery.hpp"

// Product tree (from SCALLOP code)
template <class I>
I product_tree(std::vector<I> const &leaves) {
            if (leaves.empty())
                throw std::logic_error("no leaves");
            auto prev = leaves; //<- copies leaves to not modify the original list...
            while (prev.size() > 1) {
                std::vector<I> next;
                {
                    for (size_t i = 0; i < prev.size()-1; i += 2)
                        next.push_back(prev[i] * prev[i+1]);
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
    a_poly[1] = a.second;

    NTL::conv(a_ntl, a_poly);

    return a_ntl;
}

fp2_elem fp2_conv(NTL::ZZ_pE const &a) {
    return {NTL::rep(a)[0], NTL::rep(a)[1]};
}

NTL::ZZ_pEX F0(NTL::ZZ_pE const &a) {
    NTL::ZZ_pEX X, f0;
    NTL::SetX(X);

    f0 = X - a;

    return f0*f0;
}

NTL::ZZ_pEX F1(NTL::ZZ_pE const &a, NTL::ZZ_pE const &A) {
    NTL::ZZ_pEX X, f1;
    NTL::SetX(X);

    NTL::ZZ_pEX t0 = X*a;
    f1 = (t0 + 1)*(X + a) + 2*A*t0;

    return -2*f1;
}

NTL::ZZ_pEX F2(NTL::ZZ_pE const &a) {
    NTL::ZZ_pEX X, f2;
    NTL::SetX(X);

    f2 = X*a - 1;

    return f2*f2;
}

ProjA SqrtVELU(xPoint &Q, ProjA &A, int ell, std::vector<xPoint> &evalPts) {
    
    NormalizeCoeff(A); 
    NormalizePoint(Q);
    NTL::ZZ_pE A_ntl = ntl(A.first);
    fp2_elem Del_IJ;

    std::vector<NTL::ZZ_pEX> h_I_list, D_J_list;   
    std::vector<std::vector<NTL::ZZ_pEX>> E_Js_list;
    NTL::ZZ_pEX X;
    NTL::SetX(X);

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
    h_I_list.push_back(X - ntl(R.first));

    xPoint sQ = xMUL(Q, NTL::ZZ(I_step), A);

    for (size_t i = I_start + I_step ; i < I_end ; i += I_step) {
        xPoint Rold = R;
        R = xADD(R, sQ, diff);
        NormalizePoint(R);

        h_I_list.push_back(X - ntl(R.first));
        diff = Rold;
    }

    NTL::ZZ_pEX h_I = product_tree<NTL::ZZ_pEX>(h_I_list);

    size_t num_pts = alphas.size();

    // loop through J
    R = Q;
    NTL::ZZ_pE xq = ntl(R.first);

    for (NTL::ZZ_pE const &alpha : alphas) {
        E_Js_list.push_back({F0(xq)*alpha*alpha + F1(xq, A_ntl)*alpha + F2(xq)});
    }
    D_J_list.push_back(F0(xq));

    sQ = xDBL(Q, A);
    diff = Q;

    for (size_t j = J_start + J_step; j < J_end; j += J_step) {
        xPoint Rold = R;
        R = xADD(R, sQ, diff);
        NormalizePoint(R);

        assert (PointEqual(R, xMUL(Q, NTL::ZZ(j), A)));

        NTL::ZZ_pE xq = ntl(R.first);

        NTL::ZZ_pEX F0q = F0(xq);
        NTL::ZZ_pEX F1q = F1(xq, A_ntl);
        NTL::ZZ_pEX F2q = F2(xq);

        for (size_t i = 0 ; i < num_pts ; i++) {
            NTL::ZZ_pE alpha = alphas[i];
            E_Js_list[i].push_back(F0q*alpha*alpha + F1q*alpha + F2q);
        }
        D_J_list.push_back(F0q);
        diff = Rold;
    }

    NTL::ZZ_pEX D_J = product_tree<NTL::ZZ_pEX>(D_J_list);
    std::vector<NTL::ZZ_pEX> E_Js;
    for (std::vector<NTL::ZZ_pEX> E_J_list : E_Js_list) {
        E_Js.push_back(product_tree<NTL::ZZ_pEX>(E_J_list));
    }

    Del_IJ = fp2_conv(NTL::resultant(h_I, D_J));

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
        fp2_elem Ri = fp2_conv(NTL::resultant(h_I, E_Js[i]));
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