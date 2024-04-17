#pragma once

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <vector>

#include "fp2.hpp"
#include "montgomery.hpp"
#include "poly.hpp"

/*
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
    t1 = Poly_add(t0, t1);
    fp2_elem min_two{NTL::ZZ_p(-2), NTL::ZZ_p(0)};
    return Poly_scale(min_two, t1);
}

std::vector<fp2_elem> F2(fp2_elem const &a) {
    std::vector<fp2_elem> f2{Fp2_negative(Fp2_one()), a};
    return Poly_mult(f2, f2);
}

ProjA SqrtVELU(xPoint &Q, ProjA &A, int ell, std::vector<xPoint> &evalPts) {
    
    NormalizeCoeff(A); 
    NormalizePoint(Q);
    fp2_elem Del_IJ;

    std::vector<fp2X_elem> h_I_list, D_J_list;

    std::vector<std::vector<fp2X_elem>> E_Js_list;

    std::vector<fp2_elem> eval_xi;

    eval_xi.push_back(Fp2_one());
    eval_xi.push_back(Fp2_negative(Fp2_one()));


    for (xPoint P : evalPts) {
        NormalizePoint(P);
        eval_xi.push_back(P.first);        
        eval_xi.push_back(Fp2_inv(P.first));
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


    auto start = std::chrono::steady_clock::now();
    for (size_t i = I_start + I_step ; i < I_end ; i += I_step) {
        xPoint Rold = R;
        R = xADD(R, sQ, diff);
        NormalizePoint(R);

        h_I_list.push_back(Xmin(R.first));
        diff = Rold;
    }
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "Iterating over I took: " << duration.count() << " microseconds" << std::endl;

    start = std::chrono::steady_clock::now();
    auto h_I_tree = product_tree(h_I_list);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "Product tree took: " << duration.count() << " microseconds" << std::endl;

    std::vector<fp2_elem> h_I = h_I_tree[0][0];

    size_t num_pts = eval_xi.size();

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

    start = std::chrono::steady_clock::now();
    for (size_t j = J_start + J_step; j < J_end; j += J_step) {
        xPoint Rold = R;
        R = xADD(R, sQ, diff);
        NormalizePoint(R);

        assert (PointEqual(R, xMUL(Q, NTL::ZZ(j), A)));

        fp2_elem xq = R.first;
        F0q = F0(xq);
        F1q = F1(xq, A.first);
        F2q = F2(xq);

        for (size_t i = 0 ; i < num_pts ; i++) {        
            std::vector<fp2_elem> f = Poly_scale(Fp2_sqr(eval_xi[i]), F0q);
            f = Poly_add(f, Poly_scale(eval_xi[i], F1q));
            f = Poly_add(f, F2q);
            E_Js_list[i].push_back(f);
        }
        D_J_list.push_back(F0q);
        diff = Rold;
    }

    duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "Iterating over J took: " << duration.count() << " microseconds" << std::endl;

    start = std::chrono::steady_clock::now();
    fp2X_elem D_J = product_tree(D_J_list)[0][0];
    duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "D_J product tree: " << duration.count() << " microseconds" << std::endl;

    start = std::chrono::steady_clock::now();
    fp2X_elem D_J_lin = D_J_list[0];
    for (size_t i = 1 ; i < D_J_list.size() ; i++) {
        D_J_lin = Poly_mult(D_J_lin, D_J_list[i]);
    }
    duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "D_J linear tree: " << duration.count() << " microseconds" << std::endl;

    std::vector<fp2X_elem> E_Js;

    start = std::chrono::steady_clock::now();
    for (std::vector<fp2X_elem> E_J_list : E_Js_list) {
        start = std::chrono::steady_clock::now();
        E_Js.push_back(product_tree(E_J_list)[0][0]);
        duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
        std::cout << "E_J product tree: " << duration.count() << " microseconds" << std::endl;
    }
    //Del_IJ = fp2_conv(Resultant(h_I, D_J));
    start = std::chrono::steady_clock::now();
    auto h_I_tree_rev_inv = rev_inv_tree(h_I_tree, D_J);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "Computing rev_inv_tree took: " << duration.count() << " microseconds" << std::endl;

    start = std::chrono::steady_clock::now();
    Del_IJ = Fast_Resultant(D_J, h_I_tree, h_I_tree_rev_inv);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "Resultant took: " << duration.count() << " microseconds" << std::endl;

    std::vector<fp2_elem> h_S;

    //looping through K
    R = xMUL(Q, NTL::ZZ(K_start), A);
    NormalizePoint(R);
    for (fp2_elem alpha : eval_xi) {
        h_S.push_back(Fp2_sub(alpha, R.first));
    }

    start = std::chrono::steady_clock::now();
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
    duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "Iterating over K took: " << duration.count() << " microseconds" << std::endl;

    for (size_t i = 0 ; i < num_pts ; i++) {
        start = std::chrono::steady_clock::now();
        fp2_elem Ri = Fast_Resultant(E_Js[i], h_I_tree, h_I_tree_rev_inv);
        h_S[i] = Fp2_div(Fp2_mul(h_S[i], Ri), Del_IJ);
        duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
        std::cout << "E_J resultant took: " << duration.count() << " microseconds" << std::endl;
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
*/