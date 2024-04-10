#pragma once

#include <NTL/ZZ_p.h>
#include <vector>

#include "fp2.hpp"
#include "montgomery.hpp"

ProjA ProjAToProjAC(ProjA A) {
    fp2_elem C24 = Fp2_add(A.second, A.second);
    fp2_elem A24 = Fp2_add(A.first, C24);
    C24 = Fp2_add(C24, C24);

    return {A24, C24};
}

ProjA ProjACToProjA(ProjA A24) {
	fp2_elem A = Fp2_add(A24.first, A24.first);
	A = Fp2_sub(A, A24.second);
	A = Fp2_add(A, A);

	return {A, A24.second};
}


ProjA TwoIsoCurve(xPoint const &P, xPoint &K) {
    // Compute codomain curve A
    // And constant K for eval
	assert (!(IsZero(P.first)));

	fp2_elem A24 = Fp2_sqr(P.first);
	fp2_elem C24 = Fp2_sqr(P.second);
	A24 = Fp2_sub(C24, A24);

	fp2_elem K0 = Fp2_add(P.first, P.second);
	fp2_elem K1 = Fp2_sub(P.first, P.second);

    K = {K0, K1};

	return {A24, C24};
}

xPoint TwoIsogEval(xPoint const &Q, xPoint const &K) {

	fp2_elem t0 = Fp2_add(Q.first, Q.second);
	fp2_elem t1 = Fp2_sub(Q.first, Q.second);
	fp2_elem t2 = Fp2_mul(K.first, t1);
	t1 = Fp2_mul(K.second, t0);
	t0 = Fp2_add (t2, t1);
	t1 = Fp2_sub(t2, t1);

	fp2_elem QEvalX = Fp2_mul(Q.first, t0);
	fp2_elem QEvalZ = Fp2_mul(Q.second, t1);

	return {QEvalX, QEvalZ};
}

ProjA TwoIsogChainPrecompute(xPoint const &Q, ProjA &A_in, std::vector<int> const &strategy, std::vector<xPoint> &KerGens) {

    size_t MAX = strategy.size() + 1;
    std::vector<xPoint> pts;
    pts.reserve(MAX);
    std::vector<size_t> pts_index;
    pts_index.reserve(MAX);
	size_t npts = 0;
	size_t ii = 0;

    xPoint P = Q;
    size_t index = 0;

    ProjA A = A_in;
    ProjA A24;
    xPoint K;

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
            P = xMUL(P, NTL::power(NTL::ZZ(2), m), A);
            index += m;
        }

        assert (!IsIdentity(P));
        assert (IsIdentity(xDBL(P, A)));

        A24 = TwoIsoCurve(P, K);
        KerGens.push_back(K);
        for (size_t i = 0; i < npts; i++) {
            pts[i] = TwoIsogEval(pts[i], K);
        }

        P = pts[npts-1];
        index = pts_index[npts - 1];
        npts -= 1;
        A = ProjACToProjA(A24);
    }

    A24 = TwoIsoCurve(P, K);
    KerGens.push_back(K);
    A = ProjACToProjA(A24);

    return A;
}

xPoint TwoIsogChainEvaluate(xPoint P, std::vector<xPoint> const &KerGens) {
    for (xPoint const &K : KerGens) {
        P = TwoIsogEval(P, K);
    }
    return P;
}

std::pair<fp2_elem, fp2_elem> IsomorphismConstants(ProjA &A, ProjA &Am) {
    // returns u^2 and r
    assert (Fp2_equal(jInvariant(A), jInvariant(Am)));
    NormalizeCoeff(A);
    NormalizeCoeff(Am);

    fp2_elem A_aff = A.first;
    fp2_elem Am_aff = Am.first;

    if (Fp2_equal(A_aff, Am_aff)) {
        return {Fp2_one(), Fp2_zero()};
    } else if (IsZero(Fp2_add(A_aff, Am_aff))) {
        return {Fp2_i(), Fp2_zero()};
    } else {
        fp2_elem u, r;
        fp2_elem six{NTL::ZZ_p(6), NTL::ZZ_p(0)};
        fp2_elem nine{NTL::ZZ_p(9), NTL::ZZ_p(0)};
        fp2_elem A_sqr = Fp2_sqr(A_aff);
        fp2_elem Am_sqr = Fp2_sqr(Am_aff);

        fp2_elem t0 = Fp2_add(A_sqr, Am_sqr);
        r = Fp2_sub(t0, six);
        r = Fp2_mul(r, A_aff);
        t0 = Fp2_add(t0, Am_sqr);
        t0 = Fp2_sub(t0, nine);
        r = Fp2_div(r, t0);

        fp2_elem t1 = Fp2_add(r, r);
        t1 = Fp2_add(t1, r);
        t1 = Fp2_sub(A_aff, t1);
        u = Fp2_div(Am_aff, t1);

        return {u, r};
    }
}

xPoint IsomorphismEval(xPoint &P, std::pair<fp2_elem, fp2_elem> ur) {
    NormalizePoint(P);
    return {Fp2_mul(ur.first, Fp2_add(P.first, ur.second)), Fp2_one()};
}

ProjA xISOG(xPoint const &P, ProjA const &A, int ell, std::vector<xPoint> &evalPts) {
    assert (IsOne(P.second));
    assert (IsOne(A.second));

    std::vector<fp2_elem> eval_xi;
    eval_xi.push_back(Fp2_one());
    eval_xi.push_back(Fp2_negative(Fp2_one()));

    for (xPoint &P : evalPts) {
        NormalizePoint(P);
        eval_xi.push_back(P.first);
        eval_xi.push_back(Fp2_inv(P.first));
    }

    size_t num_pts = eval_xi.size();

    std::vector<fp2_elem> h_S;

    for (size_t i = 0; i < num_pts; i++) {
        h_S.push_back(Fp2_sub(eval_xi[i], P.first));
    }

    xPoint Pprev = P;
    xPoint Pi = xDBL(P, A);
    xPoint temp;
    
    for (size_t j = 2; j <= ell/2; j++) {
        NormalizePoint(Pi);
        for (size_t i = 0; i < num_pts; i++) {
            h_S[i] = Fp2_mul(h_S[i], Fp2_sub(eval_xi[i], Pi.first));
        }
        temp = Pi;
        Pi = xADD(Pi, P, Pprev);
        Pprev = temp;
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

ProjA xISOG_Chain(xPoint &K, ProjA const &A, NTL::ZZ &L, std::vector<xPoint> &evalPts, int adj, std::vector<int> const &ells, std::vector<int> &es) {
    // Sub-routine, also updates es
    ProjA Ai = A;
    for (size_t i = 0 ; i < ells.size() ; i++) {
        xPoint Ki;
        if (L % ells[i] == 0) {
            L /= ells[i];
            Ki = xMUL(K, L, Ai);
        } else {
            continue;
        }

        if (!(IsIdentity(Ki))) {
            assert (IsIdentity(xMUL(Ki, NTL::ZZ(ells[i]), Ai)));
            if (L != 1) {
                evalPts.push_back(K);
            }
            
            NormalizePoint(Ki);
            NormalizeCoeff(Ai);
            Ai = xISOG(Ki, Ai, ells[i], evalPts);

            if (L != 1) {
                K = evalPts[evalPts.size() - 1];
                evalPts.pop_back();
            }
            es[i] -= adj;
        } else {
            continue;
        }
    }
    return Ai;
}