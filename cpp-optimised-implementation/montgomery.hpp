#pragma once

#include "fp2.hpp"

typedef std::pair<fp2_elem, fp2_elem> xPoint;
typedef std::pair<fp2_elem, fp2_elem> ProjA;

void NormalizePoint(xPoint &P) {
    //normalize in place
    if (IsZero(P.second)) {
        P.first = Fp2_one();
    } else {
        P.first = Fp2_div(P.first, P.second);
        P.second = Fp2_one();
    }
}

void NormalizeCoeff(ProjA &A) {
    //normalize in place
    assert (!(IsZero(A.second)));
    A.first = Fp2_div(A.first, A.second);
    A.second = Fp2_one();
}

bool PointEqual(xPoint P, xPoint Q) {
    NormalizePoint(P);
    NormalizePoint(Q);
    return Fp2_equal(P.first, Q.first) && Fp2_equal(P.second, Q.second);
}

bool IsIdentity(xPoint P) {
    return IsZero(P.second);
}

xPoint xDBL(xPoint const &P, ProjA const &A) {
    // Point with (X : Z), and on a curve with montgomery A
    fp2_elem t0 = Fp2_sub(P.first, P.second);
	fp2_elem t1 = Fp2_add(P.first, P.second);
	t0 = Fp2_sqr(t0);
	t1 = Fp2_sqr(t1);

	fp2_elem Z2 = Fp2_mul(A.second, t0);
	Z2 = Fp2_add(Z2, Z2);
	Z2 = Fp2_add(Z2, Z2);

	fp2_elem X2 = Fp2_mul(Z2, t1);
	t1 = Fp2_sub(t1, t0);
	t0 = Fp2_add(A.second, A.second);
	t0 = Fp2_add(t0, A.first);
	t0 = Fp2_mul(t0, t1);
	Z2 = Fp2_add(Z2, t0);
	Z2 = Fp2_mul(Z2, t1);

	return {X2, Z2};
}

std::pair<xPoint, xPoint> xDBLADD(xPoint const &P, xPoint const &Q, xPoint const &PmQ, ProjA const &A) {
    // Differential addition
    // Returns 2P, P+Q

    fp2_elem t0 = Fp2_add(P.first, P.second);
	fp2_elem t1 = Fp2_sub(P.first, P.second);
	fp2_elem X2P = Fp2_sqr(t0);
	fp2_elem t2 = Fp2_sub(Q.first, Q.second);
	fp2_elem XQP = Fp2_add(Q.first, Q.second);
	t0 = Fp2_mul(t0, t2);

	fp2_elem Z2P = Fp2_sqr(t1);
	t1 = Fp2_mul(t1, XQP);
	t2 = Fp2_sub(X2P, Z2P);

	fp2_elem t3 = Fp2_add(A.second, A.second);
	Z2P = Fp2_mul(t3, Z2P);
	t3 = Fp2_add(A.first, t3);
	Z2P = Fp2_add(Z2P, Z2P);
	X2P = Fp2_mul(X2P, Z2P);
	XQP = Fp2_mul(t3, t2);

	fp2_elem ZQP = Fp2_sub(t0, t1);
	Z2P = Fp2_add(XQP, Z2P);
	XQP = Fp2_add(t0, t1);
	Z2P = Fp2_mul(Z2P, t2);
	ZQP = Fp2_sqr(ZQP);
	XQP = Fp2_sqr(XQP);
	ZQP = Fp2_mul(PmQ.first, ZQP);
	XQP = Fp2_mul(XQP, PmQ.second);

	return {{X2P, Z2P}, {XQP, ZQP}};
}

xPoint IdentityPoint() {
    return {Fp2_one(), Fp2_zero()};
}

xPoint xMUL(xPoint const &P, NTL::ZZ const &m, ProjA A) {

    xPoint P0 = IdentityPoint();
    xPoint P1 = P;
    for (int i = NTL::NumBits(m); i >= 0; i--) {
        if (((m >> i) & 1) == 0) {
            auto output = xDBLADD(P0, P1, P, A);
            P0 = output.first;
            P1 = output.second;
        } else {
            auto output = xDBLADD(P1, P0, P, A);
            P1 = output.first;
            P0 = output.second;
        }
    }
    return P0;
}

//// Printing

void PrintPoint(xPoint const &P) {
    std::cout << "(";
    PrintFp2(P.first);
    std::cout << " : ";
    PrintFp2(P.second);
    std::cout << ")";
}