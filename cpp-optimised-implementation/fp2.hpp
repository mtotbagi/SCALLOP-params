#pragma once

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

typedef std::pair<NTL::ZZ_p, NTL::ZZ_p> fp2_elem;

fp2_elem Fp2_add(fp2_elem const &a, fp2_elem const &b) {
    return {a.first + b.first, a.second + b.second};
}

fp2_elem Fp2_sub(fp2_elem const &a, fp2_elem const &b) {
    return {a.first - b.first, a.second - b.second};
}

fp2_elem Fp2_mul(fp2_elem const &a, fp2_elem const &b) {
    NTL::ZZ_p t0 = a.first * b.first;
    NTL::ZZ_p t1 = a.second * b.second;
    NTL::ZZ_p t2 = a.first + a.second;
    NTL::ZZ_p t3 = b.first + b.second;

    NTL::ZZ_p res0 = t0 - t1;
    t0 = t0 + t1;
    t1 = t2 * t3;

    return {res0, t1 - t0};
}

fp2_elem Fp2_sqr(fp2_elem const &a) {
    NTL::ZZ_p t0 = a.first + a.second;
    NTL::ZZ_p t1 = a.first - a.second;
    NTL::ZZ_p t2 = a.first + a.first;

    return {t0*t1, t2*a.second};
}

fp2_elem Fp2_inv(fp2_elem const &a) {
    NTL::ZZ_p t0 = a.first * a.first;
    NTL::ZZ_p t1 = a.second * a.second;
    t0 = 1/(t0 + t1);

    return {a.first * t0, -(a.second * t0)};
}

void Fp2_neg(fp2_elem &a) {
    //negate in place
    a.first = -a.first;
    a.second = -a.second;
}

fp2_elem Fp2_div(fp2_elem const &a, fp2_elem const &b) {
    return Fp2_mul(a, Fp2_inv(b));
}

fp2_elem Fp2_one() {
    return {NTL::ZZ_p(1), NTL::ZZ_p(0)};
}

fp2_elem Fp2_zero() {
    return {NTL::ZZ_p(0), NTL::ZZ_p(0)};
}

bool IsZero(fp2_elem const &a) {
    return (a.first == 0) && (a.second == 0);
}

bool IsOne(fp2_elem const &a) {
    return (a.first == 1) && (a.second == 0);
}

bool Fp2_equal(fp2_elem const &a, fp2_elem const &b) {
    return IsZero(Fp2_sub(a, b));
}

//// Printing
void PrintFp2(fp2_elem const &a) {
    std::cout << a.first << " + i*" << a.second;
}