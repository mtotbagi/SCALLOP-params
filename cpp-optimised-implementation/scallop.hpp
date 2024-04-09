#pragma once

#include <random>

#include "fp2.hpp"
#include "montgomery.hpp"
#include "isog.hpp"


void printVector(std::vector<int> const &vec) {
    std::cout << "es: " << std::endl;
    for (int elem : vec) {
        std::cout << elem << ", ";
    }
    std::cout << std::endl;
}

ProjA ActionIdealStep(xPoint &P, xPoint &Q, xPoint &Qm, ProjA A, NTL::ZZ Lpos, NTL::ZZ Lneg, std::vector<int> const &ells, std::vector<int> &es)
{
    NormalizeCoeff(A);
    NTL::ZZ lampos, lamneg;
    NTL::conv(lampos, "1708298732195233441796205849078049269809774359155278336331280856185673548093108801649110316273383324791940297639831593458851060986908148784492132745287237220884212999179490292014463937711260");
    NTL::conv(lamneg, "21003476950938357369896031527264380182891980657183195527869820393272548129119708424839829416852343941922026863608741498103124098909400216982866894294803151154158272595988210522793797887274525");

    NTL::ZZ L = Lpos*Lneg;

    std::cout << "K" << std::endl;
    xPoint K = PointOfOrderDividing(A, L);

    while (IsIdentity(K)) {
        K = PointOfOrderDividing(A, L);
    }

    //Computing omega
    std::vector<int> strat{256, 128, 64, 32, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 128, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1};
    
    std::vector<xPoint> KerGens_1, KerGens_2, KerGens_3;
    ProjA Am = TwoIsogChainPrecompute(P, A, strat, KerGens_1);
    std::cout << "j-invariant away P: " << StringFp2(jInvariant(Am)) << std::endl;

    ProjA Amm = TwoIsogChainPrecompute(Q, A, strat, KerGens_2);
    std::cout << "j-invariant away Q: " << StringFp2(jInvariant(Amm)) << std::endl;

    auto ur = IsomorphismConstants(Amm, Am);
    xPoint Qmm = TwoIsogChainEvaluate(Qm, KerGens_2);

    Qmm = IsomorphismEval(Qmm, ur);
    ProjA A_start = TwoIsogChainPrecompute(Qmm, Am, strat, KerGens_3);

    std::cout << "j-invariant orig: " << StringFp2(jInvariant(A)) << std::endl;
    std::cout << "j-invariant recovered: " << StringFp2(jInvariant(A_start)) << std::endl;

    auto ur_home = IsomorphismConstants(A_start, A);

    std::cout << "Finding K and evaluating" << std::endl;

    xPoint wK = TwoIsogChainEvaluate(K, KerGens_1);
    wK = TwoIsogChainEvaluate(wK, KerGens_3);
    wK = IsomorphismEval(wK, ur_home);

    xPoint lamnegK = xMUL(K, lamneg % Lpos, A);
    xPoint lamposK = xMUL(K, lampos % Lneg, A);

    auto K_pos_both = xADDSUB(wK, lamnegK, A);
    xPoint Kp = xMUL(K_pos_both.first, Lneg, A);

    xPoint wKp = TwoIsogChainEvaluate(Kp, KerGens_1);
    wKp = TwoIsogChainEvaluate(wKp, KerGens_3);
    wKp = IsomorphismEval(wKp, ur_home);

    xPoint lamKp = xMUL(Kp, lampos, A);

    if (!(PointEqual(lamKp, wKp))) {
        std::cout << "SIGN FLIP THING ON Kp" << std::endl;
        Kp = xMUL(K_pos_both.second, Lneg, A);
    }

    auto K_neg_both = xADDSUB(wK, lamposK, A);
    xPoint Kn = xMUL(K_neg_both.first, Lpos, A);

    xPoint wKn = TwoIsogChainEvaluate(Kn, KerGens_1);
    wKn = TwoIsogChainEvaluate(wKn, KerGens_3);
    wKn = IsomorphismEval(wKn, ur_home);

    xPoint lamKn = xMUL(Kp, lamneg, A);

    if (!(PointEqual(lamKn, wKn))) {
        std::cout << "SIGN FLIP THING ON Kn" << std::endl;
        Kn = xMUL(K_neg_both.second, Lpos, A);
    }

    ProjA Ai = A;

    //Compute isogeny, and update exponent vector
    std::cout << "Computing odd isogenies..." << std::endl;
    for (size_t i = 0 ; i < ells.size() ; i++) {
        xPoint Ki;
        bool pos = true;
        if (Lpos % ells[i] == 0) {
            Lpos /= ells[i];
            Ki = xMUL(Kp, Lpos, Ai);

        } else if (Lneg % ells[i] == 0) {
            pos = false;
            Lneg /= ells[i];
            Ki = xMUL(Kn, Lneg, Ai);
        } else {
            continue;
        }

        if (!(IsIdentity(Ki))) {
            assert (IsIdentity(xMUL(Ki, NTL::ZZ(ells[i]), Ai)));

            std::vector<xPoint> evalPts{P, Q, Qm};

            if (Lpos != 1) {
                evalPts.push_back(Kp);
            }
            if (Lneg != 1) {
                evalPts.push_back(Kn);
            }
            
            NormalizePoint(Ki);
            NormalizeCoeff(Ai);
            Ai = xISOG(Ki, Ai, ells[i], evalPts);

            P = evalPts[0];
            Q = evalPts[1];
            Qm = evalPts[2];

            if (Lpos != 1) {
                Kp = evalPts[3];
                if (Lneg != 1) {
                    Kn = evalPts[4];
                }
            } else if (Lneg != 1) {
                Kn = evalPts[3];
            }

            if (pos) {
                es[i] -= 1;
            } else {
                es[i] += 1;
            }
        } else {
            continue;
        }
    }
    return Ai;
}

bool _all_zero(std::vector<int> const &vec) {
    for (int elem : vec) {
        if (elem != 0) {
            return false;
        }
    }
    return true;
}

ProjA GroupAction(xPoint &P, xPoint &Q, xPoint &Qm, ProjA const &A, std::vector<int> const &ells, std::vector<int> const &es_in)
{
    auto start = std::chrono::steady_clock::now();
    std::vector<int> es(es_in);
    
    ProjA Ai = A;

    while (!(_all_zero(es))) {
        printVector(es);
        NTL::ZZ Lpos(1);
        NTL::ZZ Lneg(1);

        for (size_t i = 0 ; i < ells.size() ; i++) {
            if (es[i] > 0) {
                Lpos *= ells[i];
            }
            if (es[i] < 0) {
                Lneg *= ells[i];
            }
        }
        std::cout << "Going into step" << std::endl;
        Ai = ActionIdealStep(P, Q, Qm, Ai, Lpos, Lneg, ells, es);
        std::cout << "finished a step" << std::endl;

        auto duration = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start);
        std::cout << "We have used " << duration.count() << "seconds" << std::endl;
    }
    return Ai;
}

std::vector<int> GenSecret(size_t length, int N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(-N, N);

    std::vector<int> vec(length);
    for (int i = 0; i < length; ++i) {
        vec[i] = dis(gen);
    }
    return vec;
}