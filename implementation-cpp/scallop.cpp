#include <chrono>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include "ec.hpp"
#include "ecp.hpp"
#include "isog.hpp"
#include "scallop.hpp"

void printVector(std::vector<int> const &vec) {
    std::cout << "es: " << std::endl;
    for (int elem : vec) {
        std::cout << elem << ", ";
    }
    std::cout << std::endl;
}

std::pair<ecp, ecp> ActionIdealStep(std::pair<ecp, ecp> PQ, ecp &Qm, NTL::ZZ Lpos, NTL::ZZ Lneg, std::vector<int> const &ells, std::vector<int> &es)
{
    ecp P = PQ.first;
    ecp Q = PQ.second;
    NTL::ZZ lampos, lamneg;
    NTL::conv(lampos, "1708298732195233441796205849078049269809774359155278336331280856185673548093108801649110316273383324791940297639831593458851060986908148784492132745287237220884212999179490292014463937711260");
    NTL::conv(lamneg, "21003476950938357369896031527264380182891980657183195527869820393272548129119708424839829416852343941922026863608741498103124098909400216982866894294803151154158272595988210522793797887274525");

    NTL::ZZ L = Lpos*Lneg;

    std::cout << "K" << std::endl;
    ecp K = P.curve().random_point_of_order_dividing(L);

    while (!(K)) {
        K = P.curve().random_point_of_order_dividing(L);
    }

    //Computing omega
    std::vector<int> strat{256,128,64,32,17,9,5,3,2,1,1,1,1,2,1,1,1,4,2,1,1,1,2,1,1,8,4,2,1,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,64,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,128,64,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,64,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1};
    std::cout << "phi_1..." << std::endl;
    isog_2_chain phi_1{P, strat};

    std::cout << "phi_2..." << std::endl;
    isog_2_chain phi_2{Q, strat};

    std::cout << "phi_2_dual..." << std::endl;
    isog_2_chain phi_2_dual{phi_2(Qm), strat};

    std::cout << "computing isomorphisms" << std::endl;
    isomorphism isom_1(std::make_shared<const ec>(phi_1.get_codomain()), std::make_shared<const ec>(phi_2.get_codomain()));
    isomorphism isom_2(std::make_shared<const ec>(phi_2_dual.get_codomain()), std::make_shared<const ec>(P.curve()));

    std::cout << "done" << std::endl;
    //Kernel generator
    ecp wK = isom_2(phi_2_dual(isom_1(phi_1(K))));
    ecp K_pos = Lneg*(wK - (lamneg % Lpos)*K);
    ecp K_neg = K_pos;

    if (K_pos) {
        if (!(isom_2(phi_2_dual(isom_1(phi_1(K_pos)))) == lampos*K_pos)) {
            //Oops, we swapped some signs
            std::cout << "sign swap thing" << std::endl;
            wK = -wK;
            K_pos += (2*Lneg)*wK;
        }
        K_neg = Lpos*(wK - (lampos % Lneg)*K);
    } else {
        K_neg = Lpos*(wK - (lampos % Lneg)*K);
        if (!(isom_2(phi_2_dual(isom_1(phi_1(K_neg)))) == lamneg*K_neg)) {
            //Oops, we swapped some signs
            std::cout << "sign swap thing" << std::endl;
            wK = -wK;
            K_neg += (2*Lpos)*wK;
        }
    }


     

    ecp newP = P;
    ecp newQ = Q;

    //Compute isogeny, and update exponent vector
    for (size_t i = 0 ; i < ells.size() ; i++) {
        ecp Ki = K_pos;
        bool pos = true;
        if (Lpos % ells[i] == 0) {
            Lpos /= ells[i];
            Ki = Lpos*K_pos;

        } else if (Lneg % ells[i] == 0) {
            pos = false;
            Lneg /= ells[i];
            Ki = Lneg*K_neg;
        } else {
            continue;
        }

        std::cout << "at ell = " << ells[i] << std::endl;
        if (Ki) {
            std::cout << "computing isogeny of degree ell = " << ells[i] << std::endl;
            assert (!(ells[i]*Ki));
            isog phi_i = isogeny(Ki, ells[i]);

            if (Lpos != 1) {
                K_pos = phi_i(K_pos);
            }
            if (Lneg != 1) {
                K_neg = phi_i(K_neg);
            }
            
            newP = phi_i(newP);
            newQ = phi_i(newQ);
            Qm = phi_i(Qm);
            if (pos) {
                es[i] -= 1;
            } else {
                es[i] += 1;
            }
        } else {
            continue;
        }
    }

    std::pair<ecp, ecp> out{newP, newQ};
    return out;
}

bool _all_zero(std::vector<int> const &vec) {
    for (int elem : vec) {
        if (elem != 0) {
            return false;
        }
    }
    return true;
}

std::pair<ecp, ecp> GroupAction(ecp const &P, ecp const &Q, std::vector<int> const &ells, std::vector<int> const &es_in)
{
    auto start = std::chrono::steady_clock::now();
    std::vector<int> es(es_in);
    std::pair<ecp, ecp> PiQi{P, Q};
    ecp Qm = Q.curve().completeBasis(2, 518, Q);

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
        PiQi = ActionIdealStep(PiQi, Qm, Lpos, Lneg, ells, es);
        std::cout << "finished a step" << std::endl;

        auto duration = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start);
        std::cout << "We have used " << duration.count() << "seconds" << std::endl;
    }
    return PiQi;
}