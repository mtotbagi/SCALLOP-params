#include <chrono>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include "ec.hpp"
#include "ecp.hpp"
#include "isog.hpp"
#include "scallop.hpp"

std::pair<ecp, ecp> ActionIdealStep(std::pair<ecp, ecp> PQ, NTL::ZZ const &Lpos, NTL::ZZ const &Lneg, std::vector<int> const &ells, std::vector<int> es)
{
    std::cout << "wtf" << std::endl;
    ecp P = PQ.first;
    ecp Q = PQ.second;
    std::cout << "PQ unwrapped" << std::endl;
    NTL::ZZ lampos, lamneg;
    NTL::conv(lampos, "1708298732195233441796205849078049269809774359155278336331280856185673548093108801649110316273383324791940297639831593458851060986908148784492132745287237220884212999179490292014463937711260");
    NTL::conv(lamneg, "21003476950938357369896031527264380182891980657183195527869820393272548129119708424839829416852343941922026863608741498103124098909400216982866894294803151154158272595988210522793797887274525");

    NTL::ZZ L = Lpos*Lneg;
    std::cout << "K" << std::endl;
    ecp K = P.curve().random_point_of_order_dividing(L);

    //Computing omega
    std::cout << "computing omega" << std::endl;
    std::vector<int> strat{256,128,64,32,17,9,5,3,2,1,1,1,1,2,1,1,1,4,2,1,1,1,2,1,1,8,4,2,1,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,64,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,128,64,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,64,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,32,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,16,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1,8,4,2,1,1,2,1,1,4,2,1,1,2,1,1};
    isog_2_chain phi_1{Q, strat};


    std::cout << "phi_1" << std::endl;
    isog_2_chain phi_2{P, strat};

    std::cout << "phi_2" << std::endl;
    ecp Pm = P.curve().completeBasis(2, 518, P);
    isog_2_chain phi_2_dual{phi_2(Pm), strat};


    std::cout << "phi_2_dual" << std::endl;

    std::cout << "computing isomorphism away" << std::endl;
    isomorphism isom_1(std::make_shared<const ec>(phi_1.get_codomain()), std::make_shared<const ec>(phi_2.get_codomain()));
    std::cout << "computing isomorphism at E" << std::endl;
    isomorphism isom_2(std::make_shared<const ec>(phi_2_dual.get_codomain()), std::make_shared<const ec>(P.curve()));

    std::cout << "done" << std::endl;
    //Kernel generator
    ecp wK = isom_2(phi_2_dual(isom_1(phi_1(K))));
    ecp K_pos = Lneg*(wK - (lampos % L)*K);
    ecp K_neg = Lpos*(wK - (lamneg % L)*K);

    ecp K_gen = K_pos + K_neg;

    ecp newP = P;
    ecp newQ = Q;

    //Compute isogeny, and update exponent vector
    for (size_t i = 0 ; i < ells.size() ; i++) {

        if (!(L % ells[i] == 0)) {
            continue;
        }

        std::cout << "computing isogeny of degree ell = " << ells[i] << std::endl;

        L /= ells[i];
        ecp Ki = L*K_gen;

        if (Ki) {
            assert (!(ells[i]*Ki));
            isog phi_i = isogeny(Ki, ells[i]);

            if (L != 1) {
                K_gen = phi_i(K_gen);
            }
            
            newP = phi_i(newP);
            newQ = phi_i(newQ);
            if (Lpos % ells[i] == 0) {
                es[i] -= 1;
            } else {
                es[i] += 1;
            }
        } else {
            std::cout << "Maybe funky stuff happens now?" << std::endl;
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

void printVector(std::vector<int> const &vec) {
    std::cout << "es: " << std::endl;
    for (int elem : vec) {
        std::cout << elem << ", ";
    }
    std::cout << std::endl;
}

std::pair<ecp, ecp> GroupAction(ecp const &P, ecp const &Q, std::vector<int> const &ells, std::vector<int> const &es_in)
{
    auto start = std::chrono::steady_clock::now();
    std::vector<int> es(es_in);
    std::pair<ecp, ecp> PiQi{P, Q};

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
        PiQi = ActionIdealStep(PiQi, Lpos, Lneg, ells, es);
        std::cout << "finished a step" << std::endl;

        auto duration = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start);
        std::cout << "We have used " << duration.count() << "seconds" << std::endl;
    }
    return PiQi;
}