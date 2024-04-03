#pragma once

#include <cassert>
#include <vector>
#include <unordered_map>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pEX.h>

#include "ec.hpp"
#include "ecp.hpp"


typedef NTL::ZZ Integer; //NTL::ZZ is probably overkill, optimise later

std::unordered_map<int, int> factor(Integer const &N);
NTL::ZZ DLP(ecp const &Q, ecp const &base, int ell, int e);
std::unordered_map<int, int> factor(Integer const &N);
std::optional<NTL::ZZ_pE> sqrt(NTL::ZZ_pE const &alpha);

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


