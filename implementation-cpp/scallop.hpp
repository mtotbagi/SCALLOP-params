#pragma once

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include "ec.hpp"
#include "ecp.hpp"
#include "isog.hpp"


std::pair<ecp, ecp> ActionIdealStep(std::pair<ecp, ecp> PQ, ecp &Qm, NTL::ZZ Lpos, NTL::ZZ Lneg, std::vector<int> const &ells, std::vector<int> &es);
std::pair<ecp, ecp> GroupAction(ecp const &P, ecp const &Q, std::vector<int> const &ells, std::vector<int> const &es);
