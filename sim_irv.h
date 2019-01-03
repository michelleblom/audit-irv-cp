#ifndef _SIM_IRV_H
#define _SIM_IRV_H

#include "model.h"

void SimIRV(const Ballots &ballots, int &winner, Candidates &cands,
	const Config &config, Ints &order_c, bool log);

void SimulateElimination(int cidx, const Ballots &ballots, Candidates &cands);

#endif
