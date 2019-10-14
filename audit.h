/*
    Copyright (C) 2018-2019  Michelle Blom

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef _AUDIT_H
#define _AUDIT_H

#include "model.h"

struct AuditSpec{
	double asn;
	int winner;
	int loser;

	Ints eliminated;
	bool wonly;

    Ints rules_out;
};

typedef std::vector<AuditSpec> Audits;


// The following two functions deal with generating audits
double EstimateASN_WONLY(const Ballots &rep_ballots, 
	const Candidates &cand, double rlimit, int winner, 
	int loser, double gamma, double lambda, int tot_auditable_ballots);

double EstimateSampleSize(const Ballots &rep_ballots, const Candidates &cand, 
	double rlimit, const Ints &tail, AuditSpec &best_audit,
	double gamma, double lambda, int tot_auditable_ballots);


#endif
