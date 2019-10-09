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

struct Result{
	int polls;
	int remaining_hypotheses;
};

typedef std::vector<AuditSpec> Audits;

double EstimateASN_WONLY(const Ballots &rep_ballots, 
	const Candidates &cand, double rlimit, int winner, 
	int loser, double gamma, double lambda);

double EstimateSampleSize(const Ballots &rep_ballots, const Candidates &cand, 
	double rlimit, const Ints &tail, AuditSpec &best_audit,
	double gamma, double lambda);

Result RunSingleWinnerLoserAudit(const Ballots &act_ballots, 
	const Ballots &rep_ballots, const Candidates &cand,
	double rlimit, int winner, int loser, const Ints &plist,
	double gamma, double lambda);

Result RunAudit(const Ballots &act_ballots, const Ballots &rep_ballots,
	const Candidates &cand, double rlimit, const Ints &winners, 
	const Ints &losers,const Ints &plist, double gamma, double lambda);  

bool LoadAudits(const char *path, Audits &audits,
	const Candidates &candidates, const Config &config);

#endif
