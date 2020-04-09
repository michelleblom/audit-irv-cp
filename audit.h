/*
    Copyright (C) 2018-2020  Michelle Blom

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

enum Assertion { VIABLE, NONVIABLE, IRV };

struct AuditSpec{
    Assertion type;
	double asn;
	
    int winner;
	int loser;

	Ints eliminated;
};

typedef std::vector<AuditSpec> Audits;

struct Parameters{
    double risk_limit;
    int tot_auditable_ballots;
};

bool RevCompareAudit(const AuditSpec &a1, const AuditSpec &a2);

double EstimateASN_VIABLE(const Contest &ctest, int c, const Ints &elim,
    const Parameters &params);

double EstimateASN_NONVIABLE(const Contest &ctest, int c, const Ints &elim,
    const Parameters &params); 

// Compute ASN to show that tail[0] beats one of tail[1..n] or i in winners
double FindBestIRV(const Contest &ctest, const Ints &tail, 
    const SInts &winners, const Parameters &params, AuditSpec &audit);

#endif
