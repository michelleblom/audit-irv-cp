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

enum Assertion { VIABLE, NONVIABLE, IRV, NEB, QSMAJ };

struct AuditSpec{
    Assertion type;
	double asn;
	
    int winner;
	int loser;

	Ints eliminated;

    double thresh;
};

typedef std::vector<AuditSpec> Audits;
typedef std::vector<Audits> Audits2d;

bool RevCompareAudit(const AuditSpec &a1, const AuditSpec &a2);

double EstimateASN_smajority(double tally, double other, double threshold_fr, 
    const Parameters &params);

double EstimateASN_cdiff(double tallyA, double tallyB, double d, 
    const Parameters &params);

double EstimateASN_VIABLE(const Contest &ctest, int c, const Ints &tallies, 
    int exhausted, const Parameters &params);

int estimate_sample_size(double margin, const Parameters &params);

double EstimateASN_NONVIABLE(const Contest &ctest, int c, const Ints &tallies,
    int exhausted, const Parameters &params); 

// Compute ASN to show that tail[0] beats one of tail[1..n] or i in winners
double FindBestIRV_NEB(const Contest &ctest, const Ints &tail, 
    const SInts &winners, const Parameters &params, const Ints &tallies, 
    const Audits2d &nebs, const Bools2d &has_neb, AuditSpec &audit);

#endif
