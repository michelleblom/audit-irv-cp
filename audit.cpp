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

#include "audit.h"
#include<fstream>
#include<iostream>
#include<cmath>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <limits>

using namespace std;

typedef boost::char_separator<char> boostcharsep;

namespace pt = boost::property_tree;

double kaplan_kolmogorov(Doubles &x, int xsize, int N, double t, double g){
    double sample_total = 0;

    double mart = (t > 0) ? (x[0]+g)/(t+g) : 1;
    double mart_max = mart;
    for(int j = 1; j < xsize; ++j){
        mart *= (x[j]+g)*(1-j/N)/(t+g - (1/N)*sample_total);
        if(mart < 0){
            break;
        }
        else {
            sample_total += x[j]+g;
        }
        mart_max = max(mart, mart_max);
    }
    return min(1.0/mart_max,1.0);
}

int estimate_sample_size(double margin, int max_ballots, double rlimit,
    double error_rate)
{
    double clean = 1.0/(2.0-margin);
    double p = 1;
    int j = 0;

    Doubles x(max_ballots, clean);
    while (p > rlimit && j <= max_ballots){
        j += 1;
        p = kaplan_kolmogorov(x, j, max_ballots, 0.5, 0);
    }

    return j;
}

bool RevCompareAudit(const AuditSpec &a1, const AuditSpec &a2){
    return a1.asn > a2.asn;
}

int GetFirstCandidateIn(const Ints &prefs, const Ints &relevant){
	for(int k = 0; k < prefs.size(); ++k){
		if(find(relevant.begin(), relevant.end(),prefs[k]) != relevant.end()){
			return prefs[k];
		}
	}
	return -1;
}

double EstimateASN_VIABLE(const Contest &ctest, int c, const Ints &tallies, 
    const Parameters &params)
{
    // We are checking that candidate 'c' is viable (tally >= 15%+1) in
    // the setting where 'elim' are eliminated.
    if(tallies[c] <= ctest.threshold)
        return -1;

    const double share = 1.0/(2*ctest.threshold_fr);
    double assorter_total = tallies[c]*share;

    double margin = 2*(assorter_total/params.tot_auditable_ballots) - 1;
    return estimate_sample_size(margin, params.tot_auditable_ballots, 
        params.risk_limit, 0);
}

double EstimateASN_NONVIABLE(const Contest &ctest, int c, const Ints &tallies,
    const Parameters &params)
{
    // We are checking that candidate 'c' is not viable (tally < 15%+1) in
    // the setting where 'elim' are eliminated. We frame this assertion as
    // having a winner (all candidates other than 'c') and a loser 'c'. 
    // We are checking whether the winner has >= 85% votes.
    if(tallies[c] >= ctest.threshold)
        return -1;

    const double share = 1.0/(2*(1 - ctest.threshold_fr));
    double assorter_total = 0;
    double total_tally = 0;
    for(int i = 0; i < tallies.size() && i < c; ++i){
        assorter_total += tallies[i] * share;
        total_tally += tallies[i];
    }
    for(int i = c+1; i < tallies.size(); ++i){
        assorter_total += tallies[i] * share;
        total_tally += tallies[i];
    }

    if(total_tally <= (1-ctest.threshold_fr)*params.tot_auditable_ballots){
        return -1;
    }

    double margin = 2*(assorter_total/params.tot_auditable_ballots) - 1;
    return estimate_sample_size(margin, params.tot_auditable_ballots, 
        params.risk_limit, 0);
}

double FindBestIRV(const Contest &ctest, const Ints &tail, 
    const SInts &winners, const Parameters &params, const Ints &tallies, 
    AuditSpec &best_audit)
{
	// Compute ASN to show that tail[0] beats one of tail[1..n] or winners. 
	const int winner = tail[0];
    best_audit.winner = winner;

    // Expand tail with candidates in winners
    Ints exp_tail(tail);
    for(SInts::const_iterator cit=winners.begin();cit!=winners.end();++cit){
        exp_tail.push_back(*cit);
    }

	double smallest = -1;

    // Estimation has been reworked to take into account
    // that we could sample ballots that do not involve this
    // contest (ie. tot_auditable_ballots >= rep_ballots.size()
	for(int i = 1; i < exp_tail.size(); ++i){
		// exp_tail[i] is the "loser"
		const int taili = exp_tail[i];
		const double V = (tallies[winner] - tallies[taili]);
		if(V <= 0) continue;

        // Note this is the diluted margin
		double margin = V/params.tot_auditable_ballots;
        double candasn = estimate_sample_size(margin, 
            params.tot_auditable_ballots, params.risk_limit, 0);

		if(smallest == -1 || candasn < smallest){
			best_audit.asn = candasn;
			best_audit.loser = taili;

			smallest = candasn;
		}
	}

	return smallest;
}

