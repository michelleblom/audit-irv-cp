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
#include<boost/property_tree/json_parser.hpp>
#include<boost/property_tree/ptree.hpp>
#include<boost/foreach.hpp>
#include<limits>
#include<algorithm>

using namespace std;

typedef boost::char_separator<char> boostcharsep;

namespace pt = boost::property_tree;

int median(vector<int> &v)
{
    size_t n = v.size() / 2;
    nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}
  
// Using Kaplan Kolgoromov
int estimate_sample_size_x(double margin, const Parameters &params, 
    mt19937_64 &gen)
{
    double clean  = 1.0/(2-margin);
    double one_vote_over = 0.5/(2-margin);

    uniform_real_distribution<> dist(0,1);

    Ints sams(params.reps, 0);

    for(int i = 0; i < params.reps; ++i){
        Doubles x(params.tot_auditable_ballots, clean);

        for(int k = 0; k < params.tot_auditable_ballots; ++k){
            double rf = dist(gen);
            if(rf < params.error_rate){
                x[k] = one_vote_over;
            }
        }

        double p = 1;
        int j = 0;

        const double t = params.t;
        const double g = params.g;
        const double N = params.tot_auditable_ballots;
        double sample_total = 0;
        double mart = (t > 0) ? (x[0]+g)/(t+g) : 1;
    
        p = min(1.0/mart,1.0);
        j += 1;
        for( ; p > params.risk_limit && j < N; ++j)
        {
            mart *= (x[j]+g)*(1-j/N)/(t+g - (1/N)*sample_total);

            if(mart < 0){
                break;
            }
            else {
                sample_total += x[j]+g;
            }
            p = min(1.0/mart,1.0);
        }

        if(p <= params.risk_limit){
            sams[i] = j;
        }
        else{
            return -1;
        }
    }


    return median(sams);
}

// Using Kaplan Kolgoromov
int estimate_sample_size(double margin, const Parameters &params)
{
    double x = 1.0/(2.0-margin);
    double p = 1;
    int j = 0;

    const double t = params.t;
    const double g = params.g;
    const double N = params.tot_auditable_ballots;

    double sample_total = 0;
    double mart = (t > 0) ? (x+g)/(t+g) : 1;
    
    p = min(1.0/mart,1.0);
    j += 1;
    for( ; p > params.risk_limit && j < N; ++j)
    {
        mart *= (x+g)*(1-j/N)/(t+g - (1/N)*sample_total);

        if(mart < 0){
            break;
        }
        else {
            sample_total += x+g;
        }
        p = min(1.0/mart,1.0);
    }

    if(p <= params.risk_limit){
        return j;
    }

    return -1;
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

double EstimateASN_cdiff(double tallyA, double tallyB, double d, 
    int exhausted, const Parameters &params, double &margin)
{
    double tvotes = params.tot_auditable_ballots;
    
    double share1 = 1.0/(1.0 + d);
    double share2 = 1.0/(2.0 + 2.0*d);
    double share3 = 0.5;

    double other = tvotes - tallyA - tallyB - exhausted;

    double assorter_total = tallyA*share1 + share2*other + share3*exhausted;
    margin = 2*(assorter_total/tvotes) - 1;
    return estimate_sample_size(margin, params);
}

double EstimateASN_smajority(double tally, double other, double threshold_fr,
    const Parameters &params, double &margin)
{
    double tvotes = params.tot_auditable_ballots;

    double share = 1.0/(2*threshold_fr);
    double assorter_total = tally*share + 0.5*other;

    margin = 2*(assorter_total/tvotes) - 1;
    return estimate_sample_size(margin, params);
}

double EstimateASN_VIABLE(const Contest &ctest, int c, const Ints &tallies, 
    int exhausted, const Parameters &params, double &margin)
{
    if(tallies[c] <= ctest.threshold)
        return -1;

    return EstimateASN_smajority(tallies[c], exhausted, ctest.threshold_fr,
        params, margin);
}

double EstimateASN_NONVIABLE(const Contest &ctest, int c, const Ints &tallies,
    int exhausted, const Parameters &params, double &margin)
{
    // We are checking that candidate 'c' is not viable (tally < 15%+1) in
    // the setting where 'elim' are eliminated. We frame this assertion as
    // having a winner (all candidates other than 'c') and a loser 'c'. 
    // We are checking whether the winner has >= 85% votes.
    if(tallies[c] >= ctest.threshold)
        return -1;

    double share = 1.0/(2*(1 - ctest.threshold_fr));
    double assorter_total = exhausted*0.5;
    double total_tally = exhausted;
    for(int i = 0; i < tallies.size() && i < c; ++i){
        assorter_total += tallies[i] * share;
        total_tally += tallies[i];
    }
    for(int i = c+1; i < tallies.size(); ++i){
        assorter_total += tallies[i] * share;
        total_tally += tallies[i];
    }

    if(total_tally <= (1-ctest.threshold_fr)*ctest.rballots.size()){
        return -1;
    }

    margin = 2*(assorter_total/params.tot_auditable_ballots) - 1;
    return estimate_sample_size(margin, params);
}

double FindBestIRV_NEB(const Contest &ctest, const Ints &tail, 
    const SInts &winners, const Parameters &params, 
    const Ints &tallies, const Audits2d &nebs, const Bools2d &has_neb,
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

        // Check NEB
        if(has_neb[winner][taili]){
            const AuditSpec &neb_wi = nebs[winner][taili];
            if(neb_wi.asn != -1 && ((smallest == -1) |
                (neb_wi.asn < smallest))){
                best_audit = neb_wi;
                smallest = neb_wi.asn;
            }
        }

        // Check IRV assertion.
        if(tallies[winner] <= tallies[taili])
            continue;

        int neither = params.tot_auditable_ballots - tallies[winner]
            - tallies[taili];

        // The assorter margin is 2 times the mean of 
        // ((winner - loser) + 1)/2 across all CVRs - 1. 
        // For each CVR, an assorter will return 1 if
        // its a vote for the winner, 0 if its a vote for the loser, and
        // 0.5 if its a vote for neither.
        double amean = (tallies[winner] + 0.5*neither)/
            params.tot_auditable_ballots;

		double margin = 2*amean - 1;
        double candasn = estimate_sample_size(margin,params);

		if(smallest == -1 || candasn < smallest){
			best_audit.asn = candasn;
			best_audit.loser = taili;
            best_audit.margin = margin;

			smallest = candasn;
		}
	}

	return smallest;
}

