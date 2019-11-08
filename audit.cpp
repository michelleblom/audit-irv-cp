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

#include "audit.h"
#include<fstream>
#include<iostream>
#include<cmath>

using namespace std;
typedef boost::char_separator<char> boostcharsep;

bool RevCompareAudit(const AuditSpec &a1, const AuditSpec &a2) {
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

double EstimateASN_WONLY(const Ballots &rep_ballots, 
	const Candidates &cand, const Parameters &params, int winner, int loser)
{
	double total_votes_loser = 0;
	double total_votes_winner = cand[winner].sim_votes;

	Ints relevant;
	relevant.push_back(winner);
	relevant.push_back(loser);
	for(Ballots::const_iterator bt = rep_ballots.begin();
		bt != rep_ballots.end(); ++bt){
		int nextcand = GetFirstCandidateIn(bt->prefs, relevant);
		if(nextcand == loser){
			total_votes_loser += 1;
		}
	}

	const double V = (total_votes_winner - total_votes_loser);
	if(V <= 0) return -1;

    // Estimatation has been reworked to take into account
    // that we could sample ballots that do not involve this
    // contest (ie. tot_auditable_ballots >= rep_ballots.size()
	const double total_votes_present=params.tot_auditable_ballots; 
    // margin is the 'diluted margin'
	const double margin=V/total_votes_present;
	
    //const double od2g=1.0/(2.0 * params.gamma);
	//const double row=-log(params.risk_limit)/(od2g+params.lambda*log(1-od2g));

    // Note: we are returning a proportion of the total number
    // of ballots involved in the audit (across all contests)
	//return ceil(row/margin)/total_votes_present;

    // NOTE THIS HAS BEEN MODIFIED FOR PHILIP's SHARPER ESTIMATION SCHEME
    // THAT DOES NOT USE LAMBDA/GAMMA PARAMETERS.
    return ceil(1.0/margin)/total_votes_present;
}

double EstimateSampleSize(const Ballots &rep_ballots, const Candidates &cand, 
	const Parameters &params, const Ints &tail, AuditSpec &best_audit){
	// Compute ASN to show that tail[0] beats one of tail[1..n]
	const int loser = tail[0];
	const int tsize = tail.size();
	Doubles total_votes_each(tsize, 0);
	
	for(int b = 0; b < rep_ballots.size(); ++b){
		const Ballot &bt = rep_ballots[b];
		int nextcand = GetFirstCandidateIn(bt.prefs, tail);
		for(int k = 0; k < tsize; ++k){
			if(nextcand == tail[k]){
				total_votes_each[k] += 1;
			}
		}
	}

	best_audit.eliminated.clear();
	for(int k = 0; k < cand.size(); ++k){
		if(find(tail.begin(), tail.end(), k) == tail.end()){
			best_audit.eliminated.push_back(k);
		}
	}

	double smallest = -1;
    // Estimatation has been reworked to take into account
    // that we could sample ballots that do not involve this
    // contest (ie. tot_auditable_ballots >= rep_ballots.size()
	double total_votes_present = params.tot_auditable_ballots; 
	for(int i = 1; i < tsize; ++i){
		// loser is the "winner" and tail[i] is the "loser"
		const int taili = tail[i];
		double total_votes_winner = total_votes_each[0];
		double total_votes_loser = total_votes_each[i];

		const double V = (total_votes_winner-total_votes_loser);
		if(V <= 0) continue;

        // Note this is the diluted margin
		const double margin = V/total_votes_present;

		//const double od2g = 1.0/(2 * params.gamma);
		//const double row = -log(params.risk_limit)/
        //    (od2g+params.lambda*log(1-od2g));

        // candasn is a proportion of the total number of ballots 
        //involved in the audit (across all contests)
		//const double candasn = ceil(row/margin)/total_votes_present;
    
        // NOTE THIS HAS BEEN MODIFIED FOR PHILIP's SHARPER ESTIMATION SCHEME
        // THAT DOES NOT USE LAMBDA/GAMMA PARAMETERS.
        const double candasn = ceil(1.0/margin)/total_votes_present;

		if(smallest == -1 || candasn < smallest){
			best_audit.asn = candasn;
			best_audit.winner = loser;
			best_audit.loser = taili;
			smallest = candasn;
		}
	}

	return smallest;
}

