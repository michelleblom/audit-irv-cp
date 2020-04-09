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

using namespace std;
typedef boost::char_separator<char> boostcharsep;


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

double EstimateASN_VIABLE(const Contest &ctest, int c, const Ints &elim,
    const Parameters &params)
{
    // We are checking that candidate 'c' is viable (tally >= 15%+1) in
    // the setting where 'elim' are eliminated.
    const double total_votes_present = params.tot_auditable_ballots; 
    int total_for_c = 0;
    for(int i = 0; i < ctest.rballots.size(); ++i){
        const Ints &prefs = ctest.rballots[i].prefs;
        for(int j = 0; j < prefs.size(); ++j){
            int pc = prefs[j];
            if(pc == c){
                total_for_c += 1;
                break;
            }
            if(find(elim.begin(), elim.end(), pc) == elim.end()){
                break;
            }
        }
    }
    if(total_for_c <= ctest.threshold)
        return -1;

    double margin = (total_for_c - ctest.threshold)/total_votes_present;
    return ceil(1.0/margin)/total_votes_present;
}

double EstimateASN_NONVIABLE(const Contest &ctest, int c, const Ints &elim,
    const Parameters &params)
{
    // We are checking that candidate 'c' is not viable (tally < 15%+1) in
    // the setting where 'elim' are eliminated.
    const double total_votes_present = params.tot_auditable_ballots; 
    int total_for_c = 0;
    for(int i = 0; i < ctest.rballots.size(); ++i){
        const Ints &prefs = ctest.rballots[i].prefs;
        for(int j = 0; j < prefs.size(); ++j){
            int pc = prefs[j];
            if(pc == c){
                total_for_c += 1;
                break;
            }
            if(find(elim.begin(), elim.end(), pc) == elim.end()){
                break;
            }
        }
    }
 
    if(total_for_c >= ctest.threshold)
        return -1;

    double margin = (ctest.threshold - total_for_c)/total_votes_present;
    return ceil(1.0/margin)/total_votes_present;
}

double FindBestIRV(const Contest &ctest, const Ints &tail, 
    const SInts &winners, const Parameters &params, AuditSpec &best_audit)
{
	// Compute ASN to show that tail[0] beats one of tail[1..n] or winners. 
	const int loser = tail[0];

    // Expand tail with candidates in winners
    Ints exp_tail(tail);
    for(SInts::const_iterator cit = winners.begin(); 
        cit != winners.end(); ++cit){
        exp_tail.push_back(*cit);
    }

    const int tailsize = exp_tail.size();
	Doubles total_votes_each(tailsize, 0);

	for(int b = 0; b < ctest.rballots.size(); ++b){
		const Ballot &bt = ctest.rballots[b];
		int nextcand = GetFirstCandidateIn(bt.prefs, exp_tail);
		for(int k = 0; k < tailsize; ++k){
			if(nextcand == exp_tail[k]){
				total_votes_each[k] += 1;
                break;
			}
		}
	}

    //cout << "FBIRV tail: ";
    //for(int i = 0; i < exp_tail.size(); ++i){
    //    cout << ctest.cands[exp_tail[i]].id << "("
    //        << total_votes_each[i] << ") ";
   // }
    //cout << endl;

	best_audit.eliminated.clear();
	for(int k = 0; k < ctest.ncandidates; ++k){
		if(find(exp_tail.begin(), exp_tail.end(), k) == exp_tail.end()){
			best_audit.eliminated.push_back(k);
		}
	}

	double smallest = -1;

    // Estimation has been reworked to take into account
    // that we could sample ballots that do not involve this
    // contest (ie. tot_auditable_ballots >= rep_ballots.size()
	double total_votes_present = params.tot_auditable_ballots; 
	for(int i = 1; i < tailsize; ++i){
		// loser is the "winner" and exp_tail[i] is the "loser"
		const int taili = exp_tail[i];
		double total_votes_winner = total_votes_each[0];
		double total_votes_loser = total_votes_each[i];

		const double V = (total_votes_winner-total_votes_loser);
		if(V <= 0) continue;

        // Note this is the diluted margin
		const double margin = V/total_votes_present;
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

