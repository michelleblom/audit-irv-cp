#include "audit.h"
#include<fstream>
#include<iostream>
#include<cmath>

const double EPSILON = 0.0000001;

using namespace std;
typedef boost::char_separator<char> boostcharsep;

void PrintBallot(const Ballot &b){
	for(int i = 0; i < b.prefs.size(); ++i){
		cout << b.prefs[i] << " ";
	}
}

bool LoadAudits(const char *path, Audits &audits,
	const Candidates &candidates, const Config &config){
	try
	{
		ifstream infile(path);

		boostcharsep sp(",");
	
		string line;	
		while(getline(infile, line))
		{
			if(!(boost::starts_with(line, "EO")||
				boost::starts_with(line,"WO"))){
				continue;
			}
	
			vector<string> columns;
			Split(line, sp, columns);

			AuditSpec spec;
			spec.wonly = false;

			if(columns[0] == "WO"){
				spec.wonly = true;
			}

			int winner = ToType<int>(columns[2]);
			int loser = ToType<int>(columns[4]);

			spec.winner = config.id2index.find(winner)->second;
			spec.loser = config.id2index.find(loser)->second;

			for(int k = 6; k < columns.size(); ++k){
				int elim = ToType<int>(columns[k]);
				spec.eliminated.push_back(config.id2index.find(elim)->second);
			} 
			audits.push_back(spec);
		}

		infile.close();
	}
	catch(exception &e)
	{
		throw e;
	}
	catch(STVException &e)
	{
		throw e;
	}
	catch(...)
	{
		cout << "Unexpected error reading in audits." << endl;
		return false;
	}

	return true;
}

int GetFirstCandidateStanding(const Candidates &cand, const Ints &prefs){
	for(Ints::const_iterator it = prefs.begin(); it != prefs.end(); ++it){
		if(cand[*it].standing == 1){
			return *it;
		}
	}
	return -1;
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
	const Candidates &cand, double rlimit, int winner, 
	int loser, double gamma, double lambda){

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

	const double total_votes_present = rep_ballots.size();
	const double margin = V/total_votes_present;
	const double od2g = 1.0/(2.0 * gamma);
	const double row = -log(rlimit)/(od2g+lambda*log(1-od2g));
	return ceil(row/margin)/total_votes_present;
}

double EstimateSampleSize(const Ballots &rep_ballots, const Candidates &cand, 
	double rlimit, const Ints &tail, AuditSpec &best_audit,
	double gamma, double lambda){
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
	double total_votes_present = rep_ballots.size();
	for(int i = 1; i < tsize; ++i){
		// loser is the "winner" and tail[i] is the "loser"
		const int taili = tail[i];
		double total_votes_winner = total_votes_each[0];
		double total_votes_loser = total_votes_each[i];

		const double V = (total_votes_winner-total_votes_loser);
		if(V <= 0) continue;

		const double margin = V/total_votes_present;

		const double od2g = 1.0/(2 * gamma);
		const double row = -log(rlimit)/(od2g+lambda*log(1-od2g));
		const double candasn = ceil(row/margin)/total_votes_present;

		if(smallest == -1 || candasn < smallest){
			best_audit.asn = candasn;
			best_audit.winner = loser;
			best_audit.loser = taili;
			smallest = candasn;
		}
	}

	return smallest;
}

Result RunSingleWinnerLoserAudit(const Ballots &act_ballots, 
	const Ballots &rep_ballots, const Candidates &cand,
	double rlimit, int winner, int loser, const Ints &plist,
	double gamma, double lambda){
	
	double total_votes_loser = 0;
	double total_votes_winner = cand[winner].sim_votes;

	Ints relevant;
	relevant.push_back(winner);
	relevant.push_back(loser);

	double total_votes_present = rep_ballots.size();
	for(Ballots::const_iterator bt = rep_ballots.begin();
		bt != rep_ballots.end(); ++bt){
		int nextcand = GetFirstCandidateIn(bt->prefs, relevant);
		if(nextcand == loser){
			total_votes_loser += 1;
		}
	}

	double pkm = 1;
	int polls = 0;
	Result r;
	r.remaining_hypotheses = 1;
	
	const double V = (total_votes_winner - total_votes_loser);
	const double margin = V/total_votes_present;

	if(margin <= 0){
		return r;
	}
	const double numerator = (1 - (1.0/((2*gamma)/margin)));
	for(int p = 0; p < plist.size(); ++p){
		++polls;
	
		const Ballot &rep_b = rep_ballots[plist[p]];
		const Ballot &act_b = act_ballots[plist[p]];
		int act_c = GetFirstCandidateIn(act_b.prefs, relevant);
		int rep_c = GetFirstCandidateIn(rep_b.prefs, relevant);

		if(rep_c == winner && rep_b.prefs[0] != winner){
			rep_c = -1;
		} 

		if(act_c == winner && act_b.prefs[0] != winner){
			act_c = -1;
		}

		cout << "Reported ";
		PrintBallot(rep_b);
		cout << "Actual ";
		PrintBallot(act_b);
		cout << "Error ";

		double opkm = pkm;

		// As V_wl == V, the denominator 1 - (error/(2*gamma/V))
		// where error given by eqn 5 in Stark's paper reduces to
		// 1 - error/2*gamma
		if(act_c == rep_c){
			pkm *= numerator;
			cout << " NA" << endl;
		}
		else if(rep_c == winner && act_c == loser){
			// error is +2 
			pkm *= (numerator/(1 - 1.0/gamma));
			cout << " +2 " << endl;
		}
		else if(rep_c == loser && act_c == winner){
			// error is -2 
			pkm *= (numerator/(1 + 1.0/gamma));
			cout << " -2 " << endl;
		}
		else if(rep_c == winner && act_c != winner && act_c != loser){
			// error is +1
			pkm *= (numerator/(1 - 1.0/(2*gamma)));
			cout << " +1 " << endl;
		} 
		else if(rep_c == loser && act_c != winner && act_c != loser){
			// error is -1
			pkm *= (numerator/(1 + 1.0/(2*gamma)));
			cout << " -1 " << endl;
		}
		else{
			pkm *= numerator;
			cout << " NA " << endl;
		}

		cout << "PKM from " << opkm << " to " << pkm << endl;
		if(pkm <= rlimit){
			r.remaining_hypotheses = 0;
			break;
		}
	}	
	r.polls = polls;
	return r;
}

Result RunAudit(const Ballots &act_ballots, const Ballots &rep_ballots,
	const Candidates &cand, double rlimit, const Ints &winners, 
	const Ints &losers,const Ints &plist, double gamma, double lambda){

	const int maxballots = rep_ballots.size();
	const int ncand = cand.size();
	
	Ints action(ncand, 0);
	double pkm = 1;

	Doubles2d margins;
	double total_votes_present = rep_ballots.size();
	for(int i = 0; i < ncand; ++i){
		margins.push_back(Doubles(ncand, 0));
	}

	// compute margin for each winner-loser pair
	double V = total_votes_present;
	for(Ints::const_iterator wt=winners.begin(); wt!=winners.end(); ++wt){
		action[*wt] = 1;
		for(Ints::const_iterator lt=losers.begin(); lt!=losers.end(); ++lt){
			double margin = cand[*wt].sim_votes-cand[*lt].sim_votes;
			margins[*wt][*lt] = margin;
			V = min(V, margin);
		}
	}
	const double mu = V/total_votes_present;
	const double numerator = (1 - 1.0/((2*gamma)/mu));
	const double denom = (2*gamma)/V;

	int polls = 0;
	Result r;
	r.remaining_hypotheses = 1;
    Ints all_relevant(winners);
	all_relevant.insert(all_relevant.end(), losers.begin(), losers.end());

	while(polls < maxballots){
		int rbidx = plist[polls];

		const Ballot &act_b = act_ballots[rbidx];
		const Ballot &rep_b = rep_ballots[rbidx];
		const int rep_c = GetFirstCandidateIn(rep_b.prefs, all_relevant);
		const int act_c = GetFirstCandidateIn(act_b.prefs, all_relevant);

		cout << "Reported ";
		PrintBallot(rep_b);
		cout << "Actual ";
		PrintBallot(act_b);
		cout << "Error ";

		double opkm = pkm;


		// Compute max error
		double max_error = -2;	
		int max_int_error = -2;	
		for(int w = 0; w < winners.size(); ++w){
			const int widx = winners[w];
			int vbw = 0, abw = 0;
			if(rep_c == widx)
				vbw = 1;
			if(act_c == widx)
				abw = 1;

			for(int l = 0; l < losers.size(); ++l){
				const int lidx = losers[l];
				// was the ballot a reported/actual vote for winner w? or
				// reported/actual vote for loser?
				int vbl = 0, abl = 0;
				if(rep_c == lidx)
					vbl = 1;
				if(act_c == lidx)
					abl = 1;
		
				double error = (vbw - abw - vbl + abl)/margins[widx][lidx];
				max_int_error = max(vbw-abw-vbl+abl,max_int_error);
				max_error = max(error, max_error);
			}
		}
		++polls;
		pkm *= (numerator/(1 - (max_error/denom)));	
		cout << "Max error " << max_error << ", int error " << max_int_error << endl;
		cout << "PKM from " << opkm << " to " << pkm << endl;

		if(pkm <= rlimit){
			r.remaining_hypotheses = 0;
			break;
		}
	}
	r.polls = polls;
	return r;	
}
