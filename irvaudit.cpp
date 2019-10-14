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


#include<iostream>
#include<fstream>
#include<string.h>
#include<stdlib.h>
#include<cstdlib>
#include<algorithm>
#include<list>
#include<cmath>
#include<boost/property_tree/ptree.hpp>
#include<boost/property_tree/json_parser.hpp>

#include "model.h"
#include "sim_irv.h"
#include "audit.h"

using namespace std;
using boost::property_tree::ptree;

void print_list(const Ints &list){
	for(int i = 0; i < list.size(); ++i){
		cout << list[i] << " ";
	}
}

struct SimpleNode{
	Ints tail;
	double estimate;
	AuditSpec best_audit;
};

struct Node{
	Ints tail;
	double estimate;
	AuditSpec best_audit;

	bool has_ancestor;
	SimpleNode best_ancestor;

	bool expandable;
};

struct OrderFacts{
	Ints after;
	Audits after_audits;
};

typedef list<Node> Frontier;

SimpleNode CreateSimpleNode(const Node &n){
	SimpleNode sn;
	sn.tail = n.tail;
	sn.estimate = n.estimate;
	sn.best_audit = n.best_audit;
	return sn;
}

Node CreateNode(const SimpleNode &sn){
	Node newn;
	newn.tail = sn.tail;
	newn.estimate = sn.estimate;
	newn.best_audit = sn.best_audit;
	newn.has_ancestor = false;
	newn.expandable = false;
	return newn;
}

bool AppearsBefore(int winner, int loser, const Ints &tail){
    int l_index = -1;
    int w_index = -1;
    for(int i = tail.size()-1; i >= 0; --i){
        if(tail[i] == loser){
            l_index = i;
        }

        if(tail[i] == winner){
            w_index = i;
        }
    }

    if(w_index < l_index)
        return true;

    return false;
}

bool Subsumes(const AuditSpec &a1, const AuditSpec &a2){
    if(!a1.wonly || a2.wonly)
        return false;

    if(a1.wonly && (a1.winner == a2.winner && a1.loser == a2.loser))
        return true;

    if(a1.wonly && AppearsBefore(a1.winner, a1.loser, a2.rules_out)){
        return true;
    }
    
    return false;
}


bool AreAuditsEqual(const AuditSpec &a1, const AuditSpec &a2){
	if(a1.winner != a2.winner || a1.loser != a2.loser)
		return false;

	if(a1.eliminated != a2.eliminated)
		return false;

	if(a1.wonly != a2.wonly)
		return false;

	return true;
}

bool AreNodesEqual(const Node &n1, const Node &n2){
	if(n1.tail[0] != n2.tail[0])
		return false;

	if(n1.tail.size() != n2.tail.size())
		return false;

	for(int i = 1; i < n1.tail.size(); ++i){
		if(find(n2.tail.begin(), n2.tail.end(), n1.tail[i]) ==
			n2.tail.end()){
			return false;
		}	
	}
	return true;
}

bool DescendantOf(const Node &d, const Node &a){
	if(d.tail.size() <= a.tail.size())
		return false;

	const int diff = d.tail.size() - a.tail.size();
	for(int s = diff; s < d.tail.size(); ++s){
		if(d.tail[s] != a.tail[s-diff])
			return false;
	}
	
	return true;
}

void InsertNode(Frontier &front, const Node &node){		
	if(!node.expandable){
		front.insert(front.end(), node);
		return;
	} 

	Frontier::iterator it = front.begin();
	if(node.estimate == -1){
		front.insert(it, node);
		return;
	}
	for( ; it != front.end(); ++it){
		if(it->estimate == -1)
			continue;

		if(it->estimate <= node.estimate)
			break;
	}
	
	front.insert(it, node);	
}

void PrintNode(const Node &n, const Candidates &cand){
	cout << cand[n.tail[0]].id << " | ";
	for(int i = 1; i < n.tail.size(); ++i){
		cout << cand[n.tail[i]].id << " ";
	}
	cout << " [";
	cout << ((n.estimate == -1) ? -1 : n.estimate*100) << "] ";
	if(n.has_ancestor){
		cout << " (Best Ancestor ";
		cout << cand[n.best_ancestor.tail[0]].id << " | ";
		for(int i = 1; i < n.best_ancestor.tail.size(); ++i){
			cout << cand[n.best_ancestor.tail[i]].id << " ";
		}
		cout << " [" << ((n.best_ancestor.estimate == -1) ? -1 : n.best_ancestor.estimate*100) << "])";
	}
}


void OutputToJSON(const Audits &aconfig, const Candidates &cand, 
    const char *json_file, int num_ballots){
    try{
        ptree pt;
        ptree children;

        double maxasn = 0;
        for(int i = 0; i < aconfig.size(); ++i){
            const AuditSpec &spec = aconfig[i];
            ptree child;
            child.put("Winner", cand[spec.winner].id);
            child.put("Loser", cand[spec.loser].id);
            ptree aelim;
            for(int j = 0; j < spec.eliminated.size(); ++j){
                ptree c;
                c.put("", cand[spec.eliminated[j]].id);
                aelim.push_back(std::make_pair("", c));
            }
            child.add_child("Already-Eliminated", aelim);
            child.put("Winner-Only", spec.wonly);

            stringstream ss;
            if(spec.wonly){
                ss << "Rules out case where " << cand[spec.winner].id << 
                    " is eliminated before " << cand[spec.loser].id;
            }
            else{
                ss << "Rules out outcomes with tail [...";
                for(int j = 0; j < spec.rules_out.size(); ++j){
                    ss << " " << cand[spec.rules_out[j]].id;
                } 
                ss << "]";
            }
            child.put("Winner-Only", spec.wonly);
            child.put("Explanation", ss.str());
            children.push_back(std::make_pair("", child));
            maxasn = max(maxasn, spec.asn);
        }
        pt.put("Expected Polls", ceil(num_ballots*min(1.0,maxasn)));
        pt.add_child("Assertions", children);
        write_json(json_file, pt);
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
        throw STVException("Something went wrong when outputing to JSON");
	}
}

void PrintAudit(const AuditSpec &audit, const Candidates &cand){
	if(audit.wonly){
		cout << "WO,";
	}
	else{
		cout << "EO,";
	}
	cout << "Winner," << cand[audit.winner].id << ",Loser," << 
		cand[audit.loser].id << ",Eliminated";

	for(int i = 0; i < audit.eliminated.size(); ++i){
		cout << "," << cand[audit.eliminated[i]].id;
	}

    if(audit.wonly){
        cout << ",Rules out case where " <<  cand[audit.winner].id 
            << " is eliminated before " << cand[audit.loser].id << endl;
    }
    else{
        cout << ",Rules out outcomes with a tail ";
        for(int i = 0; i < audit.rules_out.size(); ++i){
            cout << cand[audit.rules_out[i]].id << " ";
        }
	    cout << endl;
    }
	//cout << "," << audit.asn*100 << endl;
}

void PrintFrontier(const Frontier &front, const Candidates &cand){
	for(Frontier::const_iterator it = front.begin(); it != front.end(); ++it){
		cout << "> ";
		PrintNode(*it, cand);
		cout << endl;
	}
}


double FindBestAudit(const Candidates &candidates, const Ballots &rep_ballots, double rlimit,
	AuditSpec &best_audit, const vector<OrderFacts> &facts, const Ints &tail, 
	double gamma, double lambda, bool alglog){
	const int i = tail[0];
	double best_estimate = -1;

	for(int j = 1; j < tail.size(); ++j){ 
		// Can we show that tail[j] comes before 'i'?
		const int tt = tail[j];
		const OrderFacts &ofacts_tt = facts[tt];
		for(int k = 0; k < ofacts_tt.after.size(); ++k){
			if(ofacts_tt.after[k] == i){
				if(best_estimate == -1 || ofacts_tt.after_audits[k].asn < best_estimate){
					best_audit = ofacts_tt.after_audits[k];
					best_estimate = ofacts_tt.after_audits[k].asn;
				}	
				break;
			}
		}
	}
	for(int k = 0; k < candidates.size(); ++k){
		if(find(tail.begin(),tail.end(), k) != tail.end()){
			continue;
		}

		// Can we show that 'i' actually comes before 'k' (who is not
		// currently in the tail!)
		const OrderFacts &ofacts_i = facts[i];
		for(int l = 0; l < ofacts_i.after.size(); ++l){
			if(ofacts_i.after[l] == k){
				if(best_estimate == -1 || ofacts_i.after_audits[l].asn < best_estimate){
					best_audit = ofacts_i.after_audits[l];
					best_estimate = ofacts_i.after_audits[l].asn;
				}	
				break;
			}
		} 
	}

	// Look at alternate ways of disproving hypothesis
	AuditSpec alternate;
    alternate.rules_out = tail;
	alternate.wonly = false;
	alternate.asn = EstimateSampleSize(rep_ballots, candidates,
		rlimit, tail, alternate, gamma, lambda);

	if(alternate.asn != -1 && (best_estimate == -1 || alternate.asn < best_estimate)){
		best_audit = alternate;
		best_estimate = alternate.asn;
	} 

	return best_estimate;
}

void ReplaceWithBestAncestor(Frontier &front, const Node &newn, 
	const Candidates &candidates, bool alglog){

	Node newnode = CreateNode(newn.best_ancestor);
	if(alglog){
		cout << "Replacing descendants of: ";
		PrintNode(newnode, candidates);
		cout << endl;
	}	

	int remcntr = 0;
	Frontier::iterator ft = front.begin();
	while(ft != front.end()){
		if(!ft->expandable){
			break;
		}
		if(DescendantOf(*ft, newnode)){
			if(alglog){
				cout << "    Removing node: ";
				PrintNode(*ft, candidates);
				cout << endl;
			}
			front.erase(ft++);
			++remcntr;
		}
		else{
			++ft;	
		}
	}

	InsertNode(front, newnode);

	if(alglog){
		cout << remcntr << " nodes replaced." << endl;
	}
}


double PerformDive(const Node &toexpand, const Candidates &cands, 
	const Ballots &rep_ballots, double rlimit, const vector<OrderFacts> &facts,
	double gamma, double lambda){
	const int ncands = cands.size();
	for(int i = 0; i < ncands; ++i){
		if(find(toexpand.tail.begin(), toexpand.tail.end(), i) ==
			toexpand.tail.end()){
					
			Node newn;
			newn.tail.push_back(i);
			for(int j = 0; j < toexpand.tail.size(); ++j){
				newn.tail.push_back(toexpand.tail[j]);
			}

			newn.estimate = -1;
			newn.best_audit.wonly = false;
			newn.best_audit.asn = -1;
			newn.estimate = -1;
			newn.expandable = (newn.tail.size() == ncands) ? false : true;

			// Set new nodes best ancestor
			if(toexpand.has_ancestor){
				if((toexpand.estimate == -1) ||
					(toexpand.best_ancestor.estimate != -1 &&
					toexpand.best_ancestor.estimate <= toexpand.estimate)){
					newn.best_ancestor = toexpand.best_ancestor;
				}
				else{
					newn.best_ancestor = CreateSimpleNode(toexpand); 
				}
			} 
			else{
				newn.best_ancestor = CreateSimpleNode(toexpand); 
			}
			newn.has_ancestor = true;
			newn.estimate = FindBestAudit(cands, rep_ballots, rlimit,
				newn.best_audit, facts, newn.tail, gamma, lambda, false);

			if(!newn.expandable){
				bool replace = false;
				if(newn.estimate == -1){
					if(!newn.has_ancestor || newn.best_ancestor.estimate == -1){
						// Audit is not possible.
						return -1;
					}
					replace = true;
				}
				else if(newn.best_ancestor.estimate != -1 &&
					newn.best_ancestor.estimate <= newn.estimate){
					replace = true;
				}
				if(replace){
					return newn.best_ancestor.estimate;
				}
				else{
					return newn.estimate;
				}
			}
			else{
				return PerformDive(newn, cands, rep_ballots, rlimit, facts,
					gamma, lambda); 
			}
			break;
		}
	}
	return -2; 
}


/**
 * Summary of command line options:
 * 
 * User can design their own reported vs actual ballot files, or provide the same file for
 * both and values for -eseed and -eprob. In the latter case, errors will be automatically
 * added to the reported ballots (using eseed and eprob) after they are read from file.
 *
 * This program is designed to both find an audit specification for a given election, and 
 * run the audit -- depending on the commanline arguments specified. When used to find an
 * audit specification, it will be printed to stdout. 
 *
 * -rep_ballots FILE     File containing reported ballots (electronic records)
 * -act_ballots FILE     File containing actual ballots (for comparison with reported ballots)
 * -simlog               If present, log messages will be printed throughout initial 
 *                           simulation of election (on basis of reported ballots)
 *
 * -gamma VALUE          Gamma parameter (for comparison audit)
 * -lambda VALUE         Lambda parameter (for comparison audit)
 *
 * -agap VALUE           This program finds a set of facts that require the least 
 *                           number of anticipated ballot polls to audit (via a comparison
 *                           audit). The implemented algorithm is based on branch-and-bound, 
 *                           and will terminate once the difference between largest valued
 *                           node on the frontier and a running lower bound on required effort
 *                           (ballot polls) to audit the given election is less than
 *                           (or equal to) agap. 
 *
 * -r VALUE              Risk limit (e.g., 0.05 represents a risk limit of 5%)
 * -dive                 If present, diving optimisation will be performed during search
 *                           for best configuration of facts to audit. RECOMMEND (default).
 *
 * -alglog               If present, log messages designed to indicate how the algorithm is
 *                           progressing will be printed.
 * -s VALUE              Seed for randomly sampling of ballots (when running audits).
 * -eseed VALUE          Seed used to randomly inject errors into reported ballots.
 * -eprob VALUE          Probability (0 to 1) with which errors will be injected into a ballot.
 *
 * -runlog               If present, log messages will be printed when an audit is run.
 * -run FILE             Given a file containing an audit specification, the program will
 *                           run the audit (rather than find an audit specification).
 *
 * -json FILE           When using the program to generate an audit, this option specifies
 *                           that the audit configuration should be output to a json file.                      
 * */
int main(int argc, const char * argv[]) 
{
	try{
		Candidates candidates;
		Ballots rep_ballots; 
		Ballots act_ballots; 
		Config config;

		bool simlog = false;
		bool alglog = false;
		double rlimit = 0.05;
		bool runAudits = false;

		vector<AuditSpec> torun;
		bool runlog = false;
		double gamma = 1.1;
		double lambda = 0;
		long seed = 845737573645710;
		long error_seed = 1837584657664;
		double error_prob = 0;

		bool diving = true;

		double lowerbound = -10;
		const char *rep_blts_file = NULL;
		const char *act_blts_file = NULL;
		const char *audit_spec_file = NULL;
        const char *json_output = NULL;

		double allowed_gap = 0.005;

		for(int i = 1; i < argc; ++i){
			if(strcmp(argv[i], "-rep_ballots") == 0 && i < argc-1){
				rep_blts_file = argv[i+1];
				++i;
			}
			else if(strcmp(argv[i], "-act_ballots") == 0 && i < argc-1){
				act_blts_file = argv[i+1];
				++i;
			}
			else if(strcmp(argv[i], "-simlog") == 0){
				simlog = true;
			}
			else if(strcmp(argv[i], "-gamma") == 0 && i < argc-1){
				gamma = ToType<double>(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-lambda") == 0 && i < argc-1){
				lambda = ToType<double>(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-agap") == 0 && i < argc-1){
				allowed_gap = ToType<double>(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-r") == 0 && i < argc-1){
				rlimit = ToType<double>(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-dive") == 0){
				diving = true;
			}
			else if(strcmp(argv[i], "-alglog") == 0){
				alglog = true;
			}
			else if(strcmp(argv[i], "-s") == 0 && i < argc-1){
				seed = stol(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-eseed") == 0 && i < argc-1){
				error_seed = stol(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-eprob") == 0 && i < argc-1){
				error_prob = ToType<double>(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-runlog") == 0){
				runlog = true;
			} 
			else if(strcmp(argv[i], "-run") == 0 && i < argc - 1){
				runAudits = true;
				audit_spec_file = argv[i+1];
                ++i;
			}
			else if(strcmp(argv[i], "-json") == 0 && i < argc - 1){
				json_output = argv[i+1];
                ++i;    
			}

		}

		srand(error_seed);
		if(!ReadReportedBallots(rep_blts_file,rep_ballots,candidates,
			config,error_prob)){
			cout << "Reported ballots read error. Exiting." << endl;
			return 1;
		}
		if(act_blts_file != NULL && !ReadActualBallots(act_blts_file,
            act_ballots,candidates,config)){
			cout << "Actual ballots read error. Exiting." << endl;
			return 1;
		}

		if(runAudits && !LoadAudits(audit_spec_file, torun, candidates, config)){
			cout << "Audit specs read error. Exiting." << endl;
			return 1;
		}
		
		for(int i = 0; i < rep_ballots.size(); ++i){
            if(rep_ballots[i].prefs.empty())
                continue;

			candidates[rep_ballots[i].prefs[0]].ballots.push_back(i);
		}

		int ncands = candidates.size();

		// Simulate IRV election to determine elimination order, 
		// and max number of candidates to batch eliminate in first round.
		Ints order_c;
		int winner = -1;
		SimIRV(rep_ballots,winner,candidates,config,order_c,simlog);

		for(int i = 0; i < ncands; ++i){
			candidates[i].Reset();
		}

		if(runAudits){
			int maxpolls_req = -1;
			Ints poll_order;
			select_random_ballots(rep_ballots.size(), seed, poll_order);
			int success = 0;

			for(int j = 0; j < torun.size(); ++j){
				const AuditSpec &spec = torun[j];
				
				for(int i = 0; i < ncands; ++i){
					candidates[i].Reset();
				}

				if(spec.wonly){
					if(runlog){
						cout<<"======================WONLY======================="<<endl;
						cout << "Winning candidate: "<<candidates[spec.winner].id <<  
							" (" << spec.winner << ")" << endl;
						cout << "Against losing candidate: "<<candidates[spec.loser].id << 
							" (" << spec.loser << ")" << endl;
					}
					Result r = RunSingleWinnerLoserAudit(act_ballots, rep_ballots,
						candidates,rlimit,spec.winner,spec.loser,poll_order,gamma,lambda);
					if(runlog){
						cout << r.polls << " (" << r.polls/config.totalvotes << 
							") ballot polls required" << endl;
					}

					if(r.remaining_hypotheses > 0){
						if(runlog){
							cout << "Audit reached/exceeed maximum ballot lookups" << endl;
							cout<<"=================================================="<<endl;
						}
						maxpolls_req = rep_ballots.size();
						continue;
					}
					success += 1;
					maxpolls_req = max(maxpolls_req, r.polls);
					if(runlog){
						cout<<"=================================================="<<endl;
					}
				}
				else{
					// Simulate elimination of appropriate candidates
					for(int k = 0; k < spec.eliminated.size(); ++k){
						SimulateElimination(spec.eliminated[k], rep_ballots, candidates);
					}
					Ints winners;
					Ints losers;
					winners.push_back(spec.winner);
					losers.push_back(spec.loser);

					if(runlog){
						cout<<"================================================="<<endl;
						cout << "Auditing that "<<candidates[spec.loser].id<< " (" << 
							spec.loser << ") lost and " << candidates[spec.winner].id << 
							" (" << spec.winner << ") won " << endl;
					}

					Result r = RunAudit(act_ballots, rep_ballots, candidates, rlimit, 
						winners, losers, poll_order, gamma, lambda);
					if(runlog){
						cout << r.polls << " (" << r.polls/config.totalvotes << 
							") ballot polls required" << endl;
					}
					if(r.remaining_hypotheses > 0){
						cout << "Audit reached/exceeed maximum ballot lookups" << endl;
						cout <<"=================================================="<<endl;
						maxpolls_req = rep_ballots.size();
						continue;
					}
				
					success += 1;
					maxpolls_req = max(maxpolls_req, r.polls);
					if(runlog){
						cout<<"=================================================="<<endl;
					}
				}
			}

			cout << "TOTAL BALLOTS POLLED = " << maxpolls_req << " (" << 
				100*(maxpolls_req/config.totalvotes) << "%)" << endl;
			cout << "SUCCESSFUL AUDITS = " << success << "/" << torun.size() << endl;
			return 0;	
		}

		mytimespec tstart;
		GetTime(&tstart);

		Frontier front;

		// Compute starting facts
		vector<OrderFacts> facts(ncands);
		for(int i = 0; i < ncands-1; ++i){
			const int ci = order_c[i];
			for(int j = i+1; j < ncands; ++j){
				const int cj = order_c[j];

				// Can we prove that ci is eliminated before cj. If so,
				// add to list of ordering facts (with best audit).
				double asn=EstimateASN_WONLY(rep_ballots,candidates,rlimit,cj,ci,
					gamma, lambda);
				if(asn != -1){
					AuditSpec aspec;
                    aspec.rules_out.push_back(ci);
					aspec.asn = asn;
					aspec.winner = cj;
					aspec.loser = ci;
					aspec.wonly = true;

					facts[ci].after.push_back(cj);
					facts[ci].after_audits.push_back(aspec);
				}
			}
		}

		bool auditfailed = false;
		int nodesexpanded = 0;

		// Create initial frontier.
		if(alglog){
			cout << "Constructing initial frontier" << endl;
		}

		for(int i = 0; i < ncands; ++i){
			if(i == winner)
				continue;

			for(int j = 0; j < ncands; ++j){
				if(i == j) continue;

				Node newn;
				newn.tail.push_back(j);
				newn.tail.push_back(i);
				newn.best_audit.wonly = false;
				newn.best_audit.asn = -1;
				newn.estimate = -1;
				newn.expandable = (ncands > 2) ? true: false;
				newn.has_ancestor = false;
				if(alglog){
					cout << "TESTING ";
					cout << candidates[newn.tail[0]].id << " | ";
					for(int i = 1; i < newn.tail.size(); ++i){
						cout << candidates[newn.tail[i]].id << " ";
					}
					cout << endl;
				}

				newn.estimate = FindBestAudit(candidates, rep_ballots, rlimit,
					newn.best_audit, facts, newn.tail, gamma, lambda, alglog);

				if(alglog){
					if(newn.estimate != -1){
						cout << "   Best audit ";
						PrintAudit(newn.best_audit, candidates);
						cout << endl;
					}
					else{
						cout << "   Cannot be disproved." << endl << endl;
					}
				}
				InsertNode(front, newn);
			}
		}


		if(alglog){
			cout << "============================================" << endl;
			cout << "Initial Frontier:" << endl;
			PrintFrontier(front, candidates);
			cout << "============================================" << endl;
		}


		while(true && !auditfailed){
			if(lowerbound > 0 && allowed_gap > 0){
				double max_on_frontier = -1;
				for(Frontier::const_iterator it = front.begin(); it != front.end(); ++it){
					if(it->estimate == -1){
						max_on_frontier = -1;
						break;
					}
					else{
						max_on_frontier = max(it->estimate, max_on_frontier);
					}
				}
				if(max_on_frontier != -1 && max_on_frontier - lowerbound <= allowed_gap){
					break;
				}
			}

			// Expand node with highest ASN (-1 == infinity)
			Node toexpand = front.front();
			if(!toexpand.expandable){
				break;
			}

			front.pop_front();

			if((toexpand.has_ancestor && toexpand.best_ancestor.estimate != -1 &&
				toexpand.best_ancestor.estimate <= lowerbound)){
				// Replace all descendents of best ancestor with ancestor.
				ReplaceWithBestAncestor(front, toexpand, candidates, alglog);
				continue;
			}
			else if(toexpand.estimate != -1 && toexpand.estimate <= lowerbound){
				// Don't expand, just make "unexpandable" and move to back of list.
				toexpand.expandable = false;
				front.insert(front.end(), toexpand);
				continue;
			} 

			if(diving){
				double divelb = PerformDive(toexpand, candidates, 
					rep_ballots, rlimit, facts, gamma, lambda);
				if(divelb == -1){
					// Audit not possible
					if(alglog) cout << "Diving finds that audit is not possible." << endl;
					auditfailed = true;
					break;
				}
				else if(divelb != -2){
					if(alglog){ 
						cout << "Diving LB " << divelb << " current LB " << lowerbound << endl;
					}
					lowerbound = max(lowerbound, divelb);
				}	
				

				if((toexpand.has_ancestor && toexpand.best_ancestor.estimate != -1 &&
					toexpand.best_ancestor.estimate <= lowerbound)){
					// Replace all descendents of best ancestor with ancestor.
					ReplaceWithBestAncestor(front, toexpand, candidates, alglog);
					continue;
				}
				else if(toexpand.estimate != -1 && toexpand.estimate <= lowerbound){
					// Don't expand, just make "unexpandable" and move to back of list.
					toexpand.expandable = false;
					front.insert(front.end(), toexpand);
					continue;
				}
			}


			++nodesexpanded;
			if(alglog){ 
				cout << " Expanding node ";
				PrintNode(toexpand, candidates);
				cout << endl;
			}
			// For each candidate 'c' not in toexpand.tail, create a new node
			// with node.tail = [c] ++ toexpand.tail
			for(int i = 0; i < ncands; ++i){
				if(find(toexpand.tail.begin(), toexpand.tail.end(), i) ==
					toexpand.tail.end()){
					
					Node newn;
					newn.tail.push_back(i);
					for(int j = 0; j < toexpand.tail.size(); ++j){
						newn.tail.push_back(toexpand.tail[j]);
					}

					newn.estimate = -1;
					newn.best_audit.wonly = false;
					newn.best_audit.asn = -1;
					newn.estimate = -1;
					newn.expandable = (newn.tail.size() == ncands) ? false : true;
				
					if(alglog){
						cout << "TESTING ";
						cout << candidates[newn.tail[0]].id << " | ";
						for(int i = 1; i < newn.tail.size(); ++i){
							cout << candidates[newn.tail[i]].id << " ";
						}
						cout << endl;
					}

					// Set new nodes best ancestor
					if(toexpand.has_ancestor){
						if((toexpand.estimate == -1) ||
							(toexpand.best_ancestor.estimate != -1 &&
							toexpand.best_ancestor.estimate <= toexpand.estimate)){
							newn.best_ancestor = toexpand.best_ancestor;
						}
						else{
							newn.best_ancestor = CreateSimpleNode(toexpand); 
						}
					} 
					else{
						newn.best_ancestor = CreateSimpleNode(toexpand); 
					}
					newn.has_ancestor = true;
					newn.estimate = FindBestAudit(candidates, rep_ballots, rlimit,
						newn.best_audit, facts, newn.tail, gamma, lambda, alglog);

					if(!newn.expandable){
						bool replace = false;
						if(newn.estimate == -1){
							if(!newn.has_ancestor || newn.best_ancestor.estimate == -1){
								// Audit is not possible.
								auditfailed = true;
								break;
							}
							replace = true;
						}
						else if(newn.best_ancestor.estimate != -1 &&
							newn.best_ancestor.estimate <= newn.estimate){
							replace = true;
						}
						if(replace){
							// Replace all descendents of newn.best_ancestor with the ancestor
							// and make expandable = false
							lowerbound = max(lowerbound, newn.best_ancestor.estimate);
							ReplaceWithBestAncestor(front, newn, candidates, alglog);
						}
						else{
							if(alglog){
								cout << "   Best audit ";
								PrintAudit(newn.best_audit, candidates);
								cout << endl;
							}

							InsertNode(front, newn);
							lowerbound = max(lowerbound, newn.estimate);
						} 
					}
					else{
						if(alglog){
							if(newn.estimate != -1){
								cout << "   Best audit ";
								PrintAudit(newn.best_audit, candidates);
								cout << endl;
							}
							else{
								cout << "   Cannot be disproved." << endl;
							}
						}

						InsertNode(front, newn);
					}
				}
			}
			if(auditfailed){
				break;
			}  
			if(alglog){
				cout << endl << "Size of frontier " << front.size() << ", Nodes expanded "
					<< nodesexpanded << ", Current threshold " << lowerbound <<  endl << endl;
			}
		}

		mytimespec tend;
		GetTime(&tend);

		double maxasn = -1;
		if(!auditfailed){
			// Compile list of audits to complete.
			Audits audits;
			for(Frontier::const_iterator it = front.begin(); it != front.end(); ++it){
				bool pruneaudit = false;
				for(Audits::const_iterator at = audits.begin(); at != audits.end(); ++at){
					if(AreAuditsEqual(*at, it->best_audit)){
						pruneaudit = true;
						break;
					}
				}	
				if(!pruneaudit){
					audits.push_back(it->best_audit);
				}
			}

			cout << "============================================" << endl;
			cout << "AUDITS REQUIRED" << endl;
			maxasn = 0;
            Audits final_config;
			for(Audits::const_iterator it = audits.begin(); it != audits.end();++it){
                bool subsumed = false;
                for(Audits::const_iterator jt = audits.begin(); jt != audits.end(); ++jt){
                    if(jt == it) continue;
                    if(Subsumes(*jt, *it)){
                        subsumed = true;
                        break;
                    }
                }
                if(!subsumed){
                    final_config.push_back(*it);
				    PrintAudit(*it, candidates);
				    maxasn = max(maxasn, it->asn);
                }
			}
			maxasn *= 100;
			cout << "MAX ASN(%) " << maxasn << endl;
			cout << "============================================" << endl;
            if(json_output != NULL){
                OutputToJSON(final_config, candidates, json_output, config.totalvotes);
            }
		}
		else{
			if(alglog){
				cout << endl;
				cout << "AUDIT NOT POSSIBLE" <<endl;
			}	
		}
		cout << "TIME," << tend.seconds - tstart.seconds << ",Nodes Expanded," << nodesexpanded
			<< ",MAX ASN(%)," << maxasn << endl;
	}
	catch(exception &e)
	{
		cout << e.what() << endl;
		cout << "Exiting." << endl;
		return 1;
	}
	catch(STVException &e)
	{
		cout << e.what() << endl;
		cout << "Exiting." << endl;
		return 1;
	}	
	catch(...)
	{
		cout << "Unexpected error. Exiting." << endl;
		return 1;
	}

	return 0;
}




