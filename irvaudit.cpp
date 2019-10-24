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
		cout << " [" << ((n.best_ancestor.estimate == -1) ? 
            -1 : n.best_ancestor.estimate*100) << "])";
	}
}


void OutputToJSON(const Contests &contests, const vector<Audits> &torun,
    int tot_auditable_ballots, const char *json_file)
{
    try{
        ptree pt;
        ptree children;

        int overall_maxasn = -1;
        for(int k = 0; k < contests.size(); ++k){
            const Contest &ctest = contests[k];
            const Audits &aconfig = torun[k];

            if(aconfig.empty())
                continue;

            ptree caudit;
            ptree caudit_children;
            double maxasn = 0;
            for(int i = 0; i < aconfig.size(); ++i){
                const AuditSpec &spec = aconfig[i];
                ptree child;
                child.put("Winner", ctest.cands[spec.winner].id);
                child.put("Loser", ctest.cands[spec.loser].id);
                ptree aelim;
                for(int j = 0; j < spec.eliminated.size(); ++j){
                    ptree c;
                    c.put("", ctest.cands[spec.eliminated[j]].id);
                    aelim.push_back(std::make_pair("", c));
                }
                child.add_child("Already-Eliminated", aelim);
                if(spec.wonly)
                    child.put("assertion_type", "WINNER_ONLY");
                else
                    child.put("assertion_type", "IRV_ELIMINATION");

                stringstream ss;
                if(spec.wonly){
                    ss << "Rules out case where " << 
                        ctest.cands[spec.winner].id << 
                        " is eliminated before " << 
                        ctest.cands[spec.loser].id;
                }
                else{
                    ss << "Rules out outcomes with tail [...";
                    for(int j = 0; j < spec.rules_out.size(); ++j){
                        ss << " " << ctest.cands[spec.rules_out[j]].id;
                    }    
                    ss << "]";
                }
                child.put("Explanation", ss.str());
                caudit_children.push_back(std::make_pair("", child));
                maxasn = max(maxasn, spec.asn);
            }
            caudit.put("Contest", ctest.id);
            caudit.put("Winner", ctest.cands[ctest.winner].id);
            ptree elim;
            for(int j = 0; j < ctest.outcome.size()-1; ++j){
                ptree c;
                c.put("", ctest.cands[ctest.outcome[j]].id);
                elim.push_back(std::make_pair("", c));
            }
            caudit.add_child("Eliminated", elim);

            int maxasn_ballots = ceil(ctest.config.totalvotes*maxasn);
            int maxasn_percent = ceil(maxasn*100);

            caudit.put("Expected Polls (#)", maxasn_ballots);
            caudit.put("Expected Polls (%)", maxasn_percent);

            caudit.add_child("Assertions", caudit_children);
            overall_maxasn = max(overall_maxasn, maxasn_ballots);
            
            children.push_back(std::make_pair("", caudit));
        }
        if(overall_maxasn != -1){
            pt.put("Overall Expected Polls (#)", overall_maxasn);
            pt.put("Ballots involved in audit (#)", tot_auditable_ballots);
            pt.add_child("Audits", children);
        }
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
    cout << endl;
	//cout << "," << audit.asn*100 << endl;
}

void PrintFrontier(const Frontier &front, const Candidates &cand){
	for(Frontier::const_iterator it = front.begin(); it != front.end(); ++it){
		cout << "> ";
		PrintNode(*it, cand);
		cout << endl;
	}
}


double FindBestAudit(const Candidates &candidates, const Ballots &rep_ballots,
    double rlimit,AuditSpec &best_audit, const vector<OrderFacts> &facts, 
    const Ints &tail, double gamma, double lambda, bool alglog,
    int tot_auditable_ballots)
{
	const int i = tail[0];
	double best_estimate = -1;

	for(int j = 1; j < tail.size(); ++j){ 
		// Can we show that tail[j] comes before 'i'?
		const int tt = tail[j];
		const OrderFacts &ofacts_tt = facts[tt];
		for(int k = 0; k < ofacts_tt.after.size(); ++k){
			if(ofacts_tt.after[k] == i){
				if(best_estimate == -1 || ofacts_tt.after_audits[k].asn 
                    < best_estimate){
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
				if(best_estimate == -1 || ofacts_i.after_audits[l].asn 
                    < best_estimate){
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
		rlimit, tail, alternate, gamma, lambda, tot_auditable_ballots);

	if(alternate.asn != -1 && (best_estimate == -1 || alternate.asn 
        < best_estimate)){
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
	const Ballots &rep_ballots, double rlimit, 
    const vector<OrderFacts> &facts, double gamma, double lambda,
    int tot_auditable_ballots)
{
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
				newn.best_audit, facts, newn.tail, gamma, lambda, false,
                tot_auditable_ballots);

			if(!newn.expandable){
				bool replace = false;
				if(newn.estimate == -1){
					if(!newn.has_ancestor||newn.best_ancestor.estimate == -1){
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
					gamma, lambda, tot_auditable_ballots); 
			}
			break;
		}
	}
	return -2; 
}


/**
 * Summary of command line options:
 * 
 * User can design their own reported vs actual ballot files, or provide the
 * same file for both and values for -eseed and -eprob. In the latter case,
 * errors will be automatically added to the reported ballots (using eseed and
 * eprob) after they are read from file.
 *
 * This program is designed to both find an audit specification for a given
 * election, and run the audit -- depending on the commanline arguments
 * specified. When used to find an audit specification, it will be printed to
 * stdout. 
 *
 * -rep_ballots FILE     File containing reported ballots (electronic records)
 * -simlog               If present, log messages will be printed throughout 
 *                           initial simulation of election (on basis of 
 *                           reported ballots)
 *
 * -gamma VALUE          Gamma parameter (for comparison audit)
 * -lambda VALUE         Lambda parameter (for comparison audit)
 *
 * -agap VALUE           This program finds a set of facts that require the 
 *                           least number of anticipated ballot polls to audit
 *                           (via a comparison audit). The implemented algorithm
 *                           is based on branch-and-bound, and will terminate 
 *                           once the difference between largest valued node on 
 *                           the frontier and a running lower bound on required
 *                           effort (ballot polls) to audit the given election
 *                           is less than (or equal to) agap. 
 *
 * -r VALUE              Risk limit (e.g., 0.05 represents a risk limit of 5%)
 *
 * -alglog               If present, log messages designed to indicate how the
 *                           algorithm is progressing will be printed.
 *
 * -json FILE            When using the program to generate an audit, this 
 *                           option specifies that the audit configuration 
 *                           should be output to a json file. 
 *
 * -contests N C1 C2 ... Number of contests for which we want to audit/generate
 *                           and audit, followed by the numeric identifiers for
 *                           those contests. IF NOTE PRESENT, ASSUME ALL
 *                           CONTESTS MENTIONED IN INPUT WILL BE AUDITED.  
 * 
 * -help                 Print usage instructions.     
 * */

void PrintUsageInstructions(){
    cout << "============================================" << endl;
    cout << "USAGE" << endl;
    cout << "./irvaudit -rep_ballots FILE.raire -gamma GAMMA ";
    cout << "-lambda LAMBDA -r RISK_LIMIT -json OUTOUT.json" << endl;
    cout << endl;
    cout << "Required Arguments" << endl;
    cout << "------------------" << endl;
    cout << "-rep_ballots FILE.raire  Reported ballots in RAIRE format"<<endl; 
    cout << "-gamma GAMMA             Gamma parameter GAMMA >= 1" << endl;   
    cout << "-lambda LAMBDA           Lambda parameter 0 < LAMBDA <= 1" << endl;
    cout << "-r RISK_LIMIT            Risk Limit (e.g., 0.05 for 5%)" << endl;
    cout << "-json FILE               Generated audit will be recorded in "
        << "JSON format in FILE" << endl; 
    cout << endl;
    cout << "Optional Arguments" << endl;
    cout << "------------------" << endl;
    cout<<"-contests N C1 C2 ...    Contests to generate an audit for."<<endl;
    cout<<"                         N  = number of contests to audit" << endl;
    cout<<"                         Ci = IDs of contests to audit" << endl;
    cout<<"                         Note that if omitted, all contests in";
    cout<<" FILE.raire will be audited." << endl;
    cout << endl;
    cout<<"-agap ALLOWED_GAP        Degree of suboptimality permitted in" <<
        " generated audit." << endl;
    cout<<"                         Default is 0.005 (0.5%)." << endl; 
    cout << endl;
    cout<<"-simlog                  Print simulation of election prior to";
    cout << " generating audit." <<endl;
    cout<<"-alglog                  Print algorithm progress while";
    cout << " generating audit." <<endl;
    cout << "============================================" << endl;
}

int main(int argc, const char * argv[]) 
{
	try
    {
        Contests contests;

		bool simlog = false;
		bool alglog = false;
		double rlimit = 0.05;

		double gamma = 1.1;
		double lambda = 0;

		bool diving = true;

		const char *rep_blts_file = NULL;
        const char *json_output = NULL;

		double allowed_gap = 0.005;

		vector<Audits> audits_to_run;

        // Contests have their own unique ids, we need to know 
        // (as we are reading in ballots), the index in 
        // "contests" that each id maps to. 
        ID2IX contest_id2index;
		for(int i = 1; i < argc; ++i){
            if(strcmp(argv[i], "-help") == 0){
                PrintUsageInstructions();
                return 0;
            }
			else if(strcmp(argv[i], "-rep_ballots") == 0 && i < argc-1){
				rep_blts_file = argv[i+1];
				++i;
			}
			else if(strcmp(argv[i], "-simlog") == 0){
				simlog = true;
			}
			else if(strcmp(argv[i], "-gamma") == 0 && i < argc-1){
				gamma = atof(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-lambda") == 0 && i < argc-1){
				lambda = atof(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-agap") == 0 && i < argc-1){
				allowed_gap = atof(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-r") == 0 && i < argc-1){
				rlimit = atof(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-alglog") == 0){
				alglog = true;
			}
			else if(strcmp(argv[i], "-json") == 0 && i < argc - 1){
				json_output = argv[i+1];
                ++i;    
			}
            else if(strcmp(argv[i], "-contests") == 0){
                // If this flag is not present, we will assume that 
                // all contests mentioned in the input ballot data will
                // be audited.
                int numc = ToType<int>(argv[i+1]);
                for(int j = 0; j < numc; ++j){
                    // Create new contest structure with given id.
                    int con_id = ToType<int>(argv[i+2+j]);
                    contest_id2index.insert(pair<int,int>(con_id, j));
                    Contest newc;
                    newc.id = con_id;
                    newc.num_rballots = 0;
                    contests.push_back(newc);        
                }
                i += 1 + numc;
            }
		}

        set<string> ballot_ids;
		if(!ReadReportedBallots(rep_blts_file, contests, contest_id2index,
            ballot_ids)){
			cout << "Reported ballots read error. Exiting." << endl;
			return 1;
		}
        int tot_auditable_ballots = ballot_ids.size();       
		
        for(int i = 0; i < contests.size(); ++i){
            Contest &ctest = contests[i];
		    for(int j = 0; j < ctest.rballots.size(); ++j){
                const Ballot &bt = ctest.rballots[j];
                if(bt.prefs.empty())
                    continue;
			    ctest.cands[bt.prefs[0]].ballots.push_back(j);
		    }

            if(simlog)
                cout << "SIMULATING CONTEST " << ctest.id << endl;

		    // Simulate IRV election to determine elimination order/winner 
		    ctest.winner = -1;
		    SimIRV(ctest.rballots, ctest.winner, ctest.cands, ctest.config,
                ctest.outcome, simlog);

		    for(int j = 0; j < ctest.cands.size(); ++j){
			    ctest.cands[j].Reset();
		    }
        }

        Ints successes;
        Ints full_recounts;

        int overall_asn_ballots = -1;
        for(int k = 0; k < contests.size(); ++k){
		    const Contest &ctest = contests[k];
            if(alglog){
                cout << "GENERATING AUDIT FOR CONTEST " << ctest.id << endl;
            }
            mytimespec tstart;
		    GetTime(&tstart);

		    Frontier front;
            int ncands = ctest.cands.size();
		    double lowerbound = -10;

		    // Compute starting facts
		    vector<OrderFacts> facts(ncands);
		    for(int i = 0; i < ncands-1; ++i){
			    const int ci = ctest.outcome[i];
			    for(int j = i+1; j < ncands; ++j){
				    const int cj = ctest.outcome[j];

				    // Can we prove that ci is eliminated before cj. If so,
				    // add to list of ordering facts (with best audit).
				    double asn=EstimateASN_WONLY(ctest.rballots,ctest.cands,
                        rlimit, cj, ci, gamma, lambda, tot_auditable_ballots);

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
			    if(i == ctest.winner)
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
					    cout << ctest.cands[newn.tail[0]].id << " | ";
					    for(int i = 1; i < newn.tail.size(); ++i){
						    cout << ctest.cands[newn.tail[i]].id << " ";
					    }
					    cout << endl;
				    }

				    newn.estimate = FindBestAudit(ctest.cands, ctest.rballots,
                        rlimit, newn.best_audit, facts, newn.tail, gamma, 
                        lambda, alglog, tot_auditable_ballots);

				    if(alglog){
					    if(newn.estimate != -1){
						    cout << "   Best audit ";
						    PrintAudit(newn.best_audit, ctest.cands);
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
			    PrintFrontier(front, ctest.cands);
			    cout << "============================================" << endl;
		    }

		    while(true && !auditfailed){
			    if(lowerbound > 0 && allowed_gap > 0){
				    double max_on_frontier = -1;
				    for(Frontier::const_iterator it = front.begin(); 
                        it != front.end(); ++it){
					    if(it->estimate == -1){
						    max_on_frontier = -1;
						    break;
					    }
					    else{
						    max_on_frontier=max(it->estimate,max_on_frontier);
					    }
				    }
				    if(max_on_frontier != -1 && max_on_frontier - 
                        lowerbound <= allowed_gap){
					    break;
				    }
			    }

			    // Expand node with highest ASN (-1 == infinity)
			    Node toexpand = front.front();
			    if(!toexpand.expandable){
				    break;
			    }

			    front.pop_front();

			    if((toexpand.has_ancestor && 
                    toexpand.best_ancestor.estimate != -1 &&
				    toexpand.best_ancestor.estimate <= lowerbound)){
				    // Replace all descendents of best ancestor with ancestor.
				    ReplaceWithBestAncestor(front, toexpand, ctest.cands, 
                        alglog);
				    continue;
			    }
			    else if(toexpand.estimate != -1 && 
                    toexpand.estimate <= lowerbound){
				    // Don't expand, make "unexpandable",move to back of list.
				    toexpand.expandable = false;
				    front.insert(front.end(), toexpand);
				    continue;
			    } 

			    if(diving){
				    double divelb = PerformDive(toexpand, ctest.cands, 
					    ctest.rballots, rlimit, facts, gamma, lambda,
                        tot_auditable_ballots);
				    if(divelb == -1){
					    // Audit not possible
					    if(alglog){
                            cout << "Diving finds that audit " <<
                                "is not possible." << endl;
                        }
					    auditfailed = true;
					    break;
				    }
				    else if(divelb != -2){
					    if(alglog){ 
						    cout << "Diving LB " << divelb << 
                                " current LB " << lowerbound << endl;
					    }
					    lowerbound = max(lowerbound, divelb);
				    }	

				    if((toexpand.has_ancestor && 
                        toexpand.best_ancestor.estimate != -1 &&
					    toexpand.best_ancestor.estimate <= lowerbound)){
					    // Replace all descendents of best ancestor 
                        // with ancestor.
					    ReplaceWithBestAncestor(front, toexpand, ctest.cands, 
                            alglog);
					    continue;
				    }
				    else if(toexpand.estimate != -1 && 
                        toexpand.estimate <= lowerbound){
					    // Don't expand, just make "unexpandable" and 
                        // move to back of list.
					    toexpand.expandable = false;
					    front.insert(front.end(), toexpand);
					    continue;
				    }
			    }

    			++nodesexpanded;
	    		if(alglog){ 
		    		cout << " Expanding node ";
			    	PrintNode(toexpand, ctest.cands);
				    cout << endl;
			    }
                // For each candidate 'c' not in toexpand.tail, create a new
                // node with node.tail = [c] ++ toexpand.tail
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
					    newn.expandable = (newn.tail.size() == ncands) ? 
                            false : true;
				
					    if(alglog){
						    cout << "TESTING ";
						    cout << ctest.cands[newn.tail[0]].id << " | ";
						    for(int i = 1; i < newn.tail.size(); ++i){
							    cout << ctest.cands[newn.tail[i]].id << " ";
						    }
						    cout << endl;
					    }

					    // Set new nodes best ancestor
					    if(toexpand.has_ancestor){
						    if((toexpand.estimate == -1) ||
							    (toexpand.best_ancestor.estimate != -1 &&
							    toexpand.best_ancestor.estimate <= 
                                toexpand.estimate)){
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
					    newn.estimate = FindBestAudit(ctest.cands, 
                            ctest.rballots, rlimit, newn.best_audit,
                            facts,newn.tail, gamma, lambda, alglog,
                            tot_auditable_ballots);

					    if(!newn.expandable){
						    bool replace = false;
						    if(newn.estimate == -1){
							    if(!newn.has_ancestor || 
                                    newn.best_ancestor.estimate == -1){
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
                                // Replace all descendents of
                                // newn.best_ancestor with the ancestor and
                                // make expandable = false
                                lowerbound = max(lowerbound,
                                    newn.best_ancestor.estimate);
							    ReplaceWithBestAncestor(front, newn, 
                                    ctest.cands, alglog);
						    }
						    else{
							    if(alglog){
								    cout << "   Best audit ";
								    PrintAudit(newn.best_audit, ctest.cands);
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
								    PrintAudit(newn.best_audit, ctest.cands);
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
				    cout << endl << "Size of frontier " << front.size() << 
                        ", Nodes expanded " << nodesexpanded << 
                        ", Current threshold " << lowerbound <<  endl << endl;
			    }
		    }

		    mytimespec tend;
		    GetTime(&tend);

		    double maxasn = -1;
		    if(!auditfailed){
			    // Compile list of audits to complete.
			    Audits audits;
			    for(Frontier::const_iterator it = front.begin(); 
                    it != front.end(); ++it){
				    bool pruneaudit = false;
				    for(Audits::const_iterator at = audits.begin(); 
                        at != audits.end(); ++at){
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
			    for(Audits::const_iterator it = audits.begin(); 
                    it != audits.end();++it){
                    bool subsumed = false;
                    for(Audits::const_iterator jt = audits.begin(); 
                        jt != audits.end(); ++jt){
                        if(jt == it) continue;
                        if(Subsumes(*jt, *it)){
                            subsumed = true;
                            break;
                        }
                    }
                    if(!subsumed){
                        final_config.push_back(*it);
				        PrintAudit(*it, ctest.cands);
				        maxasn = max(maxasn, it->asn);
                    }
			    }
                int in_ballots = maxasn*ctest.config.totalvotes;
			    maxasn *= 100;
                
			    cout << "MAX ASN(%) " << maxasn << endl;
			    cout << "============================================" << endl;

                if(maxasn >= 100){
                    full_recounts.push_back(ctest.id);
                    audits_to_run.push_back(Audits());
                }
                else{
                    audits_to_run.push_back(final_config);
                    successes.push_back(ctest.id);
                    overall_asn_ballots = max(overall_asn_ballots,in_ballots);
                }
		    }
		    else{
			    if(alglog){
				    cout << endl;
				    cout << "AUDIT NOT POSSIBLE" <<endl;
			    }	
                audits_to_run.push_back(Audits());
                full_recounts.push_back(ctest.id);
		    }
		    cout << "TIME," << tend.seconds - tstart.seconds << 
                ",Nodes Expanded," << nodesexpanded
			    << ",MAX ASN(%)," << maxasn << endl;
        }

        cout << "============================================" << endl;
        cout << "SUMMARY" << endl;
        if(successes.size() > 0){
            cout << "Audit found for contests: ";
            for(int i = 0; i < successes.size(); ++i){
                cout << successes[i] << " ";
            }
            cout << "with estimated ballot polls " << overall_asn_ballots
                << endl;
        }
        if(full_recounts.size() > 0){
            cout << "Full recounts required for contests: ";
            for(int i = 0; i < full_recounts.size(); ++i){
                cout << full_recounts[i] << " ";
            }
            cout << endl;
        }
        cout << "============================================" << endl;

        if(json_output != NULL){
            OutputToJSON(contests, audits_to_run, tot_auditable_ballots,
                json_output);
        }
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




