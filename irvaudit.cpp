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
#include "audit.h"

using namespace std;
using boost::property_tree::ptree;

void print_list(const Ints &list){
    for(int i = 0; i < list.size(); ++i){
        cout << list[i] << " ";
    }
}

struct SimpleNode{
    SInts head;
    Ints tail;
    double estimate;
    AuditSpec best_audit;
};

struct Node{
    SInts head;
    Ints tail;
    double estimate;
    AuditSpec best_audit;

    bool has_ancestor;
    SimpleNode best_ancestor;

    bool expandable;
};

typedef list<Node> Frontier;

SimpleNode CreateSimpleNode(const Node &n){
    SimpleNode sn;
    sn.head = n.head;
    sn.tail = n.tail;
    sn.estimate = n.estimate;
    sn.best_audit = n.best_audit;
    return sn;
}

Node CreateNode(const SimpleNode &sn){
    Node newn;
    newn.head = sn.head;
    newn.tail = sn.tail;
    newn.estimate = sn.estimate;
    newn.best_audit = sn.best_audit;
    newn.has_ancestor = false;
    newn.expandable = false;
    return newn;
}


bool subset_of(const Ints &l1, const Ints &l2){
    if(l1.size() > l2.size())
        return false;

    Ints l1_temp(l1);
    Ints l2_temp(l2);

    sort(l1_temp.begin(), l1_temp.end());
    sort(l2_temp.begin(), l2_temp.end());

    for(int i = 0; i < l1_temp.size(); ++i)
        if(l1_temp[i] != l2_temp[i])
            return false;

    return true;
}


// Does audit a1 subsume a2
bool Subsumes(const AuditSpec &a1, const AuditSpec &a2){
    if(a1.type == VIABLE && a2.type == VIABLE && a1.winner == a2.winner){
        // If a1's eliminated set is a subset of a2's, return true
        if(subset_of(a1.eliminated, a2.eliminated)){
            return true;
        } 
    }
    else if(a1.type==NONVIABLE && a2.type==NONVIABLE && a1.winner==a2.winner){
        // if a2's eliminated set is a subset of a1's, return true
        if(subset_of(a2.eliminated, a1.eliminated)){
            return true;
        }
    }

    return false;    
}

bool AreAuditsEqual(const AuditSpec &a1, const AuditSpec &a2)
{
    if(a1.type != a2.type){
        return false;
    }

    if(a1.winner != a2.winner || a1.loser != a2.loser)
        return false;

    if(a1.eliminated != a2.eliminated)
        return false;

    return true;
}

bool DescendantOf(const Node &d, const Node &a){
    if(d.head != a.head){
        return false;
    }

    if(d.tail.size() <= a.tail.size()){
        return false;
    }

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
    if(n.tail.size() > 0){
        cout << cand[n.tail[0]].id << " | ";
        for(int i = 1; i < n.tail.size(); ++i){
            cout << cand[n.tail[i]].id << " ";
        }
    }
    cout << "( ";
    for(SInts::const_iterator cit=n.head.begin(); cit!=n.head.end(); ++cit){
        cout << cand[*cit].id << " ";
    }
    cout << ") [";
    cout << ((n.estimate == -1) ? -1 : n.estimate) << "] ";

    if(n.has_ancestor){
        cout << " (Best Ancestor ";
        if(n.best_ancestor.tail.size() > 0){
            cout << cand[n.best_ancestor.tail[0]].id << " | ";
            for(int i = 1; i < n.best_ancestor.tail.size(); ++i){
                cout << cand[n.best_ancestor.tail[i]].id << " ";
            }
        }
        cout << "( ";
        for(SInts::const_iterator cit = n.best_ancestor.head.begin(); 
            cit != n.best_ancestor.head.end(); ++cit){
            cout << cand[*cit].id << " ";
        }
        cout << ")";

        cout << " [" << ((n.best_ancestor.estimate == -1) ? 
            -1 : n.best_ancestor.estimate) << "])";
    }
}


void OutputToJSON(const Contests &contests, const vector<Audits> &torun,
    const Parameters &params, const char *json_file)
{
    try{
        ptree pt;
        ptree children;

        double overall_maxasn = -1;
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

                child.put("winner", ctest.cands[spec.winner].id);
                if(spec.loser != -1)
                    child.put("loser", ctest.cands[spec.loser].id);
                else
                    child.put("loser", -1);

                ptree aelim;
                for(int j = 0; j < spec.eliminated.size(); ++j){
                    ptree c;
                    c.put("", ctest.cands[spec.eliminated[j]].id);
                    aelim.push_back(std::make_pair("", c));
                }
                child.add_child("already_eliminated", aelim);

                if(spec.type == IRV)
                    child.put("assertion_type", "IRV_ELIMINATION");
                else if(spec.type == VIABLE)
                    child.put("assertion_type", "VIABLE");
                else if(spec.type == NONVIABLE)
                    child.put("assertion_type", "NONVIABLE");
                else 
                    child.put("assertion_type", "NEB");

                caudit_children.push_back(std::make_pair("", child));
                maxasn = max(maxasn, spec.asn);
            }

            caudit.put("contest", ctest.id);
            int maxasn_pc=ceil(100*(maxasn/params.tot_auditable_ballots));

            caudit.put("Expected Polls (#)", maxasn);
            caudit.put("Expected Polls (%)", maxasn_pc);

            caudit.add_child("assertions", caudit_children);
            overall_maxasn = max(overall_maxasn, maxasn);
            
            children.push_back(std::make_pair("", caudit));
        }
        if(overall_maxasn != -1){
            pt.put("Overall Expected Polls (#)", overall_maxasn);
            pt.put("Ballots involved in audit (#)", 
                params.tot_auditable_ballots);

            ptree parameters;
            parameters.put("risk_limit", params.risk_limit);
            pt.add_child("parameters", parameters);

            pt.add_child("audits", children);
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
    if(audit.type == VIABLE)
        cout << "V," << cand[audit.winner].id << ",Eliminated";
    else if(audit.type == NONVIABLE)
        cout << "NV," << cand[audit.winner].id << ",Eliminated";
    else if(audit.type == IRV){
        cout << "IRV," << cand[audit.winner].id << "," << 
            cand[audit.loser].id << ",Eliminated";
    }
    else{
        cout << "NEB," << cand[audit.winner].id << "," << 
            cand[audit.loser].id << ",Eliminated";
    }

    for(int i = 0; i < audit.eliminated.size(); ++i){
        cout << "," << cand[audit.eliminated[i]].id;
    }
    cout << endl;
}

void PrintFrontier(const Frontier &front, const Candidates &cand){
    for(Frontier::const_iterator it = front.begin(); it != front.end(); ++it){
        cout << "> ";
        PrintNode(*it, cand);
        cout << endl;
    }
}

int ComputeTallies(const Contest &ctest,const Ints &eliminated,Ints &tallies){
    Ints elim(ctest.ncandidates, 0);
    int exhausted = 0;
    for(int i = 0; i < eliminated.size(); ++i){
        elim[eliminated[i]] = 1;
    }

    for(int i = 0; i < ctest.rballots.size(); ++i){
        const Ints &prefs = ctest.rballots[i].prefs;
        bool ex = true;
        for(int j = 0; j < prefs.size(); ++j){
            int pc = prefs[j];
            if(elim[pc])
                continue;

            tallies[pc] += 1;
            ex = false;
            break;
        }
        if(ex) exhausted += 1;
    } 
    return exhausted;
}

int ComputeNEBTally(const Contest &ctest, int loser, int winner){

    int loser_tally = 0;

    for(int i = 0; i < ctest.rballots.size(); ++i){
        const Ints &prefs = ctest.rballots[i].prefs;

        // If loser appears before winner, then increment loser_tally
        for(int j = 0; j < prefs.size(); ++j){
            if(prefs[j] == loser){
                loser_tally += 1;
                break;
            }

            if(prefs[j] == winner){
                break;
            }
        }
    } 
    return loser_tally;
}

double FindBestAudit(const Contest &ctest, const Parameters &params,
    Node &node, const map<int,AuditSpec> &initial_viables,
    const Ints &has_init_viable, const Audits2d &nebs, 
    const Bools2d &has_neb, bool alglog) 
{
    double best_estimate = -1;

    // -------------------------------------------------------------------
    // Possible audits: 
    // -------------------------------------------------------------------
    // -- one of the winners is not viable if we treat everyone outside of
    //    the winners set as eliminated.
    // 
    // -- one of the candidates not in the winners set is viable given no
    //    one has been eliminated. 
    //
    // -- node.tail[0] is viable given all non-mentioned candidates
    //    have been eliminated.
    //
    // -- IRV assertions: node.tail[0] beats a candidate still standing
    //    at that stage. 
    //  
    // -- NEB(c1, c2) c1 cannot be eliminated before c2 as firstpref(c1)
    //    is greater than all_mentions_before_c1(c2)
    // -------------------------------------------------------------------
    Ints eliminated;
    Ints unmentioned;
    Ints tallies1(ctest.ncandidates, 0);
    Ints tallies2(ctest.ncandidates, 0);

    bool empty = node.tail.empty();

    // Checking: one of the candidates not in the winners set is viable 
    //    given no one has been eliminated. V(c, emptyset)
    // Also: define eliminated & unmentioned sets
    for(int i = 0; i < ctest.ncandidates; ++i){
        if(node.head.find(i) != node.head.end())
            continue;

        if(find(node.tail.begin(),node.tail.end(), i) == node.tail.end())
            unmentioned.push_back(i);

        eliminated.push_back(i);
        if(empty && has_init_viable[i] == 1){
            const AuditSpec &as = initial_viables.find(i)->second;
            if(best_estimate == -1 || as.asn < best_estimate){
                best_estimate = as.asn;
                node.best_audit = as;
            }
        }
    }

    // Compute tallies of candidates assuming candidates c with 
    // eliminated[c] = 1 are eliminated (tallies 1) or candidate c with
    // unmentioned[c] = 1 are eliminated (tallies 2).
    if(empty){
        int ex1 = ComputeTallies(ctest, eliminated, tallies1);

        // Checking: one of the winners is not viable if we treat everyone 
        //    outside of the winners set as eliminated. NV(c, C \setminus V)
        for(SInts::iterator cit = node.head.begin();
            cit != node.head.end(); ++cit)
        {
            double asn = EstimateASN_NONVIABLE(ctest,*cit,tallies1,ex1,params);
            if((best_estimate == -1 && asn != -1) || (asn != -1 &&
                asn < best_estimate)){
                best_estimate = asn;
                node.best_audit.asn = asn;
                node.best_audit.type = NONVIABLE;
                node.best_audit.winner = *cit;
                node.best_audit.loser = -1;

                node.best_audit.eliminated = eliminated;
            }
        }

        // Checking: that whether one of the unmentioned candidates, not in 
        // the viable set, could *not* have been eliminated before one of
        // the reportedly viable candidates. 
        for(int i = 0; i < unmentioned.size(); ++i){
            for(SInts::iterator cit = node.head.begin();
                cit != node.head.end(); ++cit)
            {
                if(has_neb[i][*cit]){
                    const AuditSpec &neb_icit = nebs[i][*cit];
                    if((best_estimate == -1 && neb_icit.asn != -1) ||
                        (neb_icit.asn != -1 && neb_icit.asn < best_estimate))
                    {
                        best_estimate = neb_icit.asn;
                        node.best_audit = neb_icit;
                    }
                }
            }
        }
    }

    // Checking: node.tail[0] is viable given all non-mentioned candidates
    //    have been eliminated. V(c, unmentioned \setminus {c})
    if(!empty){
        ComputeTallies(ctest, unmentioned, tallies2);
        double asn = EstimateASN_VIABLE(ctest, node.tail[0], tallies2, params);

        if((best_estimate == -1 && asn != -1) || (asn != -1 && 
            asn < best_estimate)){
            best_estimate = asn;
            node.best_audit.asn = asn;
            node.best_audit.type = VIABLE;
            node.best_audit.winner = node.tail[0];
            node.best_audit.loser = -1;
            node.best_audit.eliminated = unmentioned;
        } 

        // Checking: IRV assertions! 
        AuditSpec bia;
        bia.type = IRV;
        bia.eliminated = unmentioned;
        bia.winner = node.tail[0];

        double bia_asn = FindBestIRV_NEB(ctest, node.tail, node.head, 
            params, tallies2, nebs, has_neb, bia);

        if((best_estimate == -1 && bia_asn != -1) || (bia_asn != -1 && 
            bia_asn < best_estimate)){
            best_estimate = bia_asn;
            node.best_audit = bia;
        }
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

double PerformDive(const Node &toexpand, const Contest &ctest, 
    const map<int,AuditSpec> &initial_viables, const Ints &has_init_viable,
    const Audits2d &nebs, const Bools2d &has_neb, const Parameters &params)
{
    for(int i = 0; i < ctest.ncandidates; ++i){
        if(find(toexpand.tail.begin(), toexpand.tail.end(), i) ==
            toexpand.tail.end() && toexpand.head.find(i) == 
            toexpand.head.end()){
                    
            Node newn;
            newn.tail.push_back(i);
            for(int j = 0; j < toexpand.tail.size(); ++j){
                newn.tail.push_back(toexpand.tail[j]);
            }

            newn.estimate = -1;
            newn.best_audit.asn = -1;
            newn.estimate = -1;
            newn.expandable = (newn.tail.size()+newn.head.size() 
                == ctest.ncandidates) ? false : true;

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
            newn.estimate = FindBestAudit(ctest, params, newn,
                initial_viables, has_init_viable, nebs, has_neb, false);

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
                return PerformDive(newn, ctest, initial_viables, 
                    has_init_viable, nebs, has_neb, params); 
            }
            break;
        }
    }
    return -2; 
}

/**
 * Summary of command line options:
 *
 * -rep_ballots FILE     File containing reported ballots (electronic records)
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
 *                           those contests. IF NOT PRESENT, ASSUME ALL
 *                           CONTESTS MENTIONED IN INPUT WILL BE AUDITED.  
 * 
 * -help                 Print usage instructions.     
 * */

void PrintUsageInstructions(){
// TODO
}

int main(int argc, const char * argv[]) 
{
    try
    {
        Contests contests;

        bool alglog = false;

        Parameters params;
        params.risk_limit = 0.05;
        params.tot_auditable_ballots = 0;

        bool diving = true;

        const char *rep_blts_file = NULL;
        const char *rep_outc_file = NULL;
        const char *json_output = NULL;

        double allowed_gap = 0;
        double threshold_pc = 0.15;

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
            else if(strcmp(argv[i], "-rep_outcome") == 0 && i < argc-1){
                rep_outc_file = argv[i+1];
                ++i;
            }
            else if(strcmp(argv[i], "-agap") == 0 && i < argc-1){
                allowed_gap = atof(argv[i+1]);
                ++i;
            }
            else if(strcmp(argv[i], "-threshold_pc") == 0 && i < argc-1){
                threshold_pc = atof(argv[i+1]);
                ++i;
            }
            else if(strcmp(argv[i], "-r") == 0 && i < argc-1){
                params.risk_limit = atof(argv[i+1]);
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

        if(rep_blts_file == NULL || rep_outc_file == NULL){
            cout << "Reported ballots or outcome not provided." << endl;
            return 1;
        }

        set<string> ballot_ids;
        if(!ReadReportedBallots(rep_blts_file, contests, contest_id2index,
            ballot_ids, params)){
            cout << "Reported ballots read error. Exiting." << endl;
            return 1;
        }

        if(!ReadReportedOutcomes(rep_outc_file, contests, contest_id2index)){
            cout << "Reported outcomes read error. Exiting." << endl;
            return 1;
        }

        params.tot_auditable_ballots = ballot_ids.size();

        // Express allowed gap in ballots rather than as a fraction of ballots
        allowed_gap *= params.tot_auditable_ballots;       
        
        for(int i = 0; i < contests.size(); ++i){
            Contest &ctest = contests[i];
            for(int j = 0; j < ctest.rballots.size(); ++j){
                const Ballot &bt = ctest.rballots[j];
                if(bt.prefs.empty())
                    continue;
                ctest.cands[bt.prefs[0]].ballots.push_back(j);
            }
            ctest.threshold = floor(threshold_pc*ctest.rballots.size() + 1);
            ctest.threshold_fr = threshold_pc;
            if(alglog){
                cout << "Threshold (contest " << ctest.id << "): " 
                    << ctest.threshold << " ballots" << endl;
            }
        }

        Ints successes;
        Ints full_recounts;

        // NOTE: asn's are defined in ballots, not proportions/percentages
        double overall_asn_ballots = -1;
        for(int k = 0; k < contests.size(); ++k){
            const Contest &ctest = contests[k];
            if(alglog){
                cout << "GENERATING AUDIT FOR CONTEST " << ctest.id << endl;
            }
            mytimespec tstart;
            GetTime(&tstart);

            // List of audits to complete.
            Audits audits;
            double lowerbound = -10;
            bool auditfailed = false;

            // We first need assertions that will check the viability of 
            // the "reportedly viable" candidates. This can either be an
            // assertion that says "candidate c is viable even when no one 
            // else is eliminated" or "candidate c is viable when all 
            // reportedly non viable candidates have been eliminated"

            // Identify if there is a V-{} check that is feasible for
            // each reportedly viable candidate. This check is an assertion
            // that says "candidate i is viable even when no other cands
            // have been eliminated".
            map<int,AuditSpec> initial_viables;
            Ints has_init_viable(ctest.ncandidates, 0);
  
            for(SInts::const_iterator cit = ctest.winners.begin(); 
                cit != ctest.winners.end(); ++cit){
                // Can we assert that candidate *cit is viable given that 
                // no other candidates Ints() have been eliminated?
                Ints tallies1(ctest.ncandidates, 0);
                ComputeTallies(ctest, Ints(), tallies1);
                Ints tallies2(ctest.ncandidates, 0);
                ComputeTallies(ctest, ctest.eliminations, tallies2);
 
                double asn1 = EstimateASN_VIABLE(ctest,*cit,tallies1,params);
                double asn2 = EstimateASN_VIABLE(ctest,*cit,tallies2,params);

                if(asn1 == -1 && asn2 == -1){
                    cout << "Audit for contest " << ctest.id << " is not "
                        "possible, at least one reportedly viable candidate "
                        << "only *just* meets threshold." << endl;
                    auditfailed = true;
                    break;
                }

                if(asn1 != -1){
                    AuditSpec aspec;
                    aspec.type = VIABLE;
                    aspec.winner = *cit;
                    aspec.loser = -1;
                    aspec.eliminated.clear(); 
                    aspec.asn = asn1;

                    initial_viables[*cit] = aspec;
                    has_init_viable[*cit] = 1;
                    if(asn2 == -1 || asn1 <= asn2){
                        audits.push_back(aspec);
                        lowerbound = max(lowerbound, asn1);
                        if(alglog){
                            cout << "Added audit: ";
                            PrintAudit(aspec, ctest.cands);
                        }
                    }
                }

                if(asn2 != -1 && (asn1 == -1 || asn2 < asn1)){
                    AuditSpec aspec;
                    aspec.type = VIABLE;
                    aspec.winner = *cit;
                    aspec.loser = -1;
                    aspec.eliminated = ctest.eliminations;
                    aspec.asn = asn2;

                    audits.push_back(aspec);
                    lowerbound = max(lowerbound, asn2);

                    if(alglog){
                        cout << "Added audit: ";
                        PrintAudit(aspec, ctest.cands);
                    }
                }
            }    
            if(auditfailed){
                audits_to_run.push_back(Audits());
                full_recounts.push_back(ctest.id);
                continue;
            }

            // Create a matrix of NEB assertions that could be used to rule
            // out an outcome.
            if(alglog){
                cout << "Finding NEB assertions" << endl;
            }
            Bools2d has_neb;
            Audits2d nebs;
            for(int i = 0; i < ctest.ncandidates; ++i){
                Bools has_neb_i(ctest.ncandidates, false);
                Audits nebs_i(ctest.ncandidates, AuditSpec());

                // i's first preferences is equal to cand.total_votes
                const Candidate &c1 = ctest.cands[i];
                for(int j = 0; j < ctest.ncandidates; ++j){
                    if(i == j) continue;

                    // Compute j's tally including all ballots that 
                    // preference j before i
                    int c2_neb = ComputeNEBTally(ctest, j, i);
                    int neither = params.tot_auditable_ballots - c2_neb -
                        c1.total_votes;

                    if(c1.total_votes > c2_neb){
                        // Margin is 2*assorter_mean - 1
                        // assorter mean is average of (winner-loser+1)/2
                        // for each CVR where winner = 1 if vote is for 'i',
                        // loser = 1 if vote is for 'j', and result is 0.5
                        // if the vote is for neither.
                        double amean = (c1.total_votes + 0.5*neither)/
                            params.tot_auditable_ballots;

                        double margin = 2*amean - 1;
                        int ssize = estimate_sample_size(margin,
                            params.tot_auditable_ballots,params.risk_limit,0);
                        cout << i << "," << j <<"," << c1.total_votes <<
                            " " << c2_neb << " " << neither << ", "
                            << margin  << ", " << ssize << endl;
                        if(ssize < params.tot_auditable_ballots && 
                            ssize != -1){
                            AuditSpec &spec = nebs_i[j];
                            spec.winner = i;
                            spec.loser = j;
                            spec.type = NEB;
                            spec.asn = ssize;

                            has_neb_i[j] = true;
                        }
                    }
                }
                has_neb.push_back(has_neb_i);
                nebs.push_back(nebs_i);
            }

            if(alglog){
                cout << "Starting lower bound on ASN: " <<
                    lowerbound << " ballots (" <<
                    100*(lowerbound/params.tot_auditable_ballots)
                     << "%)" << endl;
            }

            Frontier front;

            // Build initial frontier by forming all 2^n (where n is the
            // number of candidates) subsets of candidates to represent 
            // possible "viable sets".
            if(alglog){
                cout << "Constructing initial frontier" << endl;
            }    

            const int NSETS = pow(2, ctest.ncandidates);

            if(alglog) cout<<NSETS<<" nodes to be added to frontier"<<endl;

            for(int i = 0; i < NSETS; i++){
                Node newn;
                for(int j = 0; j < ctest.ncandidates; j++){
                    if(i & (1 << j)){
                        newn.head.insert(j);
                    }
                }
                if(newn.head.empty())
                    continue;

                if(newn.head == ctest.winners)
                    continue;

                newn.best_audit.asn = -1;
                newn.estimate = -1;
                newn.expandable = true;
                newn.has_ancestor = false;

                // Find best audit to rule out outcome 
                mytimespec t1;
                GetTime(&t1);
  
                newn.estimate = FindBestAudit(ctest, params, newn, 
                    initial_viables, has_init_viable, nebs, has_neb, alglog);

                mytimespec t2;
                GetTime(&t2);

                if(alglog){
                    cout << "Added node with estimate " << newn.estimate
                        << ", " << t2.seconds - t1.seconds << "s" << endl;
                }
                InsertNode(front, newn);
            }

            int nodesexpanded = 0;

            if(alglog){
                cout << "========================================" << endl;
                cout << "Initial Frontier:" << endl;
                PrintFrontier(front, ctest.cands);
                cout << "========================================" << endl;
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
                            max_on_frontier = max(it->estimate,
                                max_on_frontier);
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
                    // Replace descendents of best ancestor with ancestor.
                    ReplaceWithBestAncestor(front, toexpand, 
                        ctest.cands, alglog);
                    continue;
                }
                else if(toexpand.estimate != -1 && 
                    toexpand.estimate <= lowerbound){
                    // Don't expand, make "unexpandable", move to back.
                    toexpand.expandable = false;
                    front.insert(front.end(), toexpand);
                    continue;
                }    

                if(diving){
                    double divelb = PerformDive(toexpand, ctest, 
                        initial_viables,has_init_viable,nebs,has_neb,params);
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
                        ReplaceWithBestAncestor(front, toexpand, 
                            ctest.cands, alglog);
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

                // For each candidate 'c' not in toexpand.tail or 
                // toexpand.head, create a new node with 
                // node.tail = [c] ++ toexpand.tail
                for(int i = 0; i < ctest.ncandidates; ++i){
                    if(find(toexpand.tail.begin(),toexpand.tail.end(),i) ==
                        toexpand.tail.end() && toexpand.head.find(i) == 
                        toexpand.head.end()){
                    
                        Node newn;
                        newn.head = toexpand.head;
                        newn.tail.push_back(i);
                        for(int j = 0; j < toexpand.tail.size(); ++j){
                            newn.tail.push_back(toexpand.tail[j]);
                        }

                        newn.estimate = -1;
                        newn.best_audit.asn = -1;
                        newn.expandable=(newn.tail.size()+newn.head.size() 
                            == ctest.ncandidates) ? false : true;
                
                        if(alglog){
                            cout << "TESTING ";
                            cout << ctest.cands[newn.tail[0]].id << " | ";
                            for(int i = 1; i < newn.tail.size(); ++i){
                                cout << ctest.cands[newn.tail[i]].id << " ";
                            }
                            cout << "( ";
                            for(SInts::const_iterator cit=newn.head.begin();
                                cit != newn.head.end(); ++cit){
                                cout << ctest.cands[*cit].id << " ";
                            }
                            cout << ")" << endl;
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
                                newn.best_ancestor = 
                                    CreateSimpleNode(toexpand); 
                            }
                        } 
                        else{
                            newn.best_ancestor = CreateSimpleNode(toexpand); 
                        }

                        newn.has_ancestor = true;
                        newn.estimate = FindBestAudit(ctest, params, newn, 
                            initial_viables, has_init_viable, nebs, 
                            has_neb, alglog);

                        if(alglog){
                            cout << newn.estimate << endl;
                        }

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
                                    PrintAudit(newn.best_audit,ctest.cands);
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
                                    PrintAudit(newn.best_audit,ctest.cands);
                                    cout << endl;
                                }
                                else{
                                    cout <<"   Cannot be disproved."<< endl;
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
                        ", Current threshold " << lowerbound << 
                        " ballots (" << 100*(lowerbound/
                        params.tot_auditable_ballots) << "%)" <<endl<<endl;
                }
            }
            
            mytimespec tend;
            GetTime(&tend);

            double maxasn = -1;
            if(!auditfailed){
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

                cout << "=========================================" << endl;
                cout << "AUDITS REQUIRED" << endl;
                maxasn = 0;
                Audits final_config;

                // Sort audits from largest to smallest ASN
                sort(audits.begin(), audits.end(), RevCompareAudit);

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
                double in_pc = (maxasn/params.tot_auditable_ballots)*100;
               
                cout << final_config.size() << " assertions" << endl; 
                cout << "MAX ASN(%) " << in_pc << endl;
                cout << "=========================================" << endl;

                if(maxasn >= params.tot_auditable_ballots){
                    full_recounts.push_back(ctest.id);
                    audits_to_run.push_back(Audits());
                }
                else{
                    audits_to_run.push_back(final_config);
                    successes.push_back(ctest.id);
                    overall_asn_ballots = max(overall_asn_ballots,maxasn);
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
            double in_pc = (maxasn/params.tot_auditable_ballots)*100;
            cout << "TIME," << tend.seconds - tstart.seconds << 
                ",Nodes Expanded," << nodesexpanded
                << ",MAX ASN(%)," << in_pc << endl;
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
            OutputToJSON(contests, audits_to_run, params, json_output);
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




