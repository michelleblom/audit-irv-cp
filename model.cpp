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


#include<fstream>
#include<algorithm>
#include<iostream>
#include<sstream>
#include<time.h>
#include<stdlib.h>
#include<cstdio>
#include<memory>
#include<stdexcept>
#include<string>
#include<array>

#include "model.h"

using namespace std;

typedef boost::char_separator<char> boostcharsep;

void Split2Ints(const string &line, const boostcharsep &sep, Ints &r)
{
	try
	{
		vector<string> result;
		boost::tokenizer<boostcharsep> tokens(line, sep);
		result.assign(tokens.begin(), tokens.end());

		for(int i = 0; i < result.size(); ++i)
		{
			boost::algorithm::trim(result[i]);
			int v  = ToType<int>(result[i]);
			r.push_back(v);
		}
	}
	catch(exception &e)
	{
		stringstream ss;
		ss << e.what();
		throw STVException(ss.str());
	}
	catch(STVException &e)
	{
		throw e;
	}
	catch(...)
	{
		throw STVException("Error in Split Function.");
	}
}


int RandInRange(int a, int b){
	if(a == b) return a;
	if(b + 1 == a){
		if(rand() <= 0.5)
			return a;
		else
			return b;
	}

	return a + rand() % ((b + 1) - a);
}
void TwoRandInRange(int a, int b, int &first, int &second){
	first = RandInRange(a, b);
	if(first == a){
		second = RandInRange(a+1, b);
	}
	else if(first == b){
		second = RandInRange(a, b-1);
	}
	else{
		int s1 = RandInRange(a, first-1);
		int s2 = RandInRange(first+1, b);

		if(RandInRange(0,1) == 0)
			second = s1;
		else
			second = s2;
	}
}

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}

void select_random_ballots(int nballots, long seed, Ints &r){
	stringstream ss;
    // IMPORTANT: replace python3.6 with your python version.
	ss << "python3.6 Sampler.py " << seed << " " << nballots;
	auto x = ss.str();
	cout << ss.str() << endl;
	string res = exec(x.c_str());

	boostcharsep spcom(",");
	Split2Ints(res, spcom, r);
}



void GetTime(struct mytimespec *t)
{
	#ifndef _WIN32
	struct timeval tv;
	gettimeofday(&tv, NULL);
	t->seconds = tv.tv_sec + (tv.tv_usec/1E6);
	#else
	t->seconds = clock()/CLOCKS_PER_SEC;
	#endif
}

template <typename T>
T ToType(const std::string &s) 
{ 
	try
	{
		return boost::lexical_cast<T>(s); 
	}
	catch(...)
	{
		throw STVException("Lexical cast problem");
	}
}

void Split(const string &line, const boostcharsep &sep, vector<string> &r)
{
	try
	{
		boost::tokenizer<boostcharsep> tokens(line, sep);
		r.assign(tokens.begin(), tokens.end());

		for(int i = 0; i < r.size(); ++i)
		{
			boost::algorithm::trim(r[i]);
		}
	}
	catch(exception &e)
	{
		stringstream ss;
		ss << e.what();
		throw STVException(ss.str());
	}
	catch(STVException &e)
	{
		throw e;
	}
	catch(...)
	{
		throw STVException("Error in Split Function.");
	}
}


bool ReadReportedBallots(const char *path, Contests &contests, 
    ID2IX &ct_id2index, set<string> &ballot_ids, const Parameters &params) 
{
	try
	{
		ifstream infile(path);
		boostcharsep spcom(",");

        // The first chunk of lines are related to contest
        // information.
        string line;
        getline(infile, line);
        boost::algorithm::trim(line);
	    int ncontests = ToType<int>(line);

        bool add_all_contests = contests.empty();

        for(int i = 0; i < ncontests; ++i)
        {
            getline(infile, line);

		    vector<string> columns;
		    Split(line, spcom, columns);
             
            int con_id = ToType<int>(columns[1]);
            ID2IX::const_iterator cit = ct_id2index.find(con_id);
            if(add_all_contests && cit == ct_id2index.end()){
                ct_id2index.insert(pair<int,int>(con_id,contests.size()));
                Contest newc;
                newc.id = con_id;
                newc.num_rballots = 0;
                contests.push_back(newc); 
                cit = ct_id2index.find(con_id);
            }

            if(cit == ct_id2index.end())
                continue;

            int con_index = cit->second;
            Contest &ctest = contests[con_index];
            int numcand = ToType<int>(columns[2]);

            for(int j = 0; j < numcand; ++j){
                int can_id = ToType<int>(columns[3+j]);
                Candidate c;
			    c.index = j;
			    c.id = can_id;
			    ctest.cands.push_back(c);
			    ctest.id2index.insert(pair<int,int>(can_id,j));
            }

            ctest.ncandidates = numcand;
        }

        int num_errors = 0;

        // Reading ballots rankings and mapping them to their contest.
        while(getline(infile, line))
        {
			vector<string> columns;
			Split(line, spcom, columns);

            int con_id = ToType<int>(columns[0]);
            ID2IX::const_iterator cit = ct_id2index.find(con_id);
            if(cit == ct_id2index.end())
                continue;

            int con_index = cit->second;
            Contest &ctest = contests[con_index];
 
			Ints prefs;

			Ballot b;
			b.tag = ctest.num_rballots;

            ballot_ids.insert(columns[1]);

			for(int i = 2; i < columns.size(); ++i)
			{
				int ccode = ToType<int>(columns[i]);
				int index = ctest.id2index.find(ccode)->second;
					
				if(find(prefs.begin(),prefs.end(), index) != prefs.end())
				{
					continue;
				}

				b.prefs.push_back(index);
			}

			ctest.aballots.push_back(b);

            if(params.error_prob > 0){
			    double roll = rand() / ((double)RAND_MAX);
			    if(!prefs.empty() && roll <= params.error_prob){
				    roll = rand() / ((double)RAND_MAX);

				    Ints cands_in(ctest.ncandidates, 0);
				    for(int k = 0; k < prefs.size(); ++k){
					    cands_in[prefs[k]] = 1;
				    }
				    Ints cands_out;
				    for(int k = 0; k < ctest.ncandidates; ++k){
					    if(cands_in[k] == 0)
						    cands_out.push_back(k);
					}

					if(b.prefs.size() == 1){
					    int idx = RandInRange(0, cands_out.size()-1);
						if(roll <= 0.50){
						    // Replace candidate
							b.prefs[0] = cands_out[idx];
						}
						else{	
						    // Add random candidate
						    int pos = RandInRange(0, 1);
						    if(pos == 0){
							    b.prefs.insert(b.prefs.begin(),cands_out[idx]);
						    }
						    else{
							    b.prefs.push_back(cands_out[idx]);		
						    }
						}		
					}	
					else{
					    if(roll <= 0.33){
						    // Flip two candidates
						    int f = 0, s = 1;
						    TwoRandInRange(0, b.prefs.size()-1, f, s);
						    int fval = b.prefs[f];
						    b.prefs[f] = b.prefs[s];
						    b.prefs[s] = fval;
					    }
					    else if(roll <= 0.66 && cands_out.size() != 0){
						    // Add random candidate into random position
						    int idx = RandInRange(0, cands_out.size()-1);
						    int pos = RandInRange(0, b.prefs.size()-1);
						    if(pos == 0)
							    b.prefs.insert(b.prefs.begin(),cands_out[idx]);
						    else
							    b.prefs.insert(b.prefs.begin()+pos,cands_out[idx]);
					    }
					    else{
						    // Subtract candidate from random position
						    int pos = RandInRange(0, b.prefs.size()-1);
						    if(pos == 0)
							    b.prefs.erase(b.prefs.begin());
						    else
							    b.prefs.erase(b.prefs.begin()+pos);
					    }
					}
					num_errors += 1;
                }
            }

            ctest.rballots.push_back(b);

            // Note, normally we would ignore ballots with no preferences.
            // However, we may pull them out during sampling when the audit
            // is run, so they should still be counted toward the total
            // number of ballots that are present. 
            if(!b.prefs.empty()){
			    Candidate &cand = ctest.cands[b.prefs.front()];
			    cand.total_votes += 1;
            }
			ctest.num_rballots += 1;
		}

        if(params.error_prob > 0 && params.runlog){
            cout << num_errors << " errors added to reported ballots." << endl;
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
		cout << "Unexpected error reading in reported ballots." << endl;
		return false;
	}

	return true;
}


bool ReadReportedOutcomes(const char *path, Contests &contests, 
    ID2IX &ct_id2index) 
{
	try
	{
		ifstream infile(path);
		boostcharsep spcom(",");

        string line;
        while(getline(infile, line))
        {
			vector<string> columns;
			Split(line, spcom, columns);

            int con_id = ToType<int>(columns[0]);
            ID2IX::const_iterator cit = ct_id2index.find(con_id);
            if(cit == ct_id2index.end()){
                throw STVException("ReadOutcomes: Contest doesn't exist.");
            }

            Contest &ctest = contests[cit->second];
            
            // Winners start at index 2
            int losers_start = -1;
            for(int i = 2; i < columns.size(); ++i){
                if(columns[i] == "losers"){
                    losers_start = i+1;
                    break;
                }

                int winner_id = ToType<int>(columns[i]);
                ID2IX::const_iterator cnd_it = ctest.id2index.find(winner_id);
                if(cnd_it == ctest.id2index.end()){
                    throw STVException("ReadOutcomes: Cand doesn't exist.");
                }

                int winner_idx = cnd_it->second;
                ctest.viable_order.push_back(winner_idx);
                ctest.winners.insert(winner_idx);
            }

            for(int i = losers_start; i < columns.size(); ++i){
                int loser_id = ToType<int>(columns[i]);
                ID2IX::const_iterator cnd_it = ctest.id2index.find(loser_id);
                if(cnd_it == ctest.id2index.end()){
                    throw STVException("ReadOutcomes: Cand doesn't exist.");
                }

                int loser_idx = cnd_it->second;
                ctest.eliminations.push_back(loser_idx);
            }
        }
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
		cout << "Unexpected error reading in reported outcomes." << endl;
		return false;
	}

	return true;
}
