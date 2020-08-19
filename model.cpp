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
            ctest.ndelegates = ToType<int>(columns[2]);
            // Winners start at index 4
            int losers_start = -1;
            for(int i = 4; i < columns.size(); ++i){
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
