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


bool ReadReportedBallots(const char *path, Ballots &ballots, 
	Candidates &candidates, Config &config, double errorp)
{
	try
	{
		ifstream infile(path);

		// First line is list of candidates.
		boostcharsep spcom(",");
		string line;
		getline(infile, line);

		vector<string> columns;
		Split(line, spcom, columns);

		config.ncandidates = columns.size();
		for(int i = 0; i < columns.size(); ++i)
		{
			int id = ToType<int>(columns[i]);
			Candidate c;
			c.index = i;
			c.id = id;
			candidates.push_back(c);
			config.id2index.insert(pair<int,int>(id,i));
		}

		// Ignore party affiliations
		getline(infile, line);
		// Skip next line (separator)
		getline(infile, line);
		boostcharsep sp(",():");
	
		int num_errors = 0;
	
		int cntr = ballots.size();
		while(getline(infile, line))
		{
			vector<string> columns;
			Split(line, sp, columns);

			int votes = ToType<double>(columns.back());
			Ints prefs;

			Ballot b;
			b.tag = cntr;
			b.votes = ToType<double>(columns.back());

			for(int i = 0; i < columns.size()-1; ++i)
			{
				if(columns[i] == "") continue;
				int ccode = ToType<int>(columns[i]);
				int index = config.id2index.find(ccode)->second;
					
				if(find(prefs.begin(),
					prefs.end(), index) != prefs.end())
				{
					continue;
				}

				prefs.push_back(index);
			}

			if(prefs.empty()) continue;

			for(int j = 0; j < votes; ++j){
				Ballot b;
				b.tag = cntr;
				b.prefs = prefs;
				b.votes = 1;

				double roll = rand() / ((double)RAND_MAX);
				if(roll <= errorp){
					roll = rand() / ((double)RAND_MAX);

					Ints cands_in(config.ncandidates, 0);
					for(int k = 0; k < prefs.size(); ++k){
						cands_in[prefs[k]] = 1;
					}
					Ints cands_out;
					for(int k = 0; k < config.ncandidates; ++k){
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


				ballots.push_back(b);
 
				Candidate &cand = candidates[b.prefs.front()];
				cand.sum_votes += 1;
				config.totalvotes += 1;
				
				++cntr;
			}
		}

		cout << num_errors << " errors added." << endl;
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

// Assumed input format:
// (first_id, second_id, third_id, ...) : #appears
bool ReadActualBallots(const char *path, Ballots &ballots, 
	const Candidates &candidates, const Config &config)
{
	try
	{
		ifstream infile(path);

		// First line is list of candidates.
		boostcharsep spcom(",");
		string line;
		getline(infile, line);

		// Ignore party affiliations
		getline(infile, line);

		// Skip next line (separator)
		getline(infile, line);

		boostcharsep sp(",():");
		
		int cntr = ballots.size();
		while(getline(infile, line))
		{
			vector<string> columns;
			Split(line, sp, columns);

			int votes = ToType<double>(columns.back());
			Ints prefs;

			Ballot b;
			b.tag = cntr;
			b.votes = ToType<double>(columns.back());

			for(int i = 0; i < columns.size()-1; ++i)
			{
				if(columns[i] == "") continue;
				int ccode = ToType<int>(columns[i]);
				int index = config.id2index.find(ccode)->second;
					
				if(find(prefs.begin(),
					prefs.end(), index) != prefs.end())
				{
					continue;
				}

				prefs.push_back(index);
			}

			if(prefs.empty()) continue;

			for(int j = 0; j < votes; ++j){
				Ballot b;
				b.tag = cntr;
				b.prefs = prefs;
				b.votes = 1;
				ballots.push_back(b);
 
				++cntr;
			}
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
		cout << "Unexpected error reading in ballots." << endl;
		return false;
	}

	return true;
}
