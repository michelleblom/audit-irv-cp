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


#ifndef _MODEL_H
#define _MODEL_H

#include<vector>
#include<set>
#include<exception>
#include<boost/lexical_cast.hpp>
#include<boost/filesystem.hpp>
#include<boost/tokenizer.hpp>
#include<boost/algorithm/string.hpp>
#include<map>

#ifndef _WIN32
#include<sys/time.h>
#endif

typedef std::vector<int> Ints;
typedef std::vector<double> Doubles;
typedef std::vector<Doubles> Doubles2d;
typedef std::set<int> SInts;
typedef std::set<SInts> SInts2d;
typedef std::set<SInts2d> SInts3d;
typedef std::vector<Ints> Ints2d;
typedef std::vector<SInts> VSInts;
typedef std::vector<std::string> Strings;


template <typename T>
T ToType(const std::string &s);


struct mytimespec
{
	double seconds;
};

void GetTime(struct mytimespec* t);

void Split(const std::string &line, const boost::char_separator<char> &sep, 
    std::vector<std::string> &r);

struct Config
{
	int ncandidates;
	double totalvotes;
	std::map<int,int> id2index;

	Config() : ncandidates(0), totalvotes(0) {}
};

class STVException
{
	private:
		const std::string message;

	public:
		STVException(const STVException &me) : message(me.what()){}
		STVException(const std::string &str) : message(str) {}
		const std::string& what() const { return message; }
};


struct Candidate
{
	int id;
	int index;
	double sum_votes;

	Ints ballots;
	Ints ballots_where_appear;

	// For simulation
	int sim_votes;
	Ints sim_ballots;

	int standing;

	Candidate() : id(0), index(0), sum_votes(0), 
		sim_votes(0), standing(1) {}

	void Reset(){
		standing = 1;
		sim_votes = sum_votes;
		sim_ballots.clear();
		sim_ballots.insert(sim_ballots.end(), ballots.begin(), 
			ballots.end());
	}
};

typedef std::vector<Candidate> Candidates;

struct Ballot
{
	int tag;
	// Number of votes present with this signature: Note in this 
    // codebase votes = 1 as we are not grouping ballots together.
	double votes; 
	Ints prefs;
};

typedef std::vector<Ballot> Ballots;

typedef std::map<std::vector<int>,int> I2Map;
typedef std::map<int,int> ID2IX;

struct Contest{
    int id;
    Config config;
    Candidates cands;
    Ballots rballots;

    int num_rballots;

    Ints outcome;
    int winner;

    Contest() : id(0), num_rballots(0), winner(-1) {}
};

typedef std::vector<Contest> Contests;

bool ReadReportedBallots(const char *path, Contests &contests,
	ID2IX &contest_1d2index, std::set<std::string> &ballot_ids);

#endif
