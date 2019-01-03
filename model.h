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

struct mytimespec
{
	double seconds;
};

void GetTime(struct mytimespec* t);

void Split2Ints(const std::string &line, const boost::char_separator<char> &sep, Ints &r);
void Split(const std::string &line, const boost::char_separator<char> &sep, std::vector<std::string> &r);

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
	double votes; // Number of votes present with this signature
	Ints prefs;
};

typedef std::vector<Ballot> Ballots;

template <typename T>
T ToType(const std::string &s);

typedef std::map<std::vector<int>,int> I2Map;

void select_random_ballots(int nballots, long seed, Ints &r);

bool ReadReportedBallots(const char *path, Ballots &ballots,
	Candidates &candidates, Config &config, double errorp);

bool ReadActualBallots(const char *path, Ballots &ballots,
	const Candidates &candidates, const Config &config);

#endif
