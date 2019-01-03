#ifndef _AUDIT_H
#define _AUDIT_H

#include "model.h"

struct AuditSpec{
	double asn;
	int winner;
	int loser;

	Ints eliminated;
	bool wonly;
};

struct Result{
	int polls;
	int remaining_hypotheses;
};

typedef std::vector<AuditSpec> Audits;

double EstimateASN_WONLY(const Ballots &rep_ballots, 
	const Candidates &cand, double rlimit, int winner, 
	int loser, double gamma, double lambda);

double EstimateSampleSize(const Ballots &rep_ballots, const Candidates &cand, 
	double rlimit, const Ints &tail, AuditSpec &best_audit,
	double gamma, double lambda);

Result RunSingleWinnerLoserAudit(const Ballots &act_ballots, 
	const Ballots &rep_ballots, const Candidates &cand,
	double rlimit, int winner, int loser, const Ints &plist,
	double gamma, double lambda);

Result RunAudit(const Ballots &act_ballots, const Ballots &rep_ballots,
	const Candidates &cand, double rlimit, const Ints &winners, 
	const Ints &losers,const Ints &plist, double gamma, double lambda);  

bool LoadAudits(const char *path, Audits &audits,
	const Candidates &candidates, const Config &config);

#endif
