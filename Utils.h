#ifndef UTILS_H
#define UTILS_H

#include<iostream>
#include<string>
#include<vector>
#include<iterator>
#include<sstream>
#include<stdexcept>

enum Filetype 
{	
	SAM,
	WIG,
	BED,
	BED5,
	BEDGRAPH,
	ENUM_FILETYPE_SIZE
};

enum Opt
{
	SMOOTH,
	UCSC,
	CLEAN,
	BLACKLIST,
	THREADS,
	ENUM_OPT_SIZE
};

typedef struct 
{
	int sambin;	

	std::string infile;
	std::string outfile;
	Filetype type;

	bool smooth;
	int smooth_win;
	int smooth_shift;
	bool threads;
	int threads_num;

	bool clean;
	std::string cleanfile;

	bool blacklist;
	std::string blacklistfile;

	bool ucsc;
} Opts;

void getInput(int argc, char * argv[], Opts & opt_struct);
void parseOpt(std::string option, Opts & opt_struct);
void parseOpt(std::string option, std::string arg, Opts & opt_struct);
void parseOpt(std::string option, std::string arg1, std::string arg2, Opts & opt_struct);
void initOpts(Opts & opts_struct);
#endif
