#include"Utils.h"
using namespace std;

void getInput(int argc, char * argv[], Opts & opt_struct)
{
	if (argc < 6)
		throw (1);

	string curr = "", prev = string(argv[1]);

	initOpts(opt_struct);

	for (int i = 1; i < argc; i++)
	{
		curr = string(argv[i]);
		
//		cerr << "curr is [" << curr << "] prev is [" << prev << "]" << endl;

		if (curr.find("-") == string::npos) 
		{
			if (i != argc - 1) 
			{
				string next = string (argv[i+1]);
	
				if (next.find("-") != string::npos)
					parseOpt(prev, curr, opt_struct);
				else
				{
					parseOpt(prev, curr, next, opt_struct);
					curr = next;
					i++;
				}
			}
			else
			{
//			cerr << "on incorrect opt" << endl;
				parseOpt(prev, curr, opt_struct);
			}
		}		
		else 
		{
//			cerr << "in else" << endl;
			if (i != argc - 1) 
			{
				string next = string (argv[i+1]);
//				cerr << "next is [" << next << "]" << endl;
				if (next.find("-") == string::npos) 
				{
					prev = curr;
					continue;
				}
			}
//			cerr << "on correct opt" << endl;
			parseOpt(curr, opt_struct);
		}

		prev = curr;
	}	
}

void initOpts(Opts & opt_struct) 
{
	opt_struct.infile = "";
	opt_struct.outfile = "";	
	opt_struct.type = ENUM_FILETYPE_SIZE; // NA essentially
	
	opt_struct.smooth = false;
	opt_struct.smooth_win = 0;
	opt_struct.threads = false;
	opt_struct.threads_num = 0;

	opt_struct.clean = false;
	opt_struct.cleanfile = "";

	opt_struct.blacklist = false;
	opt_struct.blacklistfile = "";

	opt_struct.ucsc = false;
}

void parseOpt(string option, Opts & opt_struct)
{
	if (option.compare("--ucsc") == 0)
		opt_struct.ucsc = true;
	else {
//		cerr << "going to throw\n"<<endl;
		throw invalid_argument("Your opt: [" + option + "] is not a valid opt!\n");
	}
}


void parseOpt(string option, string arg, Opts & opt_struct)
{
	if (option.compare("-i") == 0)
	{
		opt_struct.infile = arg;
	}
	else if (option.compare("-o") == 0)
	{
		opt_struct.outfile = arg;
	}
	else if (option.compare("-c") == 0)
	{
		if (arg.compare("wig") == 0)
			opt_struct.type = WIG;
//		else if (arg.compare("bed") == 0)
//			opt_struct.type = BED;
		else if (arg.compare("bed5") == 0)
			opt_struct.type = BED5;
		else if (arg.compare("bg") == 0)
			opt_struct.type = BEDGRAPH;
		else
			throw invalid_argument("Your outfile type [" + arg + "] is not a valid argument to -c!\n");
	}
	else if (option.compare("--sambin") == 0)
	{
		stringstream ss(arg);
		
		if (!(ss >> opt_struct.sambin))
			throw invalid_argument("Your argument to --sambin [" + arg + "] is not a valid INT!\n");
			
	}
	else if (option.compare("--clean") == 0) 
	{
		opt_struct.clean = true;
		opt_struct.cleanfile = arg;
	}
	else if (option.compare("--blacklist") == 0)
	{
		opt_struct.blacklist = true;
		opt_struct.blacklistfile = arg;
	}
	else if (option.compare("--threads") == 0)
	{
		opt_struct.threads = true;
		stringstream ss(arg);	
	
		if (!(ss >> opt_struct.threads_num))
			throw invalid_argument("Your arg to --threads [" + arg + "] is not a valid INT!\n");
	}
	else
	{
		throw invalid_argument("Could not recognize your option + argument pair: [" + option + "] + [" + arg + "]!\n");
	}
}

void parseOpt(string option, string arg1, string arg2, Opts & opt_struct)
{
	if (option.compare("--smooth") == 0)
	{
		opt_struct.smooth = true;
		stringstream ss1(arg1), ss2(arg2);
		
		if (!((ss1 >> opt_struct.smooth_win) && (ss2 >> opt_struct.smooth_shift)))
			throw invalid_argument("Your args to --smooth [" + arg1 + "] [" + arg2 + "] are not a valid INT!\n");
	}	
	else
	{
		throw invalid_argument("Could not recognize your option + arguments: [" + option + "] + [" + arg1 + "] [" + arg2 + "]!\n");
	}
}
