#include"Utils.h"
#include"Files.h"
using namespace std;

int main(int argc, char * argv[]) 
{
	try
	{
		Opts opt_struct;
		getInput(argc, argv, opt_struct);


		FileInit init(opt_struct.infile, opt_struct.outfile, opt_struct.type, opt_struct.sambin, opt_struct.threads_num);

		File * myFile = init.getFileObj();

		myFile->setOpts(opt_struct.smooth, opt_struct.smooth_win, opt_struct.smooth_shift, opt_struct.blacklist, opt_struct.blacklistfile, opt_struct.clean, opt_struct.cleanfile);

		myFile->output(opt_struct.ucsc);

	/*		
		File file1;
		Wig wig;
	
		File * fp1 = &wig;
		File& fp2 = wig;

		file1.output(); // File output
		fp1->output();	// Wig output
		fp2.output();	// Wig output
	*/
	}
	catch(int e) 
	{
		if (e == 1) 
		{
			string usagestr = "usage: convert-everything [optional] -i <STRING: input filename> -c <STRING: filetype to convert TO> -o <STRING: output filename>\n";
			usagestr += "\nIF INPUT FILETYPE IS SAM\n";
			usagestr += "THEN --sambin is REQUIRED\n";
			usagestr += "\n--sambin <INT: bin size> : bin size for sam file conversions.\n";
			usagestr += "\n-c options:\n";
			usagestr += "\twig: variableStep wig file\n";
//			usagestr += "\tbed: 3 column bed file\n";
			usagestr += "\tbed5: 5 column bed file (includes score and name (\".\"))\n";
			usagestr += "\tbg: 4 column bedgraph format\n";
			usagestr += "\n[optional]:\n";
			usagestr += "\t--clean <STRING: chromosome length file>: remove track & header lines from input, ensure all peaks are within chromosomes\n";
			usagestr += "\t--ucsc : print with UCSC track line\n";
			usagestr += "\t--smooth <INT: bp window to smooth> <INT: bp to shift> : Smooth entire dataset into bins of size given\n";
			usagestr += "\t--blacklist <STRING: blacklist filename> : Removes regions matching those in the file provided\n";
//			usagestr += "\t--threads <INT: number of threads> : max number of threads to run\n";
			// not using threads
			cerr << usagestr;
		}

		return 1;
	}

	catch(const invalid_argument &e)
	{
		cerr << "INVALID_ARGUMENT_ERROR: "  << e.what() << endl;
		return 1;		
	}

}
