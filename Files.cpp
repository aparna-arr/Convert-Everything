#include"Files.h"
using namespace std;

bool peakCmp(Peak peak1, Peak peak2)
{
	return (peak1.start < peak2.start);
}

FileInit::FileInit(string infilename, string outfilename, Filetype type, int sambinsize, int threads)
{
	outfile = outfilename;
	outfiletype = type;
	sambin = sambinsize;

	peaks = new unordered_map<string, vector<Peak>>();	

	readIn(infilename);
}

void FileInit::readIn(std::string infilename)
{
	// figure out what kind of infile this is
	// call appropriate specific read in 

	ifstream infile(infilename);

	if (!infile.is_open())
		throw(1); // throw error

	string line;

	bool detectFileType = false;
	Filetype type = ENUM_FILETYPE_SIZE; // init

	while(!detectFileType && getline(infile, line) && !infile.bad())
	{
		stringstream ss(line);
		string test;
		ss >> test;
		
		if (test.compare("variableStep") == 0)
		{
			type = WIG;
			detectFileType = true;
		}
		else if(test.find("@") != string::npos)
		{
			type = SAM;
			detectFileType = true;
		}
		else if (test.find("chr") != string::npos)
		{
			// some kind of bed format
			int start, end;	
			double value;
			string chr, name;
			// test if we have a value element, discard other columns
			if (ss >> start >> end >> value)
			{
				detectFileType = true;
				type = BEDGRAPH;
			}
			else 
			{
				stringstream ss2(line);
	
				if (ss2 >> chr >> start >> end >> name >> value)
				{
					detectFileType = true;
					type = BED;
				}
				// else throw error
			}
		}
	}

	if (detectFileType)
	{
		// rewind the file
		infile.clear();
		infile.seekg(0, infile.beg);
		
		switch(type) 
		{
			case SAM: readInSam(infile);
			break;
			case WIG: readInWig(infile);
			break;
			case BED: readInBed(infile);
			break;
			case BEDGRAPH: readInBedGraph(infile);
			break;
			default: throw(1); // throw error
		}
	}
	// else throw error
	
	infile.close();
}

void FileInit::readInWig(ifstream& fp)
{
	cerr << "readInWig()" << endl;
	string line;

	string curr_chr = "INIT";
	int span;

	while(getline(fp, line) && !fp.bad())
	{
		stringstream ss(line);

		int pos;
		double value;

		if (line.find("variableStep") == 0)
		{
			string chr, str;
		
			ss >> str; // variableStep 
			ss >> str; // chrom

			size_t found = str.find("=");
				
			chr = str.substr(found+1);

			ss >> str; // span
			found = str.find("=");
			stringstream spanstr(str.substr(found+1));
			
			spanstr >> span;

			if (chr != curr_chr)
				curr_chr = chr;
		}
		else if (ss >> pos >> value)
		{
			Peak tmp;
			tmp.start = pos;
			tmp.end = pos + span;
			tmp.value = value;
	
			(*peaks)[curr_chr].push_back(tmp);
		}
		else 
		{
			cerr << "skipping line [" << line << "]" << endl;
		} 
	}

	if ((*peaks).empty())	
		throw (1); // throw error
}

void FileInit::readInBed(ifstream& fp)
{
	cerr << "readInBed()" << endl;
	string line;

	while(getline(fp, line) && !fp.bad())
	{
		string chr, name;
		int start, end;
		double value;

		stringstream ss(line);

		if (!(ss >> chr >> start >> end >> name >> value))	
		{
			cerr << "skipping line [" << line << "]" << endl;
			continue;
		}
		
		Peak curr;
		curr.start = start;
		curr.end = end;
		curr.value = value;
	
		(*peaks)[chr].push_back(curr);	
	}

	if ((*peaks).empty())	
		throw (1); // throw error

}
void FileInit::readInBedGraph(ifstream& fp)
{
	cerr << "readInBedGraph()" << endl;
	string line;

	while(getline(fp, line) && !fp.bad())
	{
		string chr;
		int start, end;
		double value;

		stringstream ss(line);

		if (!(ss >> chr >> start >> end >> value))	
		{
			cerr << "skipping line [" << line << "]" << endl;
			continue;
		}
	
		Peak curr;
		curr.start = start;
		curr.end = end;
		curr.value = value;
		
		(*peaks)[chr].push_back(curr);	
	}

	if ((*peaks).empty())	
		throw (1); // throw error
}

void FileInit::readInSam(ifstream& fp)
{
	cerr << "readInSam()" << endl;
	string line;

	unordered_map<string, unordered_map<int, int>> sam;
	int linecount = 0;

	while(getline(fp, line) && !fp.bad())
	{
		linecount++;
		cout << "line is " << linecount << endl;
		stringstream ss(line);
		
		string trash, chr, seqstr;
		int strand, pos;	

		if (ss >> trash >> strand >> chr >> pos >> trash >> trash >> trash >> trash >> trash >> seqstr)
		{
			int binstart = pos / sambin;
			int binend = ( pos + seqstr.length() ) / sambin;

//			cerr << "binstart [" << binstart << "] binend [" << binend << "] seqlen [" << seqstr.length() << "] pos [" << pos << "]" << endl;
			for (int bin = binstart; bin <= binend; bin++)
			{
				if ((sam[chr]).find(bin) == (sam[chr]).end())
					(sam[chr])[bin] = 1;
				else
					(sam[chr])[bin]++;
			}
		}
	}

	for (auto iter = sam.begin(); iter != sam.end(); iter++)
	{
		for (auto it = sam[iter->first].begin(); it != sam[iter->first].end(); it++)
		{
			Peak tmp;
			tmp.start = it->first * sambin;
			tmp.end = tmp.start + sambin - 1;
			tmp.value = it->second;

			(*peaks)[iter->first].push_back(tmp);
		}
	}
}


File * FileInit::getFileObj(void)
{
	switch(outfiletype)
	{	
		case WIG: return new Wig(outfile, peaks);
		break;
		case BED5: return new Bed5(outfile, peaks);
		break;
		case BEDGRAPH: return new BedGraph(outfile, peaks);
		break;
		default: throw(1);
	}
}

File::File()
{
}

File::File(std::string filename, unordered_map<string,vector<Peak>>* parsed)
{
	outfile = filename;
	peaks = parsed;

	smooth_opt = false;
	blacklist_opt = false;
	clean_opt = false;
	
	cleanfile = "";
	blacklistfile = "";
	
	smooth_win = 0;
	smooth_shift = 0;

	sort_peaks();
}

void File::output(bool ucsc)
{
	if (blacklist_opt)	
		blacklist_remove();

	if (smooth_opt)
	{
		smooth();
	
		if (blacklist_opt)
			blacklist_remove(); // have to rm them again due to smooth
	}

	if (clean_opt)
		clean();
}

void File::setOpts(bool smooth_b, int smooth_w, int smooth_s, bool blacklist_b, string blacklist_f, bool clean_b, string clean_f)
{
	blacklist_opt = blacklist_b;
	blacklistfile = blacklist_f;
		
	clean_opt = clean_b;
	cleanfile = clean_f;
	
	smooth_opt = smooth_b;
	smooth_win = smooth_w;
	smooth_shift = smooth_s;
}

void File::sort_peaks(void)
{
	for (auto iter = (*peaks).begin(); iter != (*peaks).end(); iter++)
	{
		chr_list.push_back(iter->first);			

		sort((iter->second).begin(), (iter->second).end(), peakCmp);
	}

	sort(chr_list.begin(), chr_list.end());
}

void File::smooth(void)
{
	int shift_end = smooth_shift - 1;

	for (auto iter = (*peaks).begin(); iter != (*peaks).end(); iter++)
	{
		vector<Peak> window_bins;

		int chr_start, chr_end;

		chr_start = (iter->second).begin()->start;
		chr_end = ((iter->second).end()-1)->end;

		int start_win = chr_start / smooth_win * smooth_win - smooth_win;

		if (start_win < 0)
			start_win = 0;

		double circ_ar[smooth_win / smooth_shift], curr_avg = 0;
		int circ_index = 0;	
		
		for (int i = 0; i < smooth_win / smooth_shift; i++)
			circ_ar[i] = 0;	

		vector<Peak>::iterator peakIter = (iter->second).begin();	
		for (int i = start_win; i <= ((chr_end+smooth_win)/smooth_win)*smooth_win + smooth_win; i+= smooth_shift)
		{
			circ_index = (i / smooth_shift) % (smooth_win / smooth_shift);

			if (curr_avg > 0 && i - smooth_win/2 - smooth_shift/2 >= start_win && i - smooth_win / 2 + smooth_shift/2 <= ((chr_end + smooth_win)/smooth_win)*smooth_win + smooth_win)
			{
				
				Peak tmp;
				tmp.start = (i - smooth_win/2) - smooth_shift/2;
				tmp.end = tmp.start + smooth_shift - 1;
				tmp.value = (curr_avg / (smooth_win - 1));
				window_bins.push_back(tmp);
			}

			// MUST be above continue and break lines
			curr_avg -= circ_ar[circ_index];
			circ_ar[circ_index] = 0;

			while(peakIter < (iter->second).end()-1 && peakIter->end < i)
				peakIter++;
			//   ====
			//-----

			if ((peakIter == (iter->second).end()-1 && peakIter->end < i ) || i + shift_end < peakIter->start)
				continue;

			if (peakIter->start < i)
			{
				if (peakIter->end > i + shift_end)
					circ_ar[circ_index] += peakIter->value * smooth_shift;
				else 
					circ_ar[circ_index] += peakIter->value * (peakIter->end - i);
			} // if
			else
			{
				if (peakIter->end > i + shift_end)
					circ_ar[circ_index] += peakIter->value * (i + shift_end - peakIter->start);
				else
					circ_ar[circ_index] += peakIter->value * (peakIter->end - peakIter->start);
			} // else

			curr_avg += circ_ar[circ_index];
		} // for chrstart -> chrend

		peaks->at(iter->first) = window_bins;	
	} // for each chr	
}

// expect sorted & bed file
void File::blacklist_remove(void) 
{
	cerr << "start of blacklist_remove()" << endl;
	ifstream bfile(blacklistfile);
	
	if (!bfile.is_open())
		throw(1); // throw error

	string line;
		
	unordered_map<string, vector<Peak>> blacklist;

	while(getline(bfile, line) && !bfile.bad())
	{
		string chr;
		int start, end;

		stringstream ss(line);

		if (ss >> chr >> start >> end)
		{
			Peak tmp;
			tmp.start = start;
			tmp.end = end;

			blacklist[chr].push_back(tmp);
		}
	}

	bfile.close();

	for (auto iter = (*peaks).begin(); iter != (*peaks).end(); iter++)
	{
		if (blacklist[iter->first].size() == 0)
			continue;

		vector<Peak>::iterator peakIter = (iter->second).begin();
		vector<Peak>::iterator rmIter = (iter->second).end() - 1;

		for (auto it = blacklist[iter->first].begin(); it != blacklist[iter->first].end(); it++)
		{
			if (peakIter == (iter->second).end())
				break;

			if (it->end < peakIter->start)
				continue;

			while(peakIter < (iter->second).end() - 1 && it->start > peakIter->end)
				peakIter++;

			while(peakIter != (iter->second.end()) && it->end > peakIter->start && it->start < peakIter->end)
				 peakIter = (iter->second).erase(peakIter);
		}
	}
	cerr << "end of blacklist_remove()" << endl;
}

void File::clean(void)
{
	ifstream chrfile(cleanfile);
	
	if (!chrfile.is_open())
		throw (1); // throw error

	string line;

	unordered_map<string, int> chromLengths;

	while(getline(chrfile, line) && !chrfile.bad())
	{
		string chr;
		int len;

		stringstream ss(line);

		if (ss >> chr >> len)
			chromLengths[chr] = len;		
	}

	chrfile.close();	
	
	// if chromLengths size is 0, throw error

	for (auto iter = (*peaks).begin(); iter != (*peaks).end(); iter++)
	{
		vector<Peak>::iterator eraseStart = (iter->second).end();

		for (vector<Peak>::iterator pIter = (iter->second).begin(); pIter != (iter->second).end(); pIter++)
		{
			if (pIter->end > chromLengths[iter->first])
			{
				eraseStart = pIter;		
				break;
			}		
		}

		(iter->second).erase(eraseStart, (iter->second).end());
	}		
}

void Wig::output(bool ucsc)
{		
	File::output(ucsc);
	cerr << "Wig output" << endl;

	ofstream out(outfile);

	if(!out.is_open())
		throw(1); // throw error

// FIXME this is not outputting
	if (ucsc)
	{
		out << "track type=wiggle_0 name=\"" << outfile << "\" description=\"" << outfile << "\"" << endl;
	}

	for (vector<string>::iterator iter = chr_list.begin(); iter != chr_list.end(); iter++)
	{
		int span = -1;
//		cerr << "total peaks for chr " << *iter << " are " << (*peaks)[*iter].size() << endl;
		for (vector<Peak>::iterator peakIt = (*peaks)[*iter].begin(); peakIt != (*peaks)[*iter].end(); peakIt++)
		{
			if (peakIt->end - peakIt->start != span)
			{
				span = peakIt->end - peakIt->start;
				out << "variableStep chrom=" << *iter << " span=" << span << endl;
			}
			
			out << peakIt->start << "\t" << peakIt->value << endl;
		}
	}	

	out.close(); 
}

void Bed5::output(bool ucsc)
{
	File::output(ucsc);
	cerr << "Bed5 output" << endl; 

	ofstream out(outfile);

	if(!out.is_open())
		throw(1); // throw error

	if (ucsc)
	{
		out << "track name=\"" << outfile << "\" description=\"" << outfile << "\"" << endl;
	}

	for (vector<string>::iterator iter = chr_list.begin(); iter != chr_list.end(); iter++)
	{
		for (vector<Peak>::iterator peakIt = (*peaks)[*iter].begin(); peakIt != (*peaks)[*iter].end(); peakIt++)
		{
			out << *iter << "\t" << peakIt->start << "\t" << peakIt->end << "\t.\t" << peakIt->value << endl;
		}
	}
	
	out.close();
}

void BedGraph::output(bool ucsc)
{
	File::output(ucsc);
	cerr << "BedGraph output" << endl; 

	ofstream out(outfile);

	if(!out.is_open())
		throw(1); // throw error

	if (ucsc)
	{
		out << "track type=bedGraph name=\"" << outfile << "\" description=\"" << outfile << "\"" << endl;
	}

	for (vector<string>::iterator iter = chr_list.begin(); iter != chr_list.end(); iter++)
	{
		for (vector<Peak>::iterator peakIt = (*peaks)[*iter].begin(); peakIt != (*peaks)[*iter].end(); peakIt++)
		{
			out << *iter << "\t" << peakIt->start << "\t" << peakIt->end << "\t" << peakIt->value << endl;
		}
	}
	
	out.close();
}
