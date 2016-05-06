#ifndef FILES_H
#define FILES_H

#include"Utils.h"
#include<fstream>
#include<unordered_map>
#include<algorithm>

typedef struct 
{
	int start;
	int end;
	double * value;
} Peak;

bool peakCmp(Peak peak1, Peak peak2);

class File;

class FileInit
{
public:
	FileInit(std::string infilename, std::string outfilename, Filetype type, int threads = 0);
	File * getFileObj(void);
private: 
	void readIn(std::string infilename);

	void readInWig(std::ifstream& fp);
	void readInBed(std::ifstream& fp);
	void readInBedGraph(std::ifstream& fp);
	void readInSam(std::ifstream& fp);

	std::string outfile;
	Filetype outfiletype;
	std::unordered_map<std::string, std::vector<Peak>>* peaks;
};

class File 
{
public:
	File();
	File(std::string outfilename, std::unordered_map<std::string,std::vector<Peak>>* parsed);
	void setOpts(bool smooth_b, int smooth_w, int smooth_s, bool blacklist_b, std::string blacklist_f, bool clean_b, std::string clean_f);
	virtual void output(bool ucsc);
protected:
	void sort_peaks(void);
	void smooth(void);
	void blacklist_remove(void);
	void clean(void);

	std::unordered_map<std::string, std::vector<Peak>>* peaks;
	std::vector<std::string> chr_list; // sorted list

	bool smooth_opt;	
	bool blacklist_opt;
	bool clean_opt;
		
	std::string cleanfile;
	std::string blacklistfile;
	std::string outfile;
	
	int smooth_win;	
	int smooth_shift;	
};

class Wig : public File 
{
	public:
	Wig(std::string outfilename, std::unordered_map<std::string,std::vector<Peak>>* parsed) : File(outfilename, parsed) { };
	virtual void output(bool ucsc);
	
	private:
};

class Bed : public File 
{
public:
	Bed(std::string outfilename, std::unordered_map<std::string,std::vector<Peak>>* parsed) : File(outfilename, parsed) { };
	virtual void output(bool ucsc);
private:

};

class Bed5 : public File 
{
public:
	Bed5(std::string outfilename, std::unordered_map<std::string,std::vector<Peak>>* parsed) : File(outfilename, parsed) { };
	virtual void output(bool ucsc);
private:

};

class BedGraph : public File 
{
public:
	BedGraph(std::string outfilename, std::unordered_map<std::string,std::vector<Peak>>* parsed) : File(outfilename, parsed) { };
	virtual void output(bool ucsc);
private:

};

#endif
