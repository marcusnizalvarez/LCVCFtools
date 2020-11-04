#ifndef LCVCFFILTER_H
#define LCVCFFILTER_H
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <chrono>
#include <ctime>
#include <cmath>
/************************************************************
Copyright Joe Coder 2004 - 2006.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file BOOST.license or copy at
https://www.boost.org/LICENSE_1_0.txt)
************************************************************/
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#define THIS_VERSION "1.0.0"
struct STATSstruct{
    std::vector<unsigned> DP, GQ;
    double GCR, NMR;
    std::string Name;
};
struct SAMPLEstruct{
    std::vector<std::string> FORMAT;
    std::string *GT, *PL, *DP, *AD, *GQ;
};
class LCVCFtools
{
public:
    LCVCFtools(int argc, char* argv[]);
    void Run();
private:
    /*************
        FUNCTIONS
    *************/
    void CheckParameters();
    void CheckARG(std::string Argument, bool IsUnique = true);
    void Log(std::string Msg);
    void OutputHeader();
    void OutputLine();
    void ReadData();
    void ReadHeader();
    void ReadRemoveList();
    void SetInputMode();
    void SetParameters(std::vector<std::string>& args);
    void ShowHelp();
    void ShowStatus();
    void OutputSampleStatistics();
    void OutputOtherStatistics();
    void Terminate(std::string Msg = "", bool Help = false, int ReturnValue = -1);
    bool CheckRate(const std::vector<int> &vec, int val, double qnt);
    bool Filter();
    bool GetLine(std::string& TmpString);
    bool StringToVcf(const std::string& tmpLineString);
    int StrToData(const std::string& String);
    std::string NowString();
    std::vector<std::string> SplitString(const std::string&, char Delim);
    /*************
        VARIABLES
    *************/
    std::set<std::string> DefinedArguments;
    std::string Parameters;
    std::string StartingTimeStr;
    /*************
        INPUT
    *************/
    std::ifstream file;
    boost::iostreams::filtering_istream in;
    std::string InputFilename;
    std::string RemoveFilename;
    bool IsGzipped;
    bool IsFile;
    /*************
        FILTER PARAMETERS
    *************/
    int    minDP;
    int    minGQ;
    double minGCR;
    std::vector<int>    DPRlevel, GQRlevel;
    std::vector<double> DPRvalue, GQRvalue;
    double MAF;
    /*************
        OTHER PARAMETERS
    *************/
    bool IsOtherStats;
    bool IsID;
    bool IsSampleStats;
    std::set<std::string> RemoveSamples;
    std::set<size_t> RemoveIndex;
    /*************
        VCF HEADER
    *************/
    std::vector<std::string> CommentLines, HeaderColumns, HeaderSamples;
    /*************
        VCF LINE DATA
    *************/
    std::string CHR, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMATstr, lastFORMATstr;
    std::vector<SAMPLEstruct> SAMPLES;
    std::map<std::string,size_t> FORMATtagsMap;
    std::vector<std::string> FORMATtags;
    std::vector<std::pair<short,double>> Alleles;
    /*************
        STATUS
    *************/
    size_t LogFreq;
    size_t tmpCounter;
    size_t OutputtedLines;
    size_t ReadedLines;
    size_t RemovedDueGCR;
    size_t RemovedDueMAF;
    std::vector<size_t> RemovedDueDPR, RemovedDueGQR;
    /*************
        OTHER STATS
    *************/
    double GCR;
    /*************
        SAMPLE STATS
    *************/
    std::ofstream StatisticsFile;
    std::vector<STATSstruct> StatsVector;
    int SampleStatsYLim;
};
#endif
