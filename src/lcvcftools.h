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
#include <numeric>
/************************************************************
Copyright Joe Coder 2004 - 2006.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file BOOST.license or copy at
https://www.boost.org/LICENSE_1_0.txt)
************************************************************/
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#define THIS_VERSION "1.0.4"
class LCVCFtools
{
public:
    LCVCFtools(int argc, char* argv[]);
    LCVCFtools();
    void Run();
private:
    struct SampleStatsStruct{
        std::vector<int> Depth, Quality;
        double NonMissingRate, MeanDepth, MeanQuality;
        std::string SampleId;
    };
    struct SampleDataStruct{
        std::vector<std::string> FORMAT;
        std::string *GT, *PL, *DP, *AD, *GQ;
    };
    struct SnpDataStruct{
        std::string CHR, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMATstr;
        std::vector<SampleDataStruct> SampleDataVector;
        std::vector<std::pair<int,double>> AlleleCountVector;
        double GCR=0;
    };
    /*************
        FUNCTIONS
    *************/
    void CheckARG(std::string Argument, bool IsUnique = true);
    void Log(std::string Msg);
    void OutputHeader();
    void OutputLine();
    void ReadData();
    void ReadHeader();
    void ReadKeepList(std::string Filename);
    void ReadRemoveList(std::string Filename);
    void OpenInputStream();
    void SetParameters(std::vector<std::string>& args);
    void ShowHelp();
    void ShowProgress();
    void OutputSampleStatistics();
    void Terminate(std::string Msg);
    bool CheckRate(const std::vector<int> &vec, int val, double qnt);
    bool Filter();
    bool GetLine(std::string& TmpString);
    bool StringToVcf(const std::string& tmpLineString);
    int StringToInt(const std::string& String);
    std::string GetParametersString();
    std::string NowString();
    /*************
        VARIABLES
    *************/
    std::set<std::string> DefinedArguments;
    std::string StartingTimeStr;
    /*************
        INPUT
    *************/
    std::ifstream file;
    size_t filesize;
    boost::iostreams::filtering_istream in;
    std::string InputFilename;
    bool IsGzipped = false;
    bool IsFile = true;
    /*************
        FILTER PARAMETERS
    *************/
    int    minDP = 5;
    int    minGQ = 20;
    double minGCR = 0;
    double MAF = 0.1;
    std::vector<int>    DPRlevel, GQRlevel;
    std::vector<double> DPRvalue, GQRvalue;
    /*************
        OTHER PARAMETERS
    *************/
    bool IsVerbose = false;
    bool IsID = false;
    bool IsSampleStats = false;
    bool IsRemoveMultiallelic = true;
    std::set<std::string> RemoveSamples, KeepSamples;
    std::set<size_t> RemoveIndex;
    /*************
        VCF DATA
    *************/
    std::vector<std::string> CommentLines, HeaderColumns, HeaderSamples;
    std::string lastFORMATstr;
    std::map<std::string,size_t> FORMATtagsMap;
    std::vector<std::string> FORMATtagsVector;
    /*************
        TEMP SNP DATA
    *************/
    SnpDataStruct tmpSnpData;
    /*************
        STATUS
    *************/
    size_t InputCounter = 0;
    size_t OutputCounter = 0;
    size_t RemovedGenotypeCallRate = 0;
    size_t RemovedMultiallelic = 0;
    size_t RemovedIndels = 0;
    size_t RemovedMAF = 0;
    size_t Verbosity = 10000;
    std::vector<size_t> RemovedDepthRate;
    std::vector<size_t> RemovedQualityRate;
    /*************
        STATS
    *************/
    std::ofstream SampleStatsFile;
    std::vector<SampleStatsStruct> SampleStatsVector;
    int YLim = 100;
    double YLimThreshold = 0.001;
};
#endif
