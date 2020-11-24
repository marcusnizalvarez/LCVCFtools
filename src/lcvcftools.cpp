#include "lcvcftools.h"
void LCVCFtools::ShowHelp(){
    std::clog << "LCVCFtools Version " << THIS_VERSION << " is under GNU GPLv3.0\n" <<
                 "Documentation and source code available at https://github.com/marcusnizalvarez/LCVCFtools\n"
                 "If you would like more information, please contact me at marcus.alvarez@unesp.br\n"
                 "\n"
                 "[Input mode] \n"
                 "--vcf    <STRING>       Read from VCF file. Use - to read from stdin.\n"
                 "--gzvcf  <STRING>       Read from Gzip compressed VCF file. Use - to read from stdin.\n"
                 "\n"
                 "[Filter parameters] \n"
                 "--minGQ  <INT>          Minimum genotype quality in PhredScale. [Default=10]\n"
                 "--minDP  <INT>          Minimum depth. [Default=1]\n"
                 "--MAF    <FLOAT>        Minor allele frequency, based on allele depth (AD). [Default=0.01]\n"
                 "--minGCR <FLOAT>        Minimum genotype call rate [Default=0].\n"
                 "--minDPR <INT> <FLOAT>  Minimum DP rate. Can be defined multiple times.\n"
                 "--minGQR <INT> <FLOAT>  Minimum GQ rate. Can be defined multiple times.\n"
                 "\n"
                 "[Other arguments] \n"
                 "--remove <STRING>       Remove samples from a list file (One sample name per line).\n"
                 "--keep <STRING>         Keep samples from a list file (One sample name per line).\n"
                 "--sample-stats          Output sample statistics to 'stats1.tsv'.\n"
                 "--other-stats           Output other statistics to 'stats2.tsv'. (AD-based freq, GCR,...)\n"
                 "--ID                    Generate generic ID, useful for programs like Plink.\n"
                 "--verbose               Verbose mode.\n"
                 "--help                  Print this message.\n"
              << std::endl;
}
LCVCFtools::LCVCFtools(int argc, char* argv[]){
    std::clog << std::fixed << std::setprecision(2);
    /*************
        Initialize Variables
    *************/
    VerboseFreq = 10000;
    SampleStatsYLim = 30;
    tmpCounter = 0;
    ReadedLines = 0;
    OutputtedLines = 0;
    RemovedDueGCR = 0;
    RemovedDueMAF = 0;
    /*************
        Initialize filter parameters
    *************/
    minGQ = 10;
    minDP = 1;
    MAF = 0.01;
    minGCR = 0;
    /*************
        Initialize other parameters
    *************/
    IsVerbose = false;
    IsFile = true;
    IsGzipped = false;
    IsID = false;
    IsOtherStats = false;
    IsSampleStats = false;
    /*************
        Read Args
    *************/
    std::vector<std::string> args;
    for(int i(1); i < argc; i++) args.push_back(std::string(argv[i]));
    SetParameters(args);
    SetInputMode();
    /*************
        Open statistics file
    *************/
    if(IsSampleStats){
        SampleStatsFile.open("stats1.tsv");
        SampleStatsFile << std::scientific << std::setprecision(2);
    }
    if(IsOtherStats){
        OtherStatsFile.open("stats2.tsv");
        OtherStatsFile << std::scientific << std::setprecision(2);
    }
}
void LCVCFtools::ReadRemoveList(){
    if(RemoveFilename.empty()) return;
    std::ifstream File(RemoveFilename);
    if(!File.is_open()) Terminate("Remove list file does not exist or is not readable.", false);
    while(true){
        std::string tmpString;
        if(!std::getline(File, tmpString)) break;
        if(!tmpString.empty()) RemoveSamples.insert(tmpString);
    }
}
void LCVCFtools::ReadKeepList(){
    if(KeepFilename.empty()) return;
    std::ifstream File(KeepFilename);
    if(!File.is_open()) Terminate("Keep list file does not exist or is not readable.", false);
    while(true){
        std::string tmpString;
        if(!std::getline(File, tmpString)) break;
        if(!tmpString.empty()) KeepSamples.insert(tmpString);
    }
}
void LCVCFtools::CheckParameters(){
    std::ostringstream sso;
    sso << std::fixed;
    sso << std::setprecision(2);
    sso << "minGQ=" << minGQ << ";";
    sso << "minDP=" << minDP << ";";
    sso << "minGCR=" << minGCR << ";";
    sso << "MAF=" << MAF << ";";
    for(size_t i(0); i < DPRlevel.size(); i++) sso  << "minDPR" << i+1 << "=[" << DPRlevel[i] << ';' << DPRvalue[i] << "];";
    for(size_t i(0); i < GQRlevel.size(); i++) sso  << "minGQR" << i+1 << "=[" << GQRlevel[i] << ';' << GQRvalue[i] << "];";
    Parameters = sso.str();
}
void LCVCFtools::Run(){
    StartingTimeStr = NowString();
    Log("Running LCVCFtools v" + std::string(THIS_VERSION));
    CheckParameters();
    Log("Running with parameters: " + Parameters);
    Log("Starting...");
    ReadHeader();
    ReadData();
    ShowStatus();
    OutputSampleStatistics();
    Log("Finished.");
}
void LCVCFtools::SetParameters(std::vector<std::string>& args){
    if(args.empty()) Terminate("",true,0);
    for(size_t i(0); i < args.size(); i++){
        if(args[i]=="--vcf"){
            CheckARG("input");
            if(++i >= args.size()) Terminate("Missing argument value for input file", true);
            InputFilename = args[i];
            if(InputFilename=="-") IsFile = false;
            continue;
        }
        if(args[i]=="--gzvcf"){
            CheckARG("input");
            if(++i >= args.size()) Terminate("Missing argument value for input file", true);
            InputFilename = args[i];
            if(InputFilename=="-") IsFile = false;
            IsGzipped = true;
            continue;
        }
        if(args[i]=="--minGQ"){
            CheckARG("minGQ");
            if(++i >= args.size()) Terminate("Missing argument value for minGQ", true);
            minGQ = std::stoi(args[i]);
            if(minGQ <= 0) Terminate("minGQ must be greater than 0", true);
            continue;
        }
        if(args[i]=="--minDP"){
            CheckARG("minDP");
            if(++i >= args.size()) Terminate("Missing argument value for minDP", true);
            minDP = std::stoi(args[i]);
            if(minDP <= 0) Terminate("minDP must be greater than 0", true);
            continue;
        }
        if(args[i]=="--minGCR"){
            CheckARG("minGCR");
            if(++i >= args.size()) Terminate("Missing argument value for minGCR", true);
            minGCR = std::stod(args[i]);
            if(minGCR < 0 || minGCR > 1) Terminate("minGCR must be between 0 and 1", true);
            continue;
        }
        if(args[i]=="--MAF"){
            CheckARG("MAF");
            if(++i >= args.size()) Terminate("Missing argument value for MAF", true);
            MAF = std::stod(args[i]);
            if(MAF < 0 || MAF > 1) Terminate("MAF must be between 0 and 1", true);
            continue;
        }
        if(args[i]=="--minDPR"){
            CheckARG("minDPR", false);
            if(i+2 >= args.size()) Terminate("Missing arguments value for minDPR", true);
            if(std::stoi(args[i+1]) < 0) Terminate("Level for minDPR must be greater than 0", true);
            if(std::stod(args[i+2]) < 0 || std::stod(args[i+2]) > 1) Terminate("Rate for minDPR must be between 0 and 1", true);
            DPRlevel.push_back(std::stoi(args[++i]));
            DPRvalue.push_back(std::stod(args[++i]));
            RemovedDueDPR.push_back(0);
            continue;
        }
        if(args[i]=="--minGQR"){
            CheckARG("minGQR", false);
            if(i+2 >= args.size()) Terminate("Missing arguments value for minGQR", true);
            if(std::stoi(args[i+1]) < 0) Terminate("Level for minGQR must be greater than 0", true);
            if(std::stod(args[i+2]) < 0 || std::stod(args[i+2]) > 1) Terminate("Rate for minGQR must be between 0 and 1", true);
            GQRlevel.push_back(std::stoi(args[++i]));
            GQRvalue.push_back(std::stod(args[++i]));
            RemovedDueGQR.push_back(0);
            continue;
        }
        if(args[i]=="--sample-stats"){
            CheckARG("sample-stats");
            IsSampleStats = true;
            continue;
        }
        if(args[i]=="--other-stats"){
            CheckARG("other-stats");
            IsOtherStats = true;
            continue;
        }
        if(args[i]=="--ID"){
            CheckARG("ID");
            IsID = true;
            continue;
        }
        if(args[i]=="--remove"){
            CheckARG("remove");
            if(!KeepFilename.empty()) Terminate("Argument remove cannot be used with keep");
            if(++i >= args.size()) Terminate("Missing argument value for remove", true);
            RemoveFilename = args[i];
            continue;
        }
        if(args[i]=="--keep"){
            CheckARG("keep");
            if(!RemoveFilename.empty()) Terminate("Argument keep cannot be used with remove");
            if(++i >= args.size()) Terminate("Missing argument value for keep", true);
            KeepFilename = args[i];
            continue;
        }
        if(args[i]=="--verbose"){
            CheckARG("verbose");
            IsVerbose=true;
            continue;
        }
        if(args[i]=="--help") Terminate("",true,0);
        Terminate(args[i] + " is an invalid argument", true);
    }
    if(DefinedArguments.find("input")==DefinedArguments.end())
        Terminate("Missing input mode argument", true);
}
void LCVCFtools::SetInputMode(){
    //          Copyright Joe Coder 2004 - 2006.
    // Distributed under the Boost Software License, Version 1.0.
    //    (See accompanying file BOOST.license or copy at
    //          https://www.boost.org/LICENSE_1_0.txt)
    if(IsFile){
        if(IsGzipped){
            file.open(InputFilename, std::ios_base::in | std::ios_base::binary);
            if(!file.is_open()) Terminate("VCF file does not exist or is not readable", true);
            in.push(boost::iostreams::gzip_decompressor());
            in.push(file);
        }
        else{
            file.open(InputFilename);
            if(!file.is_open()) Terminate("VCF file does not exist or is not readable", true);
        }
    }
    else{
        if(IsGzipped){
            in.push(boost::iostreams::gzip_decompressor());
            in.push(std::cin);
        }
    }
}
void LCVCFtools::Terminate(std::string Msg, bool Help, int ReturnValue){
    if(Help) ShowHelp();
    if(!Msg.empty()) Log("\033[1;31m**ERROR** " + Msg + "\033[0m");
    throw ReturnValue;
}
void LCVCFtools::ShowStatus(){
    std::clog << NowString() << "Input=" << ReadedLines << ";Filtered:{"
              << "GCR=" << (static_cast<double>(RemovedDueGCR)/ReadedLines)*100 << "%;"
              << "MAF=" << (static_cast<double>(RemovedDueMAF)/ReadedLines)*100 << "%;";
    for(size_t i(0);i < RemovedDueDPR.size(); i++)
        std::clog << "DPR[" << DPRlevel[i] << "," << DPRvalue[i] << "]="
                  << (static_cast<double>(RemovedDueDPR[i])/ReadedLines)*100
                  << "%;";
    for(size_t i(0);i < RemovedDueGQR.size(); i++)
        std::clog << "GQR[" << GQRlevel[i] << "," << GQRvalue[i] << "]="
                  << (static_cast<double>(RemovedDueGQR[i])/ReadedLines)*100
                  << "%;";
    std::clog << "};Output=" << OutputtedLines << "(" << (static_cast<double>(OutputtedLines)/ReadedLines)*100 << "%);..."
              << std::endl;
}
void LCVCFtools::OutputSampleStatistics(){
    if(!IsSampleStats) return;
    Log("Calculating sample statistics...");
    for(size_t i(0); i < HeaderSamples.size(); i++) StatsVector[i].Name = HeaderSamples[i];
    std::sort(StatsVector.begin(),StatsVector.end(),[](const STATSstruct& a, const STATSstruct& b)->bool{return a.NMR > b.NMR;});
    SampleStatsFile << "## Sample Statistics Table" << std::endl;
    SampleStatsFile << "## NMR=Mean non-missing rate; GCR=Mean genotype call rate; DP=Depth at a given level; GQ=Genotype at a given level" << std::endl;
    SampleStatsFile << "Sample\tVariable\tLevel\tValue" << std::endl;
    for(STATSstruct& x : StatsVector){
        int tmpDPsum(x.DP[0]), tmpGQsum(x.GQ[0]);
        for(int i(1); i < SampleStatsYLim+1; i++){
            tmpDPsum += x.DP[i];
            tmpGQsum += x.GQ[i];
            for(int j(i+1); j < SampleStatsYLim+1; j++){
                x.DP[i] += x.DP[j];
                x.GQ[i] += x.GQ[j];
            }
        }
        SampleStatsFile << x.Name << "\t" << "NMR" << "\t.\t" << x.NMR/OutputtedLines  << std::endl;
        SampleStatsFile << x.Name << "\t" << "GCR" << "\t.\t" << x.GCR/OutputtedLines  << std::endl;
        for(int i(1); i < SampleStatsYLim+1; i++)
            SampleStatsFile << x.Name << "\t" << "DP" << "\t" << i << "\t" << static_cast<double>(x.DP[i])/tmpDPsum  << std::endl;
        for(int i(1); i < SampleStatsYLim+1; i++)
            SampleStatsFile << x.Name << "\t" << "GQ" << "\t" << i << "\t" << static_cast<double>(x.GQ[i])/tmpGQsum  << std::endl;
    }
}
void LCVCFtools::OutputOtherStatistics(){
    if(!OutputtedLines) {
        OtherStatsFile << "## Other Statistics Table" << std::endl;
        OtherStatsFile << "## ID=Identity; AFO=Allele Order (Decreasing); AF=Allele Frequency; GCR=Genotype Call Rate" << std::endl;
        OtherStatsFile << "ID\tAO\tAF\tGCR" << std::endl;
    }
    OtherStatsFile << CHR + ":" + POS + ':' + REF + ':' + ALT;
    OtherStatsFile << '\t' << Alleles[0].first;
    for(size_t i(1); i < Alleles.size(); i++) OtherStatsFile << ';' << Alleles[i].first;
    OtherStatsFile << '\t' << Alleles[0].second;
    for(size_t i(1); i < Alleles.size(); i++) OtherStatsFile << ';' << Alleles[i].second;
    OtherStatsFile << '\t' << GCR;
    OtherStatsFile << std::endl;
}
bool LCVCFtools::GetLine(std::string& TmpString){
    if(IsFile){
        if(IsGzipped){
            if(!std::getline(in, TmpString, '\n')) return false;
        }
        else{
            if(!std::getline(file, TmpString, '\n')) return false;
        }
    }
    else{
        if(IsGzipped){
            if(!std::getline(in, TmpString, '\n')) return false;
        }
        else{
            if(!std::getline(std::cin, TmpString, '\n')) return false;
        }
    }
    return true;
}
void LCVCFtools::ReadHeader(){
    std::string CantReadMsg("Can't read VCF file, check the file format.");
    while(true){
        std::string temp_string;
        if(!GetLine(temp_string)) Terminate(CantReadMsg,false);
        if(temp_string[0]=='#'){
            if(temp_string.substr(0,6)=="#CHROM"){
                std::istringstream temp_ssin(temp_string);
                {
                    std::vector<std::string> tmpColumnCheck ({"#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"});
                    std::string tmpString2;
                    for(int b(0); b < 9; b++){
                        if(!std::getline(temp_ssin, tmpString2, '\t')) Terminate(CantReadMsg,false);
                        if(tmpString2!=tmpColumnCheck[b]) Terminate(tmpString2 + " column must be " + tmpColumnCheck[b], false);
                        HeaderColumns.push_back(tmpString2);
                    }
                }
                /******* Apply remove/keep OPT *******/
                ReadRemoveList();
                ReadKeepList();
                size_t i(0);
                while(true){
                    std::string tmpName;
                    if(!std::getline(temp_ssin, tmpName, '\t')) break;
                    bool Rm=false;
                    if(!RemoveSamples.empty())
                        if(RemoveSamples.find(tmpName)!=RemoveSamples.end()){
                            RemoveIndex.insert(i);
                            Rm=true;
                        }
                    if(!KeepSamples.empty())
                        if(KeepSamples.find(tmpName)==KeepSamples.end()) {
                            RemoveIndex.insert(i);
                            Rm=true;
                        }
                    if(!Rm){
                        HeaderSamples.push_back(tmpName);
                        STATSstruct tmpSTATS;
                        tmpSTATS.DP.resize(SampleStatsYLim+1,0);
                        tmpSTATS.GQ.resize(SampleStatsYLim+1,0);
                        StatsVector.push_back(std::move(tmpSTATS));
                    }
                    i++;
                }
                Log(std::to_string(HeaderSamples.size()) + " samples identified...");
                if(RemoveSamples.size()>0)
                    Log("Removing " + std::to_string(RemoveIndex.size()) + " samples...");
                break;
            }
            CommentLines.push_back(temp_string);
        }
        else Terminate(CantReadMsg, false);
    }
    if(!HeaderSamples.size()) Terminate("No samples in VCF file.",false);
    OutputHeader();
}
void LCVCFtools::ReadData(){
    while(true){
        std::string tmpLineString;
        if(!GetLine(tmpLineString)) break;
        if(tmpLineString[0]=='#') Terminate("Comment at line "+std::to_string(ReadedLines)+".",false);
        ReadedLines++;
        if(!StringToVcf(tmpLineString)) Terminate("Failed to read data at line " + std::to_string(ReadedLines)+", check the file format", false);
        if(Filter()){
            OutputtedLines++;
            OutputLine();
        }
        if(IsVerbose && (++tmpCounter >= VerboseFreq)){
            tmpCounter = 0;
            ShowStatus();
        }
    }
}
std::vector<std::string> LCVCFtools::SplitString(const std::string& String, char Delim){
    std::vector<std::string> tmpVector;
    std::istringstream ssin(String);
    while(true){
        std::string tmpString;
        if(!std::getline(ssin, tmpString, Delim)) break;
        tmpVector.push_back(tmpString);
    }
    return tmpVector;
}
bool LCVCFtools::StringToVcf(const std::string& tmpLineString){
    std::istringstream ssin (tmpLineString);
    if(!std::getline(ssin, CHR, '\t')) return false;
    if(!std::getline(ssin, POS, '\t')) return false;
    if(!std::getline(ssin, ID, '\t')) return false;
    if(!std::getline(ssin, REF, '\t')) return false;
    if(!std::getline(ssin, ALT, '\t')) return false;
    if(!std::getline(ssin, QUAL, '\t')) return false;
    if(!std::getline(ssin, FILTER, '\t')) return false;
    if(!std::getline(ssin, INFO, '\t')) return false;
    if(!std::getline(ssin, FORMATstr, '\t')) return false;
    if(lastFORMATstr!=FORMATstr){
        FORMATtagsMap.clear();
        FORMATtags.clear();
        size_t i = 0;
        std::set<std::string> RequiredTags = {"GQ","DP","AD","GT","PL"};
        for(const std::string& String : SplitString(FORMATstr,':')){
            FORMATtags.push_back(String);
            FORMATtagsMap[String] = i++;
            RequiredTags.erase(String);
        }
        if(!RequiredTags.empty()){
            std::string tmpMsg = "Missing VCF required tag(s): ";
            for(std::string s : RequiredTags) tmpMsg += s + "; ";
            Terminate(tmpMsg, false);
        }
    }
    lastFORMATstr = FORMATstr;
    SAMPLES.clear();
    size_t i = 0;
    while(true){
        SAMPLEstruct tmpSample;
        std::string tmpString;
        if(!std::getline(ssin, tmpString, '\t')) break;
        if(tmpString.empty()) return false;
        if(RemoveIndex.find(i++)!=RemoveIndex.end()) continue;
        for(const auto& x : SplitString(tmpString,':')) tmpSample.FORMAT.push_back(x);
        if(tmpSample.FORMAT.size()!=FORMATtags.size()) Terminate("Incorrect number of fields for FORMAT", false);
        tmpSample.GT = &tmpSample.FORMAT[FORMATtagsMap["GT"]];
        tmpSample.PL = &tmpSample.FORMAT[FORMATtagsMap["PL"]];
        tmpSample.DP = &tmpSample.FORMAT[FORMATtagsMap["DP"]];
        tmpSample.AD = &tmpSample.FORMAT[FORMATtagsMap["AD"]];
        tmpSample.GQ = &tmpSample.FORMAT[FORMATtagsMap["GQ"]];
        SAMPLES.push_back(std::move(tmpSample));
    }
    Alleles.clear();
    for(size_t i(0); i < SplitString(ALT,',').size()+1; i++)
        Alleles.push_back(std::pair<short,double>(i,0));
    return true;
}
void LCVCFtools::OutputHeader(){
    for(const auto& tmpString : CommentLines) std::cout << tmpString << '\n';
    std::cout << "##LCVCFtools_v"
              << THIS_VERSION
              << " "
              << Parameters
              << "Date="
              << StartingTimeStr
              << '\n';
    std::string tmpHeaderString;
    for(const auto& tmpString : HeaderColumns) tmpHeaderString += tmpString + '\t';
    for(const auto& tmpString : HeaderSamples) tmpHeaderString += tmpString + '\t';
    std::cout << tmpHeaderString.erase(tmpHeaderString.size()-1) << '\n';
}
void LCVCFtools::OutputLine(){
    std::string tmpString =
                 CHR    + '\t' +
                 POS    + '\t' +
                 ID     + '\t' +
                 REF    + '\t' +
                 ALT    + '\t' +
                 QUAL   + '\t' +
                 FILTER + '\t' +
                 INFO   + '\t';
    bool tmpBl = true;
    for(const auto& TAG : FORMATtags){
        if(tmpBl) tmpBl = false;
        else tmpString += ':';
        tmpString += TAG;
    }
    for(const auto& SAMPLE : SAMPLES){
        tmpBl = true;
        tmpString += '\t';
        for(const auto& FORMAT : SAMPLE.FORMAT){
            if(tmpBl) tmpBl = false;
            else tmpString += ':';
            tmpString += FORMAT;
        }
    }
    std::cout << tmpString << '\n';
}
int LCVCFtools::StrToData(const std::string &String){
    if(String.find('.')!=std::string::npos) return -1;
    return std::stoi(String);
}
bool LCVCFtools::CheckRate(const std::vector<int> &vec, int val, double qnt){
    double tmp(0);
    for(const auto& i : vec) if(i>=val) tmp++;
    tmp /= vec.size();
    return(tmp<qnt);
}
void LCVCFtools::CheckARG(std::string Argument, bool IsUnique){
    if(IsUnique)
        if(DefinedArguments.find(Argument)!=DefinedArguments.end())
            Terminate("Multiple " + Argument + " definition", true);
    DefinedArguments.insert(Argument);
}
void LCVCFtools::Log(std::string Msg){
    std::clog << NowString() << Msg << std::endl;
}
std::string LCVCFtools::NowString(){
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeinfo);
    return(std::string("[") + buffer + "] ");
}
bool LCVCFtools::Filter(){
    GCR=0;
    std::vector<int> tmpGQ, tmpDP;
    for(SAMPLEstruct& Sample : SAMPLES){
        int DP = StrToData(*Sample.DP);
        int GQ = StrToData(*Sample.GQ);
        tmpDP.push_back(DP);
        tmpGQ.push_back(GQ);
        if(DP==0){
            *Sample.GT = "./.";
            //*Sample.DP = "0";
            *Sample.GQ = "0";
            tmpGQ[tmpGQ.size()-1]=0;
            int PLsep(0), ADsep(0);
            for(const char& c : *Sample.PL) if(c==',') PLsep++;
            *Sample.PL = '0';
            while(PLsep--) *Sample.PL += ",0";
            for(const char& c : *Sample.AD) if(c==',') ADsep++;
            *Sample.AD = '0';
            while(ADsep--) *Sample.AD += ",0";
            continue;
        }
        else{
            /******* Apply minDP FILTER *******/
            if(DP<minDP) *Sample.GT = "./.";
            /******* Apply minGQ FILTER *******/
            if(GQ<minGQ) *Sample.GT = "./.";
            if(DP>=minDP && GQ>=minGQ) GCR++;
            std::vector<int> AD;
            int tmpADsum(0);
            for(std::string ADstr : SplitString(*Sample.AD,',')) {
                AD.push_back(StrToData(ADstr));
                tmpADsum += AD[AD.size()-1];
            }
            if(AD.size()!=Alleles.size()) Terminate("Fatal error at AD.size()!=AlleleCount.size()", false);
            if(tmpADsum>0)
                for(size_t i(0); i<AD.size(); i++)
                    Alleles[i].second += static_cast<double>(AD[i])/tmpADsum;
        }
    }// for Sample END_HERE
    GCR /= SAMPLES.size();
    /******* Apply minGCR FILTER *******/
    if(GCR < minGCR) {
        RemovedDueGCR++;
        return false;
    }
    /******* Apply minDPR FILTER *******/
    for(size_t i(0); i < DPRlevel.size();i++){
        if(CheckRate(tmpDP,DPRlevel[i],DPRvalue[i])){
            RemovedDueDPR[i]++;
            return false;
        }
    }
    /******* Apply minGQR FILTER *******/
    for(size_t i(0); i < GQRlevel.size();i++){
        if(CheckRate(tmpGQ,GQRlevel[i],GQRvalue[i])){
            RemovedDueGQR[i]++;
            return false;
        }
    }
    /******* Apply MAF FILTER *******/
    double AlleleSum(0);
    for(const auto& A : Alleles) AlleleSum += A.second;
    if(AlleleSum==0) return false;
    for(auto& A : Alleles) A.second /= AlleleSum;
    std::sort(Alleles.begin(),Alleles.end(),
              [](const std::pair<short,double>& a, const std::pair<short,double>& b)->bool{return a.second > b.second;});
    if((1-Alleles[0].second)<MAF){
        RemovedDueMAF++;
        return false;
    }
    /******* Apply ID *******/
    if(IsID) ID = CHR + ":" + POS + ':' + REF + ':' + ALT;
    /******* Apply Other Stats *******/
    if(IsOtherStats) OutputOtherStatistics();
    /******* Apply Sample Stats *******/
    if(IsSampleStats)
        for(size_t i(0); i<StatsVector.size(); i++){
            if(tmpDP[i]>SampleStatsYLim) tmpDP[i] = SampleStatsYLim;
            if(tmpGQ[i]>SampleStatsYLim) tmpGQ[i] = SampleStatsYLim;
            StatsVector[i].DP[tmpDP[i]]++;
            StatsVector[i].GQ[tmpGQ[i]]++;
            StatsVector[i].NMR += (tmpDP[i]>0)?1:0;
            StatsVector[i].GCR += (tmpDP[i]>=minDP && tmpGQ[i]>=minGQ)?1:0;
        }
    return true;
}
int main(int argc, char* argv[]){
    try{
        LCVCFtools a(argc, argv);
        a.Run();
    }
    catch(int i){return i;}
    return 0;
}
