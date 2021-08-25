#include "lcvcftools.h"
LCVCFtools::LCVCFtools(){
}
LCVCFtools::LCVCFtools(int argc, char* argv[]){
    /*************
        Read Args
    *************/
    std::vector<std::string> args;
    for(int i(1); i < argc; i++) args.push_back(std::string(argv[i]));
    SetParameters(args);
    OpenInputStream();
    /*************
        Open statistics file
    *************/
    if(IsSampleStats){
        SampleStatsFile.open("stats1.tsv");
        SampleStatsFile << std::scientific << std::setprecision(5);
    }
}
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
                 "--minGQ  <INT>          Minimum genotype quality in PhredScale. [Default=20]\n"
                 "--minDP  <INT>          Minimum depth. [Default=5]\n"
                 "--MAF    <FLOAT>        Minor allele frequency, based on allele depth (AD). [Default=0.1]\n"
                 "--minGCR <FLOAT>        Minimum Genotype Call rate. [Default=0]\n"
                 "--minDPR <INT> <FLOAT>  Minimum Depth rate.\n"
                 "--minGQR <INT> <FLOAT>  Minimum Genotype Quality rate.\n"
                 "\n"
                 "[Other arguments] \n"
                 "--remove <STRING>       Remove samples listed in a file (One ID per line).\n"
                 "--sample-stats          Output sample statistics to 'stats1.tsv'.\n"
                 "--keep-multiallelic     Don't skip multiallelic variants.\n"
                 "--ID                    Generate generic ID, useful for programs like Plink.\n"
                 "--verbose               Verbose mode.\n"
                 "--help                  Print this message.\n"
              << std::endl;
}
void LCVCFtools::SetParameters(std::vector<std::string>& args){
    if(args.empty()) {
        ShowHelp();
        throw 0;
    }
    for(size_t i(0); i < args.size(); i++){
        if(args[i]=="--vcf"){
            CheckARG("input");
            if(++i >= args.size()) Terminate("Missing argument value for input file");
            InputFilename = args[i];
            if(InputFilename=="-") IsFile = false;
            continue;
        }
        if(args[i]=="--gzvcf"){
            CheckARG("input");
            if(++i >= args.size()) Terminate("Missing argument value for input file");
            InputFilename = args[i];
            if(InputFilename=="-") IsFile = false;
            IsGzipped = true;
            continue;
        }
        if(args[i]=="--minGQ"){
            CheckARG("minGQ");
            if(++i >= args.size()) Terminate("Missing argument value for minGQ");
            minGQ = std::stoi(args[i]);
            if(minGQ <= 0) Terminate("minGQ must be greater than 0");
            continue;
        }
        if(args[i]=="--minDP"){
            CheckARG("minDP");
            if(++i >= args.size()) Terminate("Missing argument value for minDP");
            minDP = std::stoi(args[i]);
            if(minDP <= 0) Terminate("minDP must be greater than 0");
            continue;
        }
        if(args[i]=="--minGCR"){
            CheckARG("minGCR");
            if(++i >= args.size()) Terminate("Missing argument value for minGCR");
            minGCR = std::stod(args[i]);
            if(minGCR < 0 || minGCR > 1) Terminate("minGCR must be between 0 and 1");
            continue;
        }
        if(args[i]=="--MAF"){
            CheckARG("MAF");
            if(++i >= args.size()) Terminate("Missing argument value for MAF");
            MAF = std::stod(args[i]);
            if(MAF < 0 || MAF > 1) Terminate("MAF must be between 0 and 1");
            continue;
        }
        if(args[i]=="--minDPR"){
            CheckARG("minDPR", false);
            if(i+2 >= args.size()) Terminate("Missing arguments value for minDPR");
            if(std::stoi(args[i+1]) < 0) Terminate("Level for minDPR must be greater than 0");
            if(std::stod(args[i+2]) < 0 || std::stod(args[i+2]) > 1) Terminate("Rate for minDPR must be between 0 and 1");
            DPRlevel.push_back(std::stoi(args[++i]));
            DPRvalue.push_back(std::stod(args[++i]));
            RemovedDepthRate.push_back(0);
            continue;
        }
        if(args[i]=="--minGQR"){
            CheckARG("minGQR", false);
            if(i+2 >= args.size()) Terminate("Missing arguments value for minGQR");
            if(std::stoi(args[i+1]) < 0) Terminate("Level for minGQR must be greater than 0");
            if(std::stod(args[i+2]) < 0 || std::stod(args[i+2]) > 1) Terminate("Rate for minGQR must be between 0 and 1");
            GQRlevel.push_back(std::stoi(args[++i]));
            GQRvalue.push_back(std::stod(args[++i]));
            RemovedQualityRate.push_back(0);
            continue;
        }
        if(args[i]=="--sample-stats"){
            CheckARG("sample-stats");
            IsSampleStats = true;
            continue;
        }
        if(args[i]=="--keep-multiallelic"){
            CheckARG("keep-multiallelic");
            IsRemoveMultiallelic = false;
            continue;
        }
        if(args[i]=="--ID"){
            CheckARG("ID");
            IsID = true;
            continue;
        }
        if(args[i]=="--remove"){
            CheckARG("remove");
            if(++i >= args.size()) Terminate("Missing argument value for remove");
            ReadRemoveList(args[i]);
            continue;
        }
        if(args[i]=="--verbose"){
            CheckARG("verbose");
            IsVerbose=true;
            continue;
        }
        if(args[i]=="--help"){
            ShowHelp();
            throw 0;
        }
        Terminate(args[i] + " is an invalid argument");
    }
    if(DefinedArguments.find("input")==DefinedArguments.end())
        Terminate("Missing input mode argument");
}
int LCVCFtools::StrToData(const std::string &String){
    if(String.find('.')!=std::string::npos) return -1;
    return std::stoi(String);
}
void LCVCFtools::CheckARG(std::string Argument, bool IsUnique){
    if(IsUnique)
        if(DefinedArguments.find(Argument)!=DefinedArguments.end())
            Terminate("Multiple " + Argument + " definition");
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
bool LCVCFtools::StringToVcf(const std::string& tmpLineString){
    std::vector<std::string> tmpStrings;
    boost::split(tmpStrings, tmpLineString, boost::algorithm::is_any_of("\t"));
    if(tmpStrings.size()<9) Terminate("Invalid number of columns");
    tmpSnpData.CHR=tmpStrings[0];
    tmpSnpData.POS=tmpStrings[1];
    tmpSnpData.ID=tmpStrings[2];
    tmpSnpData.REF=tmpStrings[3];
    tmpSnpData.ALT=tmpStrings[4];
    tmpSnpData.QUAL=tmpStrings[5];
    tmpSnpData.FILTER=tmpStrings[6];
    tmpSnpData.INFO=tmpStrings[7];
    tmpSnpData.FORMATstr=tmpStrings[8];
    if(lastFORMATstr!=tmpSnpData.FORMATstr){
        FORMATtagsMap.clear();
        FORMATtagsVector.clear();
        size_t i = 0;
        std::set<std::string> RequiredTags = {"GQ","DP","AD","GT","PL"};
        std::vector<std::string> tmpStrings2;
        boost::split(tmpStrings2, tmpSnpData.FORMATstr, boost::algorithm::is_any_of(":"));
        for(const std::string& String : tmpStrings2){
            FORMATtagsVector.push_back(String);
            FORMATtagsMap[String] = i++;
            RequiredTags.erase(String);
        }
        if(!RequiredTags.empty()){
            std::string tmpMsg = "Missing VCF required tag(s): ";
            for(std::string s : RequiredTags) tmpMsg += s + "; ";
            Terminate(tmpMsg);
        }
    }
    lastFORMATstr = tmpSnpData.FORMATstr;
    tmpStrings.erase(tmpStrings.begin(), tmpStrings.begin()+9);
    tmpSnpData.SampleDataVector.reserve(tmpStrings.size());
    size_t i = 0;
    for(const std::string& tmpString : tmpStrings){
        SampleDataStruct tmpSample;
        if(tmpString.empty()) return false;
        if(RemoveIndex.find(i++)!=RemoveIndex.end()) continue;
        boost::split(tmpSample.FORMAT,tmpString,boost::algorithm::is_any_of(":"));
        if(tmpSample.FORMAT.size()!=FORMATtagsVector.size()) Terminate("Incorrect number of fields for FORMAT");
        tmpSample.GT = &tmpSample.FORMAT[FORMATtagsMap["GT"]];
        tmpSample.PL = &tmpSample.FORMAT[FORMATtagsMap["PL"]];
        tmpSample.DP = &tmpSample.FORMAT[FORMATtagsMap["DP"]];
        tmpSample.AD = &tmpSample.FORMAT[FORMATtagsMap["AD"]];
        tmpSample.GQ = &tmpSample.FORMAT[FORMATtagsMap["GQ"]];
        tmpSnpData.SampleDataVector.push_back(std::move(tmpSample));
    }
    {
        std::vector<std::string> tmpStrings2;
        boost::split(tmpStrings2, tmpSnpData.ALT,boost::algorithm::is_any_of(","));
        for(size_t i(0); i < tmpStrings2.size()+1; i++)
            tmpSnpData.AlleleCountVector.push_back(std::pair<short,double>(i,0));
    }
    return true;
}
void LCVCFtools::ReadRemoveList(std::string Filename){
    if(Filename.empty()) return;
    std::ifstream File(Filename);
    if(!File.is_open()) Terminate("Remove list file does not exist or is not readable.");
    while(true){
        std::string tmpString;
        if(!std::getline(File, tmpString)) break;
        if(!tmpString.empty()) RemoveSamples.insert(tmpString);
    }
}
std::string LCVCFtools::GetParametersString(){
    std::string tmpString;
    std::ostringstream sso;
    sso << std::fixed;
    sso << std::setprecision(2);
    sso << "minGQ=" << minGQ << ";";
    sso << "minDP=" << minDP << ";";
    sso << "minGCR=" << minGCR << ";";
    sso << "MAF=" << MAF << ";";
    for(size_t i(0); i < DPRlevel.size(); i++) sso  << "minDPR" << i+1 << "=[" << DPRlevel[i] << ';' << DPRvalue[i] << "];";
    for(size_t i(0); i < GQRlevel.size(); i++) sso  << "minGQR" << i+1 << "=[" << GQRlevel[i] << ';' << GQRvalue[i] << "];";
    tmpString = sso.str();
    return tmpString;
}
void LCVCFtools::OpenInputStream(){
    //          Copyright Joe Coder 2004 - 2006.
    // Distributed under the Boost Software License, Version 1.0.
    //    (See accompanying file BOOST.license or copy at
    //          https://www.boost.org/LICENSE_1_0.txt)
    if(IsFile){
        if(IsGzipped){
            file.open(InputFilename, std::ios_base::in | std::ios_base::binary);
            if(!file.is_open()) Terminate("VCF file does not exist or is not readable");
            in.push(boost::iostreams::gzip_decompressor());
            in.push(file);
        }
        else{
            file.open(InputFilename);
            if(!file.is_open()) Terminate("VCF file does not exist or is not readable");
        }
    }
    else{
        if(IsGzipped){
            in.push(boost::iostreams::gzip_decompressor());
            in.push(std::cin);
        }
    }
}
void LCVCFtools::Terminate(std::string Msg){
    if(!Msg.empty()) Log("\033[1;31m**ERROR** " + Msg + "\033[0m");
    throw 1;
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
    while(true){
        std::string temp_string;
        if(!GetLine(temp_string))
            Terminate("Can't read VCF file, check the file format.");
        if(temp_string[0]=='#'){
            if(temp_string.substr(0,6)=="#CHROM"){
                std::vector<std::string> tmpHeaderStrings;
                boost::split(tmpHeaderStrings,temp_string,boost::algorithm::is_any_of("\t"));
                std::vector<std::string> tmpHeaderStringsCheck ({"#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"});
                for(size_t i(0); i<tmpHeaderStringsCheck.size();i++){
                    if(tmpHeaderStrings[i]!=tmpHeaderStringsCheck[i])
                        Terminate(tmpHeaderStrings[i] + " column must be " + tmpHeaderStringsCheck[i]);
                    HeaderColumns.push_back(tmpHeaderStrings[i]);
                }
                tmpHeaderStrings.erase(tmpHeaderStrings.begin(),tmpHeaderStrings.begin()+9);
                size_t i(0);
                Log(std::to_string(tmpHeaderStrings.size()) + " samples identified...");
                for(const std::string& tmpName : tmpHeaderStrings){
                    bool Remove=false;
                    if(!RemoveSamples.empty())
                        if(RemoveSamples.find(tmpName)!=RemoveSamples.end()){
                            RemoveIndex.insert(i);
                            Remove=true;
                        }
                    if(!Remove){
                        HeaderSamples.push_back(tmpName);
                        SampleStatsStruct tmpSTATS;
                        tmpSTATS.DP.resize(SampleStatsYLim+1,0);
                        tmpSTATS.GQ.resize(SampleStatsYLim+1,0);
                        SampleStatsVector.push_back(std::move(tmpSTATS));
                    }
                    i++;
                }
                if(RemoveSamples.size()>0)
                    Log(std::to_string(HeaderSamples.size()) + " samples remaining after 'remove' applied...");
                break;
            }
            CommentLines.push_back(temp_string);
        }
        else Terminate("VCF header without '#' starting character.");
    }
    if(!HeaderSamples.size()) Terminate("No samples in VCF file.");
    OutputHeader();
}
void LCVCFtools::ReadData(){
    size_t tmpCounter = 0;
    while(true){
        std::string tmpLineString;
        if(!GetLine(tmpLineString)) break;
        if(tmpLineString[0]=='#') Terminate("Invalid line at "+std::to_string(InputCounter)+".");
        tmpSnpData=SnpDataStruct();
        InputCounter++;
        if(!StringToVcf(tmpLineString)) Terminate("Failed to read data at line " + std::to_string(InputCounter)+", check the file format");
        if(Filter()){
            OutputCounter++;
            OutputLine();
        }
        if(IsVerbose && (++tmpCounter >= STATUS_FREQ)){
            tmpCounter = 0;
            ShowStatus();
        }
    }
}
void LCVCFtools::ShowStatus(){
    if(InputCounter==0) return;
    std::clog << NowString()
              << "Input=" << InputCounter
              << ";Output=" << OutputCounter << "(" << (static_cast<double>(OutputCounter)/InputCounter)*100 << "%);"
              << "Filtered:{";
    if(IsRemoveMultiallelic)
        std::clog << "MAL=" << (static_cast<double>(RemovedMultiallelic)/InputCounter)*100 << "%;";
    if(minGCR>0)
        std::clog << "GCR=" << (static_cast<double>(RemovedGenotypeCallRate)/InputCounter)*100 << "%;";
    if(MAF>0)
        std::clog << "MAF=" << (static_cast<double>(RemovedMAF)/InputCounter)*100 << "%;";
    for(size_t i(0);i < RemovedDepthRate.size(); i++)
        std::clog << "DPR[" << DPRlevel[i] << "," << DPRvalue[i] << "]="
                  << (static_cast<double>(RemovedDepthRate[i])/InputCounter)*100
                  << "%;";
    for(size_t i(0);i < RemovedQualityRate.size(); i++)
        std::clog << "GQR[" << GQRlevel[i] << "," << GQRvalue[i] << "]="
                  << (static_cast<double>(RemovedQualityRate[i])/InputCounter)*100
                  << "%;";
    std::clog << "}..."
              << std::endl;
}
void LCVCFtools::OutputSampleStatistics(){
    if(!IsSampleStats) return;
    Log("Calculating sample statistics...");
    for(size_t i(0); i < HeaderSamples.size(); i++) SampleStatsVector[i].Name = HeaderSamples[i];
    std::sort(SampleStatsVector.begin(),SampleStatsVector.end(),[](const SampleStatsStruct& a, const SampleStatsStruct& b)->bool{return a.NMR > b.NMR;});
    SampleStatsFile << "## Sample Statistics Table" << std::endl;
    SampleStatsFile << "## NMR=Mean non-missing rate; GCR=Mean genotype call rate; DP=Depth at a given level; GQ=Genotype at a given level" << std::endl;
    SampleStatsFile << "Sample\tVariable\tLevel\tValue" << std::endl;
    for(SampleStatsStruct& x : SampleStatsVector){
        int tmpDPsum(x.DP[0]), tmpGQsum(x.GQ[0]);
        for(int i(1); i < SampleStatsYLim+1; i++){
            tmpDPsum += x.DP[i];
            tmpGQsum += x.GQ[i];
            for(int j(i+1); j < SampleStatsYLim+1; j++){
                x.DP[i] += x.DP[j];
                x.GQ[i] += x.GQ[j];
            }
        }
        SampleStatsFile << x.Name << "\t" << "NMR" << "\t.\t" << x.NMR/OutputCounter  << std::endl;
        SampleStatsFile << x.Name << "\t" << "GCR" << "\t.\t" << x.GCR/OutputCounter  << std::endl;
        for(int i(1); i < SampleStatsYLim+1; i++)
            SampleStatsFile << x.Name << "\t" << "DP" << "\t" << i << "\t" << static_cast<double>(x.DP[i])/tmpDPsum  << std::endl;
        for(int i(1); i < SampleStatsYLim+1; i++)
            SampleStatsFile << x.Name << "\t" << "GQ" << "\t" << i << "\t" << static_cast<double>(x.GQ[i])/tmpGQsum  << std::endl;
    }
}
void LCVCFtools::OutputHeader(){
    for(const auto& tmpString : CommentLines) std::cout << tmpString << '\n';
    std::cout << "##LCVCFtools_v"
              << THIS_VERSION
              << " "
              << GetParametersString()
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
                 tmpSnpData.CHR    + '\t' +
                 tmpSnpData.POS    + '\t' +
                 tmpSnpData.ID     + '\t' +
                 tmpSnpData.REF    + '\t' +
                 tmpSnpData.ALT    + '\t' +
                 tmpSnpData.QUAL   + '\t' +
                 tmpSnpData.FILTER + '\t' +
                 tmpSnpData.INFO   + '\t';
    bool tmpBl = true;
    for(const auto& TAG : FORMATtagsVector){
        if(tmpBl) tmpBl = false;
        else tmpString += ':';
        tmpString += TAG;
    }
    for(const auto& SAMPLE : tmpSnpData.SampleDataVector){
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
bool LCVCFtools::CheckRate(const std::vector<int> &vec, int val, double qnt){
    double tmp(0);
    for(const auto& i : vec) if(i>=val) tmp++;
    tmp /= vec.size();
    return(tmp<qnt);
}
bool LCVCFtools::Filter(){
    if(IsRemoveMultiallelic && tmpSnpData.AlleleCountVector.size()>2) {
        RemovedMultiallelic++;
        return false;
    }
    std::vector<int> tmpGQ, tmpDP;
    tmpGQ.reserve(tmpSnpData.SampleDataVector.size());
    tmpDP.reserve(tmpSnpData.SampleDataVector.size());
    for(SampleDataStruct& Sample : tmpSnpData.SampleDataVector){
        int DP = StrToData(*Sample.DP);
        int GQ = StrToData(*Sample.GQ);
        tmpDP.push_back(DP);
        tmpGQ.push_back(GQ);
        if(DP==0){
            *Sample.GT = "./.";
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
            if(DP>=minDP && GQ>=minGQ) tmpSnpData.GCR++;
            std::vector<int> AD;
            int tmpADsum(0);
            std::vector<std::string> tmpADStrings;
            boost::split(tmpADStrings,*Sample.AD,boost::algorithm::is_any_of(","));
            for(std::string ADstr : tmpADStrings) {
                AD.push_back(StrToData(ADstr));
                tmpADsum += AD[AD.size()-1];
            }
            if(AD.size()!=tmpSnpData.AlleleCountVector.size())
                Terminate("Fatal error at AD.size()!=AlleleCount.size()");
            if(tmpADsum>0)
                for(size_t i(0); i<AD.size(); i++)
                    tmpSnpData.AlleleCountVector[i].second += static_cast<double>(AD[i])/tmpADsum;
        }
    }// for Sample END_HERE
    tmpSnpData.GCR /= tmpSnpData.SampleDataVector.size();
    /******* Apply minGCR FILTER *******/
    if(tmpSnpData.GCR < minGCR) {
        RemovedGenotypeCallRate++;
        return false;
    }
    /******* Apply minDPR FILTER *******/
    for(size_t i(0); i < DPRlevel.size();i++){
        if(CheckRate(tmpDP,DPRlevel[i],DPRvalue[i])){
            RemovedDepthRate[i]++;
            return false;
        }
    }
    /******* Apply minGQR FILTER *******/
    for(size_t i(0); i < GQRlevel.size();i++){
        if(CheckRate(tmpGQ,GQRlevel[i],GQRvalue[i])){
            RemovedQualityRate[i]++;
            return false;
        }
    }
    /******* Apply MAF FILTER *******/
    double AlleleSum(0);
    for(const auto& A : tmpSnpData.AlleleCountVector) AlleleSum += A.second;
    if(AlleleSum==0) return false;
    for(auto& A : tmpSnpData.AlleleCountVector) A.second /= AlleleSum;
    std::sort(tmpSnpData.AlleleCountVector.begin(),tmpSnpData.AlleleCountVector.end(),
              [](const std::pair<short,double>& a, const std::pair<short,double>& b)->bool{return a.second > b.second;});
    if((1-tmpSnpData.AlleleCountVector[0].second)<MAF){
        RemovedMAF++;
        return false;
    }
    /******* Apply ID *******/
    if(IsID) tmpSnpData.ID = tmpSnpData.CHR + ":" + tmpSnpData.POS + ':' + tmpSnpData.REF + ':' + tmpSnpData.ALT;
    /******* Update Sample Stats *******/
    if(IsSampleStats)
        for(size_t i(0); i<SampleStatsVector.size(); i++){
            if(tmpDP[i]>SampleStatsYLim) tmpDP[i] = SampleStatsYLim;
            if(tmpGQ[i]>SampleStatsYLim) tmpGQ[i] = SampleStatsYLim;
            SampleStatsVector[i].DP[tmpDP[i]]++;
            SampleStatsVector[i].GQ[tmpGQ[i]]++;
            SampleStatsVector[i].NMR += (tmpDP[i]>0)?1:0;
            SampleStatsVector[i].GCR += (tmpDP[i]>=minDP && tmpGQ[i]>=minGQ)?1:0;
        }
    return true;
}
void LCVCFtools::Run(){
    std::clog << std::fixed << std::setprecision(1);
    StartingTimeStr = NowString();
    Log("Running LCVCFtools v" + std::string(THIS_VERSION));
    Log("Running with parameters: " + GetParametersString());
    Log("Starting...");
    ReadHeader();
    ReadData();
    ShowStatus();
    OutputSampleStatistics();
    Log("Finished.");
}
