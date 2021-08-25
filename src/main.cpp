#include "lcvcftools.h"
int main(int argc, char* argv[]){
    try{
        LCVCFtools a(argc, argv);
        a.Run();
    }
    catch(int i){return i;}
    return 0;
}
