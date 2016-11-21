#include "Consolidator.cpp"
#include <iostream>
#include <string>


class ErrorFinderManager
{
 public:
       ErrorFinderManager();
       void performConsolidation(int argc, char *argv[]);
       void displayError(std::string argv);
       void initiateErrorFinder();
 private:
       int WINDOW, MIN_SNP,GAP, MA_SNP_END, TRUESNP;
       float MIN_CM,MA_ERR_THRESHOLD_START, MA_ERR_THRESHOLD_END,
             PCT_ERR_THRESHOLD,HO_THRESHOLD, TRUECM, PIELENGTH, MA_THRESHOLD,EMPIRICAL_MA_RESULT, EMPIRICAL_PIE_RESULT,CUT_VALUE,CM_CUT_VALUE;
       bool ISMOL, COUNTGAPERR;
       std::string BMIDFILE,BMATCHFILE,BSIDFILE,PEDFILE,
                     HPEDFILE,HMAPFILE,OPTION, HO_MISSING, LOGFILE;
       ErrorCalculator eCalculator;
       Consolidator consolidator;
       std::string wrongParam;
};


