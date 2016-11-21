#ifndef ERRORFINDERMANAGER_H
#define ERRORFINDERMANAGER_H
//#include "Consolidator.cpp"
#include "Consolidator.hpp"
#include "ErrorCalculator.hpp"
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
       int WINDOW, MIN_SNP,GAP, MA_SNP_END, TRUESNP,MIN_WINDOW,EXTENDSNP,SNIPEXTEND;//<piyush> SNIPEXTEND added  remove min_window and extendsnp
       float MIN_CM,MA_ERR_THRESHOLD_START, MA_ERR_THRESHOLD_END,
             PCT_ERR_THRESHOLD,HO_THRESHOLD, TRUECM, PIELENGTH, MA_THRESHOLD,EMPIRICAL_MA_RESULT, EMPIRICAL_PIE_RESULT;
       bool ISMOL, COUNTGAPERR;
       std::string BMIDFILE,BMATCHFILE,BSIDFILE,PEDFILE,
                     HPEDFILE,HMAPFILE,OPTION, HO_MISSING, LOGFILE, SNPWEIGHTFILE;
       ErrorCalculator eCalculator;
       Consolidator consolidator;
       std::string wrongParam;
};


#endif // ndef ERRORFINDERMANAGER_H
