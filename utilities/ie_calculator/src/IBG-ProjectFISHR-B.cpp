
//============================================================================
// Name        : IBG-FISHR_B.cpp
// Author      : Piyush
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "HandleFlags.hpp"
#include "ReadFiles.hpp"
#include "Compute.hpp"

void exitprogram(std::string);	//exitprogram handler

int main(int argc, char *argv[]) {

	HandleFlags flags;
	if (flags.setFlagValues(argc,argv) == 1){
		exitprogram("1");
		}
	ReadFiles readfiles;
	if (readfiles.readpedfile(flags.getpedfilename()) == 1){
		exitprogram("2");
		}

	if (readfiles.readibdfile(flags.getibdfilename()) == 1){
		exitprogram("3");
		}

	if (readfiles.readbmidfile(flags.getbmidfilename()) == 1){
		exitprogram("4");
		}

	Compute compute(flags,readfiles);	//	Opens file to write
	compute.calculate(flags,readfiles);


return 0;
}

void exitprogram(std::string status)
{
	std::cerr<<"Program is exiting with status: "<<status<<std::endl;
	exit(0);
}
