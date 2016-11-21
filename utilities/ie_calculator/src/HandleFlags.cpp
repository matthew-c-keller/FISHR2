/*
 * HandleFlags.cpp
 *
 *  Created on: May 24, 2016
 *      Author: root
 */
#include <string.h>
#include <iostream>
#include <stdlib.h>     /* atoi */

#include "HandleFlags.hpp"

// Defintion for the static variable
//	http://stackoverflow.com/questions/37553581/assign-value-to-private-static-variable-in-a-class
std::string HandleFlags::pedfile;
std::string HandleFlags::bmidfile;
std::string HandleFlags::ibdfile;
std::string HandleFlags::outfile;
bool HandleFlags::reduced;
bool HandleFlags::silent;
double HandleFlags::trim;
int HandleFlags::windowsize;

HandleFlags::HandleFlags(){
}

int HandleFlags::setFlagValues(int argc,char *argv[]){
	bool pedflag = false;
	bool bmidflag = false;
	bool ibdflag = false;
	bool winflag = false;
	bool outflag = false;
	bool trimflag = false;
	reduced = false;
	silent = false;


		for (int i=1;i<argc;i++)
			{
				if (strcmp("-ped",argv[i]) == 0)
					{
					pedfile = argv[++i];
					pedflag = true;
					}
				else if (strcmp("-bmid",argv[i]) == 0)
					{
					bmidfile = argv[++i];
					bmidflag = true;
					}
				else if (strcmp("-ibd",argv[i]) == 0)
					{
					ibdfile = argv[++i];
					ibdflag = true;
					}
				else if (strcmp("-window",argv[i]) == 0 )
					{
					windowsize = atoi(argv[++i]);
					winflag = true;
					}
				else if (strcmp("-out",argv[i]) == 0)
					{
					outfile = argv[++i];
					outflag = true;
					}
				else if (strcmp("-emp_ma_threshold",argv[i]) == 0)//also called -trim
					{
					trim = atof(argv[++i]);
					trimflag = true;
					}
				else if (strcmp("-reduced",argv[i]) == 0)
					{
					reduced = true;
					}
				else if (strcmp("-silent",argv[i]) == 0)
					{
					silent = true;
					}
				else
					{
					std::cerr<<"Wrong parameter encountered:\t"<<argv[i]<<std::endl;
					return 1;
					}
			}

		if (!(pedflag	&&	bmidflag	&&	ibdflag  &&  winflag  &&  outflag	&&	trimflag))
		{
			std::cerr<<"flag missing.. Program is exiting";
			exit(0);
		}
		return 0;
	}

void HandleFlags::displayFlagValues()
{
		std::cout<<"pedfile "<<pedfile<<std::endl;
		std::cout<<"bmidfile "<<bmidfile<<std::endl;
		std::cout<<"ibdfile "<<ibdfile<<std::endl;
}


std::string HandleFlags::getpedfilename()
{
	return pedfile;
}

std::string HandleFlags::getbmidfilename()
{
	return bmidfile;
}

std::string HandleFlags::getibdfilename()
{
	return ibdfile;
}
std::string HandleFlags::getoutputfilename()
{
	return outfile;
}
double HandleFlags::gettrimvalue()
{
	return trim;
}

int HandleFlags::getwindowsize()
{
	return windowsize;
}

bool HandleFlags::getisreduced()
{
	return reduced;
}

bool HandleFlags::getissilent()
{
	return silent;
}

