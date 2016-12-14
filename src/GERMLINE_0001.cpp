#include "GERMLINE.h"
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cstring>
#include <vector>
#include <stdio.h>      /* printf */
#include <stdlib.h>
#include "Haps.hpp"
#include "Compute.hpp"
#include "ErrorFinderManager.hpp"

double MIN_MATCH_LEN = 3;
int MARKER_SET_SIZE = 128;
bool PRINT_MATCH_HAPS = false;
bool ROI = false;
bool HAP_EXT = false;
bool WIN_EXT = false;
bool ALLOW_HOM = false;
bool HOM_ONLY = false;
bool HAPLOID = false;
bool SILENT = false;
bool DEBUG = false;
bool BINARY_OUT = false;
int MAX_ERR_HOM = 4;
int MAX_ERR_HET = 1;
int ERR_W=0;
bool NO_SUFFIX=false;
bool REDUCE=false;
bool GERMLINE_OUTPUT = false;
std::string PEDFILE ="";
std::string OUTFILE="";
std::string MAPFILE="";

std::string HAPFILE = "";
std::string SAMPLEFILE="";
std::string GENFILE="";

double IBD_THRESHOLD = -1.0;

bool isInteger(std::string s);
bool isPath(std::string s);
bool IBD = false; // To check if the file is IBD2 or IBD4. Changes done in Match,cpp file to compute .match along with .bmatch
void helpShowParameters();

int main(int argc, char* argv[])
{
//http://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
/*	char** argvcopy = new char*[argc];
	for (int i = 0 ; i < argc ; i++)
	{
		argvcopy[i] = new char[strlen(argv[i])+1];
		strcpy(argvcopy[i],argv[i]);
	}*/
	std::map <std::string, std::string> grade_list;
	//grade_list["-min_m"] = "NULL";

	grade_list["-min_cm_initial"] = "NULL";
	grade_list["-err_hom"] = "NULL";
	grade_list["-err_het"] = "NULL";
	grade_list["-from_snp"] = "NULL";
	grade_list["-to_snp"] = "NULL";
	grade_list["-print"] = "NULL";
	grade_list["-silent"] = "NULL";
	grade_list["-debug"] = "NULL";
	grade_list["-map"] = "NULL";
	grade_list["-pedfile"] = "NULL";
	grade_list["-mapfile"] = "NULL";
	grade_list["-bits"] = "NULL";
	grade_list["-err_w"] = "NULL";
	grade_list["-homoz-only"] = "NULL";
	grade_list["-homoz"] = "NULL";
	//grade_list["-outfile"] = "NULL";
	//grade_list["-bin_out"] = "NULL";
	grade_list["-haploid"] = "NULL";
	grade_list["-h_extend"] = "NULL";
	grade_list["-w_extend"] = "NULL";
	grade_list["-no_suffix"] = "NULL";
	//grade_list["-reduced"] = "NULL";
	grade_list["-hapsfile"] = "NULL";
	grade_list["-geneticmapfile"] = "NULL";
	grade_list["-samplefile"] = "NULL";
	grade_list["-germline_output"] = "NULL";
	grade_list["-ibd2"] = "NULL";

	//-mapfile ./src/Test.map  -pedfile ./src/Test.ped -outfile BEAGLE_OUT -bin_out -bits 20 -err_hom 0 -err_het 0 -min_m 3 -homoz  -w_extend -h_extend
	std::map <std::string, std::string> grade_list2;
	grade_list2["-extendSNP"] = "NULL";
//	grade_list2["-bmatch"] = "NULL";
//	grade_list2["-bmid"] = "NULL";
//	grade_list2["-bsid"] = "NULL";
//	grade_list2["-reduced"] = "NULL";	//never use it, it's replaced by min_snp and min_cm
//	grade_list2["-ped-file"] = "NULL";
	grade_list2["-min_cm_final"] = "NULL";
	grade_list2["-min_snp"] = "NULL";
	grade_list2["-holdout-ped"] = "NULL";
	grade_list2["-holdout-map"] = "NULL";
	grade_list2["-window"] = "NULL";
	grade_list2["-holdout-threshold"] = "NULL";
	grade_list2["-trueCM"] = "NULL";
	grade_list2["-trueSNP"] = "NULL";
	grade_list2["-holdout-missing"] = "NULL";
	grade_list2["-gap"] = "NULL";
	grade_list2["-ma-snp"] = "NULL";
	grade_list2["-pct-err-threshold"] = "NULL";
	grade_list2["-emp-pie-threshold"] = "NULL";
	grade_list2["-output-type"] = "NULL";
	grade_list2["-snpfile"] = "NULL";
	grade_list2["-log-file"] = "NULL";
	grade_list2["-ma-threshold"] = "NULL";
	grade_list2["-empirical-ma-threshold"] = "NULL";
	grade_list2["-PIE.dist.length"] = "NULL";
	grade_list2["-count.gap.errors"] = "NULL";
	bool help = false;
	bool bad_param = false;
	//bool ibd = false;	//Defined globally s othat it can be used inthe Match.cpp file

	std::vector<std::string> flagPool;
	for(map<std::string, std::string>::iterator it = grade_list.begin(); it != grade_list.end(); ++it) {
		flagPool.push_back(it->first);
		//std::cout<<it->first<<endl;
	}
	for(map<std::string, std::string>::iterator it = grade_list2.begin(); it != grade_list2.end(); ++it) {
		flagPool.push_back(it->first);
		//std::cout<<it->first<<endl;
	}



	for (size_t i = 0 ; i < argc  ; i++){
		string temp = argv[i];
		if (argc == 1)
		{
				cerr<<endl;
				cerr<<"***********************************************************************************"<<endl;
				cerr<<"*                        FISHR2 (Genetic Analysis Program)                        *"<<endl;
				cerr<<"*                                  				                                 *"<<endl;
				cerr<<"*                (C) 2014  Matthew C. Keller, Doug Bjelland, Piyush Sudip Patel   *"<<endl;
				cerr<<"*                    Institute for Behavioral Genetics (IBG)                      *"<<endl;
				cerr<<"*                       University of Colorado at Boulder                         *"<<endl;
				cerr<<"***********************************************************************************"<<endl;
				cerr<<"****************************\"./FISHR2 -help\" for Usage *****************************"<<endl;
				cerr<<"******For a complete worked through example, see \"FISHR2.complete.example.R\"*****"<<endl;
				cerr<<endl;
			    exit(0);
		}
		if( strcmp(argv[i], "-help") == 0 )
		{
			help = true;
			helpShowParameters();
			//bad_param = true;
			//continue;
		}
/*		if( strcmp(argv[i], "-ibd") == 0 )
		{
			IBD = true;
			continue;
		}*/

		if (temp[0] != '-')
			 continue;


		bool flag = false;
		for (size_t j = 0; j < flagPool.size(); j++ )
			{
				//cout<<temp<<"\t"<<flagPool[j]<<endl;
				if (flagPool[j] == argv[i])
				{
					flag = true;
					break;
				}

			}
		if (flag == false)
		{

			cerr<<"ERROR: Wrong flag. "<<argv[i]<<endl;
			exit(0);
		}
	}

	char** glcopy = new char*[argc];
	int glcopycount = 0;
	for (int i = 0; i < argc ; i++)
	{
		//cout<<argv[i]<<endl;
		//std::cout<<"glcopycount = "<<glcopycount<<std::endl;
		if (i == 0)
		{
			glcopy[i] = new char[strlen(argv[i])+1];	strcpy(glcopy[i],argv[i]);	glcopycount++;
			continue;
		}
		  if (grade_list.count(argv[i]))
		  {
			  //std::cout<<argv[i]<<std::endl;
			  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  if( strncmp(argv[i], "-min_cm_initial", strlen("-min_cm_initial")) == 0)
			  {
				  i++;
				  //std::cout<<argv[i]<<std::endl;
				  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if( strncmp(argv[i], "-germline_output", strlen("-germline_output")) == 0)
			  {
				  i++;
				  //std::cout<<argv[i]<<std::endl;
				  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if( strncmp(argv[i], "-ibd2", strlen("-ibd2")) == 0)
			  {
				  i++;
				  //std::cout<<argv[i]<<std::endl;
				  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if( strncmp(argv[i], "-err_hom", strlen("-err_hom")) == 0)
			  {
				  i++;
				  //std::cout<<argv[i]<<std::endl;
				  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if( strncmp(argv[i], "-err_het", strlen("-err_het")) == 0 )
			  {
				  i++;
				  //std::cout<<argv[i]<<std::endl;
				  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if( strncmp(argv[i], "-from_snp", strlen("-from_snp")) == 0 )
			  {	  i++;
			  	  //std::cout<<argv[i]<<std::endl;
			  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }

			  else if( strncmp(argv[i], "-to_snp", strlen("-to_snp")) == 0)
			  {i++;
			  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if( strcmp(argv[i], "-map") == 0 )
			  {i++;
			  //std::cout<<argv[i]<<std::endl;
			  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if( strncmp(argv[i], "-pedfile", strlen("-pedfile"))==0)
			  {i++;
			  //std::cout<<argv[i]<<std::endl;
			  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if( strncmp(argv[i], "-mapfile" ,strlen("-mapfile")) == 0)
			  {i++;
			  //std::cout<<argv[i]<<std::endl;
			  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if( strncmp(argv[i], "-bits", strlen("-bits")) == 0)
			  {i++;
			  //std::cout<<argv[i]<<std::endl;
			  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if( strncmp(argv[i], "-err_w", strlen("-err_w")) == 0)
			  {i++;
			  //std::cout<<argv[i]<<std::endl;
			  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if (strncmp(argv[i], "-hapsfile", strlen("-hapsfile")) == 0)
			  {i++;
			 // std::cout<<argv[i]<<std::endl;
			  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;
			  }
			  else if (strncmp(argv[i], "-geneticmapfile", strlen("-geneticmapfile")) == 0)
			  {i++;
			  //std::cout<<argv[i]<<std::endl;
			  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;			  }
			  else if (strncmp(argv[i], "-samplefile", strlen("-samplefile")) == 0)
			  {i++;
			  //std::cout<<argv[i]<<std::endl;
			  glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;			  }
		  }
		  else
			  {
			  	  //glcopy[i] = new char[strlen(" ")+1];	strcpy(glcopy[i]," ");	glcopycount++;
			  	  continue;
			  	  //std::cout<<"NOT FOUND"<<"\t"<<argv[i]<<std::endl;
			  }
	}
/*Add the outfile flaghere and later in the program, make sure the OUTFILE variable is handled well with
	 * either the hapfile name or the pedfilename
	 *
	 * */
	glcopy[glcopycount ] = new char[strlen("-outfile")+1];
	strcpy(glcopy[glcopycount],"-outfile");
	glcopycount++;

	glcopy[glcopycount ] = new char[strlen("-bin_out")+1];
	strcpy(glcopy[glcopycount],"-bin_out");
	glcopycount++;

	glcopy[glcopycount ] = new char[strlen("-reduced")+1];
	strcpy(glcopy[glcopycount],"-reduced");
	glcopycount++;

/*	for (int i = 0; i < argc ; i++)
	{
		if (strcmp(argv[i],"-help")==0)
		{
			glcopy[glcopycount ] = new char[strlen("-help")+1];
			strcpy(glcopy[glcopycount],"-help");
			glcopycount++;
		}
	}*/



//exit(0);


/*
std::cout<<"here1"<<std::endl;

	for (int i = 0;	i < argc ; i++ )
	{
		std::cout<<argv[i]<<std::endl;
	}

std::cout<<"--------"<<std::endl;


std::cout<<glcopycount<<std::endl;
	for (int i = 0;	i < glcopycount ; i++ )
	{
		std::cout<<"i="<<i<<"\t"<<glcopy[i]<<std::endl;
	}

exit(0);
*/




/****************************FISHR*****************************/
/*
	-pedfile ./src/Beagle.Phased.Group2.1k.ped
	-bmatch ./src/Beagle.Phased.Group2.1k.bmatch
	-bsid ./src/Beagle.Phased.Group2.1k.bsid
	-bmid ./src/Beagle.Phased.Group2.1k.bmid
	-reduced 64 3 -window 50 -gap 100 	-output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -ma-threshold 0.2 -log-file logs
*/

/*	std::map <std::string, std::string> grade_list2;
	grade_list2["-extendSNP"] = "NULL";
//	grade_list2["-bmatch"] = "NULL";
//	grade_list2["-bmid"] = "NULL";
//	grade_list2["-bsid"] = "NULL";
//	grade_list2["-reduced"] = "NULL";	//never use it, it's replaced by min_snp and min_cm
//	grade_list2["-ped-file"] = "NULL";

	grade_list2["-min_cm_final"] = "NULL";
	grade_list2["-min_snp"] = "NULL";
	grade_list2["-holdout-ped"] = "NULL";
	grade_list2["-holdout-map"] = "NULL";
	grade_list2["-window"] = "NULL";
	grade_list2["-holdout-threshold"] = "NULL";
	grade_list2["-trueCM"] = "NULL";
	grade_list2["-trueSNP"] = "NULL";
	grade_list2["-holdout-missing"] = "NULL";
	grade_list2["-gap"] = "NULL";
	grade_list2["-ma-snp"] = "NULL";
	grade_list2["-pct-err-threshold"] = "NULL";
	grade_list2["-emp-pie-threshold"] = "NULL";
	grade_list2["-output-type"] = "NULL";
	grade_list2["-snpfile"] = "NULL";
	grade_list2["-log-file"] = "NULL";
	grade_list2["-ma-threshold"] = "NULL";
	grade_list2["-empirical-ma-threshold"] = "NULL";
	grade_list2["-PIE.dist.length"] = "NULL";
	grade_list2["-count.gap.errors"] = "NULL";*/


	char** fishrcopy = new char*[argc];
	int fishrcopycount = 0;
	int lowramflag = 0;//0 for normal version of fishr, 1 for low ram

	//glcopy[glcopycount] = new char[strlen(argv[i])+1];	strcpy(glcopy[glcopycount],argv[i]);	glcopycount++;

for (int i = 0; i < argc ; i++)
	{

		if ((strcmp(argv[i],"-low_ram")==0))
			{
				lowramflag = 1;
			}

	if (i==0)
		{
			fishrcopy[fishrcopycount] = new char[strlen(argv[i])+1];	strcpy(fishrcopy[fishrcopycount],argv[i]);	fishrcopycount++;
		}


	if (grade_list2.count(argv[i]))
	{
		//std::cout<<argv[i]<<std::endl;
		fishrcopy[fishrcopycount] = new char[strlen(argv[i])+1];	strcpy(fishrcopy[fishrcopycount],argv[i]);	fishrcopycount++;

	  if (	(strcmp(argv[i],"-extendSNP")==0) || (strcmp(argv[i],"-bmatch")==0) || (strcmp(argv[i],"-bmid")==0) ||
			(strcmp(argv[i],"-bsid")==0) || (strcmp(argv[i],"-ped-file")==0) || (strcmp(argv[i],"-holdout-ped")==0)	||
			(strcmp(argv[i],"-holdout-map")==0) || (strcmp(argv[i],"-window")==0) || (strcmp(argv[i],"-holdout-threshold")==0) ||
			(strcmp(argv[i],"-trueCM")==0) || ( strcmp( argv[i],"-trueSNP" )==0) || (strcmp(argv[i],"-holdout-missing")==0) ||
			(strcmp(argv[i],"-gap")==0) || (strcmp(argv[i],"-ma-snp")==0) || (strcmp(argv[i],"-pct-err-threshold")==0) ||
			(strcmp(argv[i],"-emp-pie-threshold")==0) || (strcmp(argv[i],"-output-type")==0) || (strcmp(argv[i], "-snpfile") ==0) ||
			(strcmp(argv[i],"-log-file")==0) || (strcmp(argv[i],"-ma-threshold")==0) ||(strcmp(argv[i],"-empirical-ma-threshold")==0 ) ||
			(strcmp(argv[i],"-PIE.dist.length")==0) || (strcmp(argv[i],"-count.gap.errors")==0 ) ||  (strcmp(argv[i],"-min_cm_final")==0  ) ||
			(strcmp(argv[i],"-min_snp")==0))
				  {
					  i++;
					  //std::cout<<argv[i]<<std::endl;
					  //fishrcopy[i] = new char[strlen(argv[i])+1];	strcpy(fishrcopy[i],argv[i]);	fishrcopycount++;
					  fishrcopy[fishrcopycount] = new char[strlen(argv[i])+1];	strcpy(fishrcopy[fishrcopycount],argv[i]);	fishrcopycount++;

				  }
	}
	else
	{
		//fishrcopy[i] = new char[strlen(argv[i])+1];	strcpy(fishrcopy[i]," ");	fishrcopycount++;
		continue;
	}
}
//exit(0);


/*for (int i = 0;	i < argc ; i++ )
{
	std::cout<<argv[i]<<std::endl;
}
std::cout<<"--------------------"<<std::endl;
for (int i = 0;	i < fishrcopycount ; i++ )
{
	std::cout<<fishrcopy[i]<<std::endl;
}*/


/*Before the execution of Fishr begins Add the bmid sid bmatch and pedfile. preferably handle this just before the execution of fishr*/
//exit(0);






/*	// parse arguments
	string rs_range[2] , map; map = rs_range[0] = rs_range[1] = "";
	string params = argv[0];
        string badparam = "";
	bool bad_param = false;*/
// parse arguments
string rs_range[2] , map; map = rs_range[0] = rs_range[1] = "";
string params = glcopy[0];
string badparam = "";
string germline_output = "./";

	for(int i=1;i<glcopycount;i++){
		//std::cout<<glcopy[i]<<std::endl;
		params += " " + string(glcopy[i]);
		if( strncmp(glcopy[i], "-min_cm_initial", strlen("-min_cm_initial")) == 0 && i < glcopycount-1)				MIN_MATCH_LEN = atof(glcopy[++i]);
		else if( strncmp(glcopy[i], "-err_hom", strlen("-max_err")) == 0 && i < glcopycount-1)		{ MAX_ERR_HOM = atoi(glcopy[++i]); }
		else if( strncmp(glcopy[i], "-err_het", strlen("-err_het")) == 0 && i < glcopycount-1)		{ MAX_ERR_HET = atoi(glcopy[++i]); }
		else if( strncmp(glcopy[i], "-from_snp", strlen("-from_snp")) == 0 && i < glcopycount-1 )	rs_range[0] = glcopy[++i];
		else if( strncmp(glcopy[i], "-to_snp", strlen("-to_snp")) == 0 && i < glcopycount-1 )		rs_range[1] = glcopy[++i];
		else if( strncmp(glcopy[i], "-print", strlen("-print")) == 0 )						PRINT_MATCH_HAPS = true;
		else if( strncmp(glcopy[i], "-silent", strlen("-silent")) == 0 )						SILENT = true;
		else if( strncmp(glcopy[i], "-debug", strlen("-debug")) == 0 )						DEBUG = true;
		else if( strcmp(glcopy[i], "-map") == 0 && i < glcopycount-1)				map = glcopy[++i];
		else if( strncmp(glcopy[i], "-pedfile", strlen("-pedfile")) == 0 && i < glcopycount-1)                           PEDFILE = glcopy[++i];
		else if( strcmp(glcopy[i], "-mapfile") == 0 && i < glcopycount-1)	MAPFILE = glcopy[++i];
		else if( strncmp(glcopy[i], "-outfile",strlen("-outfile")) == 0 && i < glcopycount-1) {}//OUTFILE = glcopy[++i];
		else if( strncmp(glcopy[i], "-bits", strlen("-bits")) == 0 && i < glcopycount-1)				MARKER_SET_SIZE = atoi(glcopy[++i]);

		else if( strncmp(glcopy[i], "-err_w", strlen("-err_w")) == 0 && i < glcopycount-1)                  ERR_W = atoi(glcopy[++i]);

		else if( strncmp(glcopy[i], "-homoz-only", strlen("-homoz-only")) == 0 )				{ ALLOW_HOM = true; HOM_ONLY = true; }
		else if( strncmp(glcopy[i], "-homoz", strlen("-homoz")) == 0 )						ALLOW_HOM = true;
		else if( strncmp(glcopy[i], "-bin_out", strlen("-bin_out")) == 0 )						BINARY_OUT = true;
		else if( strncmp(glcopy[i], "-haploid", strlen("-haploid")) == 0 )					{ HAPLOID = true; HAP_EXT = true; }
		else if( strncmp(glcopy[i], "-h_extend", strlen("-h_extend")) == 0 )					HAP_EXT = true;
		else if( strncmp(glcopy[i], "-w_extend", strlen("-w_extend")) == 0 )					WIN_EXT = true;
		else if( strncmp(glcopy[i], "-no_suffix", strlen("-no_suffix")) == 0 )                                      NO_SUFFIX = true;
		else if( strncmp(glcopy[i], "-reduced", strlen("-reduced")) == 0 )                                      REDUCE = true;
		else if ( strncmp(glcopy[i], "-germline_output", strlen("-germline_output")) == 0 )
		{
			GERMLINE_OUTPUT = true;
			germline_output = glcopy[++i];
		}
		else if ( strncmp(glcopy[i], "-ibd2", strlen("-ibd2")) == 0 )
		{
			IBD = true;
			IBD_THRESHOLD = atof(glcopy[++i]);
		}

		else if (strncmp(glcopy[i], "-hapsfile", strlen("-hapsfile")) == 0)			HAPFILE = glcopy[++i];
		else if (strncmp(glcopy[i], "-geneticmapfile", strlen("-geneticmapfile")) == 0)			GENFILE = glcopy[++i];
		else if (strncmp(glcopy[i], "-samplefile", strlen("-samplefile")) == 0)	SAMPLEFILE = glcopy[++i];
		else
                {
					std::cerr<<glcopy[i]<<std::endl;
                   badparam = badparam + " " + glcopy[ i ];
                   bad_param = true;

                }
	}
if (IBD && (IBD_THRESHOLD == -1.0))
{
	bad_param = true;
	std::cerr<<"-ibd2 <value>";
}
	//BINARY_OUT = true;	//Should always be true



/*	for(int i=1;i<argc;i++){
		//std::cout<<argv[i]<<std::endl;
		params += " " + string(argv[i]);
		if( strncmp(argv[i], "-min_m", strlen("-min_m")) == 0 && i < argc-1)				MIN_MATCH_LEN = atof(argv[++i]);
		else if( strncmp(argv[i], "-err_hom", strlen("-max_err")) == 0 && i < argc-1)		{ MAX_ERR_HOM = atoi(argv[++i]); }
		else if( strncmp(argv[i], "-err_het", strlen("-max_err")) == 0 && i < argc-1)		{ MAX_ERR_HET = atoi(argv[++i]); }
		else if( strncmp(argv[i], "-from_snp", strlen("-from_snp")) == 0 && i < argc-1 )	rs_range[0] = argv[++i];
		else if( strncmp(argv[i], "-to_snp", strlen("-to_snp")) == 0 && i < argc-1 )		rs_range[1] = argv[++i];
		else if( strncmp(argv[i], "-print", strlen("-print")) == 0 )						PRINT_MATCH_HAPS = true;
		else if( strncmp(argv[i], "-silent", strlen("-silent")) == 0 )						SILENT = true;
		else if( strncmp(argv[i], "-debug", strlen("-debug")) == 0 )						DEBUG = true;
		else if( strcmp(argv[i], "-map") == 0 && i < argc-1)				map = argv[++i];
		else if( strncmp(argv[i], "-pedfile", strlen("-pedfile")) == 0 && i < argc-1)                           PEDFILE = argv[++i];
		else if( strcmp(argv[i], "-mapfile") == 0 && i < argc-1)	MAPFILE = argv[++i];
		else if( strncmp(argv[i], "-outfile",strlen("-outfile")) == 0 && i < argc-1) OUTFILE = argv[++i];
		else if( strncmp(argv[i], "-bits", strlen("-bits")) == 0 && i < argc-1)				MARKER_SET_SIZE = atoi(argv[++i]);

		else if( strncmp(argv[i], "-err_w", strlen("-err_w")) == 0 && i < argc-1)                  ERR_W = atoi(argv[++i]);

		else if( strncmp(argv[i], "-homoz-only", strlen("-homoz-only")) == 0 )				{ ALLOW_HOM = true; HOM_ONLY = true; }
		else if( strncmp(argv[i], "-homoz", strlen("-homoz")) == 0 )						ALLOW_HOM = true;
		else if( strncmp(argv[i], "-bin_out", strlen("-bin_out")) == 0 )						BINARY_OUT = true;
		else if( strncmp(argv[i], "-haploid", strlen("-haploid")) == 0 )					{ HAPLOID = true; HAP_EXT = true; }
		else if( strncmp(argv[i], "-h_extend", strlen("-h_extend")) == 0 )					HAP_EXT = true;
		else if( strncmp(argv[i], "-w_extend", strlen("-w_extend")) == 0 )					WIN_EXT = true;
		else if( strncmp(argv[i], "-no_suffix", strlen("-no_suffix")) == 0 )                                      NO_SUFFIX = true;
		else if( strncmp(argv[i], "-reduced", strlen("-reduced")) == 0 )                                      REDUCE = true;

		else if (strncmp(argv[i], "-hapsfile", strlen("-hapsfile")) == 0)			HAPFILE = argv[++i];
		else if (strncmp(argv[i], "-geneticmapfile", strlen("-geneticmapfile")) == 0)			GENFILE = argv[++i];
		else if (strncmp(argv[i], "-samplefile", strlen("-samplefile")) == 0)	SAMPLEFILE = argv[++i];
		else
                {
                   badparam = badparam + " " + argv[ i ];
                   bad_param = true;
                }
	}
	*/

//std::cout<<params<<std::endl;
	//exit(0);

/*
	if (PEDFILE!="" && MAPFILE!="")
	{

	}
*/

/*

	else if (HAPFILE!="" &&  MAPFILE!="" && SAMPLEFILE!=""  )	//if hapfiel and mapfile is provided
	{
		OUTFILE =  HAPFILE.substr(0,HAPFILE.find_last_of("."));
		Compute compute(HAPFILE, SAMPLEFILE, GENFILE,true );
		PEDFILE = HAPFILE.substr(0,HAPFILE.find_last_of("."))+".ped";
	}


	else	if (HAPFILE == "" && GENFILE == ""	 && SAMPLEFILE == "")
	{

	}
	else
	{
	else if ( !(HAPFILE == "" || GENFILE == ""	 || SAMPLEFILE == "" || PEDFILE!= "" || MAPFILE!="" ))
	{
		std::cerr<<" If providing input file in .hap format, must include .hap and .sample file. If providing .ped file, must also include .map file .Program exiting"<<std::endl;
		exit(0);
	}

	else
	{
		OUTFILE =  HAPFILE.substr(0,HAPFILE.find_last_of("."));
		Compute compute(HAPFILE, SAMPLEFILE, GENFILE );
		PEDFILE = HAPFILE.substr(0,HAPFILE.find_last_of("."))+".ped";
		MAPFILE = HAPFILE.substr(0,HAPFILE.find_last_of("."))+".map";
	}






	//}


	if (HAPFILE == "" && GENFILE == ""	 && SAMPLEFILE == "" && PEDFILE== "" && MAPFILE=="" )
	{
		std::cerr<<"No Input file provided.. Exiting"<<std::endl;
		exit(0);
	}
*/
bool input = false; // (HAPFILE,GENFILE,SAMPLEFILE) or (HAPFILE, MAPFILE) or (PEDFILE, MAPFILE)


/*
if (HAPFILE == "" && GENFILE == ""	 && SAMPLEFILE == "")
{


}
else
{
	if (HAPFILE == "" || GENFILE == ""	 || SAMPLEFILE == "" || PEDFILE!= "" || MAPFILE!="" )
	{
		std::cerr<<" If providing input file in .hap format, must include .hap, .gen and .sample file. If providing .ped file, must also include .map file .Program exiting"<<std::endl;
		exit(0);
	}

	else
	{
		OUTFILE =  HAPFILE.substr(0,HAPFILE.find_last_of("."));
		Compute compute(HAPFILE, SAMPLEFILE, GENFILE );
		PEDFILE = HAPFILE.substr(0,HAPFILE.find_last_of("."))+".ped";
		MAPFILE = HAPFILE.substr(0,HAPFILE.find_last_of("."))+".map";
	}
}
*/








if (HAPFILE == "" && GENFILE == ""	 && SAMPLEFILE == "" && PEDFILE== "" && MAPFILE=="" )
{
	std::cerr<<"No Input file provided.. Exiting"<<std::endl;
	exit(0);
}

if (HAPFILE!="" && MAPFILE!="" && SAMPLEFILE!="")
{
	OUTFILE =  HAPFILE.substr(0,HAPFILE.find_last_of("."));
	Compute compute(HAPFILE, SAMPLEFILE );
	PEDFILE = HAPFILE.substr(0,HAPFILE.find_last_of("."))+".ped";

}
else if (HAPFILE!="" && SAMPLEFILE!="" && GENFILE!="")
{
	OUTFILE =  HAPFILE.substr(0,HAPFILE.find_last_of("."));
	Compute compute(HAPFILE, SAMPLEFILE, GENFILE );
	PEDFILE = HAPFILE.substr(0,HAPFILE.find_last_of("."))+".ped";
	MAPFILE = HAPFILE.substr(0,HAPFILE.find_last_of("."))+".map";

}
else
{

}













if (OUTFILE == "")
{

	OUTFILE =  PEDFILE.substr(0,PEDFILE.find_last_of("."));
	//std::cout<<OUTFILE<<std::endl;
	//exit(0);
}


	if (OUTFILE == "")
	{

		OUTFILE =  PEDFILE.substr(0,PEDFILE.find_last_of("."));
		//std::cout<<OUTFILE<<std::endl;
		//exit(0);
	}


	if(MIN_MATCH_LEN < 0)
	{
		cerr << "-min_cm_initial must be non-negative" << endl << endl;
		bad_param = true;
	} else if(MAX_ERR_HOM < 0 || MAX_ERR_HET < 0 )
	{
		cerr << "-err_hom,-err_het must be non-negative" << endl << endl;
		bad_param = true;
	}
//	 else if(MAX_ERR_HOM == 0 || MAX_ERR_HET == 0 )
       /* {
                MAX_ERR_HET=2; MAX_ERR_HOM=2;
        }*/

	if(bad_param)
	{    
		if (!help)
		{

		cerr << badparam<<endl 
               <<" above options which you have provided are not allowed"<<endl
                      << "usage: " << glcopy[0]
                << " -pedfile [ped file ] -mapfile [map file ]"
                << " -outfile [ out file ]"<<endl
                << "<flags (optional)>" << endl;
		}
		helpShowParameters();
		/*cerr<< "flags:" << endl
		<< '\t' << "-silent" << '\t' << "Suppress all output except for warnings and prompts." << endl
		<< '\t' << "-bin_out" << '\t' << "Output in binary format to save space." << endl
		<< '\t' << "-min_cm_initial" << '\t' << "Minimum length for match to be used for imputation (in cM or MB)." << endl
		<< '\t' << "-err_hom" << '\t' << "Maximum number of mismatching homozygous markers (per slice)." << endl
		<< '\t' << "-err_het" << '\t' << "Maximum number of mismatching heterozygous markers (per slice)." << endl
		<< '\t' << "-from_snp" << '\t' << "Start SNP (rsID)." << endl
		<< '\t' << "-to_snp" << '\t' << "End SNP (rsID)." << endl
		<< '\t' << "-haps" << '\t' << "Print the resolved haplotypes in a seperate HAPS file." << endl
		<< '\t' << "-map" << '\t' << "Genetic distance map." << endl
		<< '\t' << "-bits" << '\t' << "Slice size." << endl
		<< '\t' << "-homoz" << '\t' << "Allow self matches (homozygosity)" << endl
		<< '\t' << "-homoz-only" << '\t' << "Look for autozygous/homozygous segments only, does not detect IBD" << endl
		<< '\t' << "-haploid" << '\t' << "Treat input individual as two fully phased chromosomes with no recombination\n\toutput IDs with 0/1 suffix for chromosome destinction" << endl
		<< '\t' << "-h_extend" << '\t' << "Extend from seeds if *haplotypes* match" << endl
		<< '\t' << "-w_extend" << '\t' << "Extend, one marker at a time, beyong the boundaries of a found match" << endl
		<< '\t' << "-reduced" << '\t' << "output only reduced elements" << endl
		<< '\t' << "-no_suffix" << '\t' << "use with -hapoloid to outputthe cell with no .0 or .1 suffix" << endl
		<< '\t' << "--err_w" << '\t' << "use with -wextend to allow error matches like 1,2,3 etc in extending" << endl
		<<'\t'	<< "-hap"	 <<'\t'	 << "if using hap file, use this flag. usuage: -hap <hapfilename.hap>  Must also supply \".sample\" and \".gen\" file using flags -sample and -gen"<< endl
		<<'\t'	<< "-sample"	 <<'\t'	 << "if using hap file, use this flag for the corresponding sample file. usuage: -sample <samplefilename.sample>  Must also supply \".hap\" and \".gen\" file using flags -hap and -gen"<< endl
		<<'\t'	<< "-gen"	 <<'\t'	 << "if using hap file, use this flag to supply the corresponding gen file. usuage: -gen <genfile.gen>  Must also supply \".hap\" and \".sample\" file using flags -hap and -sample"<< endl;
		*/
		return 0;
	}
	if( rs_range[0] != "" && rs_range[1] != "" )
	{
		ROI = true;
		ALL_SNPS.setROI(rs_range);
	}

	if(map != "")
	{
		ALL_SNPS.loadGeneticDistanceMap( map );
	}

	GERMLINE germline;
    germline.mine( params );


    /*Add these t o the command line argument for fishr,since these are must to be able to run fishr*/
    std::string fishrpedfile = PEDFILE;
    std::string fishrbsidfile = OUTFILE+".bsid";
    std::string fishrbmidfile = OUTFILE+".bmid";
    std::string fishrbmatchfile = OUTFILE+".bmatch";

/*	glcopy[glcopycount ] = new char[strlen("-outfile")];
	strcpy(glcopy[glcopycount],"-outfile");*/

/*Add the pedfile*/
	fishrcopy[fishrcopycount ] = new char[strlen("-pedfile")+1];
	strcpy(fishrcopy[fishrcopycount],"-pedfile");
	fishrcopycount++;


	fishrcopy[fishrcopycount ] = new char[strlen(fishrpedfile.c_str())+1];
	strcpy(fishrcopy[fishrcopycount],fishrpedfile.c_str());
	fishrcopycount++;

/*^^ Added*/



/*Add the BMID file*/
	fishrcopy[fishrcopycount ] = new char[strlen("-bmid")+1];
	strcpy(fishrcopy[fishrcopycount],"-bmid");
	fishrcopycount++;


	fishrcopy[fishrcopycount ] = new char[strlen(fishrbmidfile.c_str())+1];
	strcpy(fishrcopy[fishrcopycount],fishrbmidfile.c_str());
	fishrcopycount++;

/*^^ Added*/


/*Add the BSID file*/
	fishrcopy[fishrcopycount ] = new char[strlen("-bsid")+1];
	strcpy(fishrcopy[fishrcopycount],"-bsid");
	fishrcopycount++;


	fishrcopy[fishrcopycount ] = new char[strlen(fishrbsidfile.c_str())+1];
	strcpy(fishrcopy[fishrcopycount],fishrbsidfile.c_str());
	fishrcopycount++;

/*^^ Added*/


/*Add the BMATCH file*/
	fishrcopy[fishrcopycount ] = new char[strlen("-bmatch")+1];
	strcpy(fishrcopy[fishrcopycount],"-bmatch");
	fishrcopycount++;


	fishrcopy[fishrcopycount ] = new char[strlen(fishrbmatchfile.c_str())+1];
	strcpy(fishrcopy[fishrcopycount],fishrbmatchfile.c_str());
	fishrcopycount++;

/*^^ Added*/



/*std::cout<<"FISH HERE"<<std::endl;
	for (int i = 0 ; i < fishrcopycount;i++ )
	{
		std::cout<<fishrcopy[i]<<std::endl;
	}*/


	ErrorFinderManager manager;
	manager.performConsolidation(fishrcopycount,fishrcopy);




    for(int i = 0; i < glcopycount; ++i) {
        delete [] glcopy[i];
    }
    delete [] glcopy;


    for(int i = 0; i < fishrcopycount; ++i) {
        delete [] fishrcopy[i];
    }
    delete [] fishrcopy;

    if (GERMLINE_OUTPUT == false)
		{
    		string delfile = PEDFILE.substr(0,PEDFILE.find_last_of('.'));


    		for (int i = 0 ; i <4; i ++)
    		{

    			if (i==0)
    			{
    				delfile = PEDFILE.substr(0,PEDFILE.find_last_of('.'));
    				delfile+=".bsid";
    			}
    			else if (i==1)
				{
					delfile = PEDFILE.substr(0,PEDFILE.find_last_of('.'));
					delfile+=".bmid";
				}
    			else if (i==2)
				{
					delfile = PEDFILE.substr(0,PEDFILE.find_last_of('.'));
					delfile+=".bmatch";
				}
    			else if (i==3)
				{
					delfile = PEDFILE.substr(0,PEDFILE.find_last_of('.'));
					delfile+=".log";
				}


    				if( remove( delfile.c_str())!=0)
					{
						cerr<< "Error deleting file" <<endl;
						cerr<<delfile<<endl;
					}
    		}

		}

    else//if true
		{
    		if (germline_output.length() < 3)//No path provided
    		{

    		}

    		else
    		{
    			string makedir="";
    			makedir = "mkdir -p " + germline_output;
    			system(makedir.c_str());

    			string movedir="";
    			for (int i = 0 ; i <4; i ++)
    			    {

    					if (i == 0)
    						movedir = "mv "+ PEDFILE.substr(0,PEDFILE.find_last_of('.')) + ".bmid "+ germline_output;
    					else if(i == 1)
							movedir = "mv "+ PEDFILE.substr(0,PEDFILE.find_last_of('.')) + ".bsid  "+ germline_output;
    					else if(i == 2)
							movedir = "mv "+ PEDFILE.substr(0,PEDFILE.find_last_of('.')) + ".bmatch " + germline_output;
    					else if(i == 3)
							movedir = "mv "+ PEDFILE.substr(0,PEDFILE.find_last_of('.')) + ".log " + germline_output;
    					/*else if ((i == 4) && IBD	)
    						movedir = "mv "+ PEDFILE.substr(0,PEDFILE.find_last_of('.')) + ".match " + germline_output;*/
    					else
    					{}

    					system(movedir.c_str());
    			    }
    		}
		}

   return 0;	//Piyush Debug
}

bool isPath(const std::string s)
{
	if ((s.find('.')==0)   &&  (s.find('/')==1))
		return true;
	else
		return false;
}



bool isInteger(const std::string s)
{
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+')))
	   return false ;

   char * p ;
   strtol(s.c_str(), &p, 10) ;

   return (*p == 0) ;
}
void helpShowParameters()
{
	cerr<<endl;
					cerr<<"***********************************************************************************"<<endl;
					cerr<<"*                        FISHR2 (Genetic Analysis Program)                        *"<<endl;
					cerr<<"*                                  				                                 *"<<endl;
					cerr<<"*                (C) 2014  Matthew C. Keller, Doug Bjelland, Piyush Sudip Patel   *"<<endl;
					cerr<<"*                    Institute for Behavioral Genetics (IBG)                      *"<<endl;
					cerr<<"*                       University of Colorado at Boulder                         *"<<endl;
					cerr<<"******For a complete worked through example, see \"FISHR2.complete.example.R\"*****"<<endl;
					cerr<<"***********************************************************************************"<<endl;
					cerr<<endl<<endl<<endl;
cerr<<"*For sample inputs please refer the README.md file*"<<endl<<endl;
cerr<< "List of commands:" << endl;
cerr<<"-pedfile"<<"\n\t-pedfile [pedfile ]"<<endl;
cerr<<"-mapfile"<<"\n\t-mapfile [mapfile ]"<<endl;
cerr<<"-geneticmapfile"<<"\n\t-geneticmapfile [geneticmapfile ]"<<endl;
cerr<<"-samplefile"<<"\n\t-samplefile [samplefile ]"<<endl;
cerr<<"-hapsfile"<<"\n\t-hapsfile [hapsfile ]"<<endl;

cerr<<"\n***************************************************************\n"<<endl;

cerr<<"-silent"<<"\n\tSuppress all output except for warnings and prompts."<<endl;
cerr<<"-bin_out"<<"\n\tOutput in binary format to save space."<<endl;
cerr<<"-min_cm_initial"<<"\n\tMinimum length for match to be used for imputation (in cM or MB)."<<endl;
cerr<<"-err_hom"<<"\n\tMaximum number of mismatching homozygous markers (per slice)."<<endl;
cerr<<"-err_het"<<"\n\tMaximum number of mismatching heterozygous markers (per slice)."<<endl;
cerr<<"-from_snp"<<"\n\tStart SNP (rsID)."<<endl;
cerr<<"-to_snp"<<"\n\tEnd SNP (rsID)."<<endl;
cerr<<"-map"<<"\n\tGenetic distance map."<<endl;
cerr<<"-bits"<<"\n\tSlice size."<<endl;
cerr<<"-homoz"<<"\n\tAllow self matches (homozygosity)."<<endl;
cerr<<"-homoz-only"<<"\n\tLook for autozygous/homozygous segments only, does not detect IBD."<<endl;
cerr<<"-haploid"<<"\n\tTreat input individual as two fully phased chromosomes with no recombination\n\t\toutput IDs with 0/1 suffix for chromosome destinction."<<endl;
cerr<<"-h_extend"<<"\n\tExtend from seeds if *haplotypes* match."<<endl;
cerr<<"-w_extend"<<"\n\tExtend, one marker at a time, beyond the boundaries of a found match."<<endl;
cerr<<"-reduced"<<"\n\tOutput only reduced elements."<<endl;
cerr<<"-no_suffix"<<"\n\tUse with -haploid to output the cell with no .0 or .1 suffix."<<endl;
cerr<<"-err_w"<<"\n\tUse with -wextend to allow error matches like 1,2,3 etc in extending"<<endl;
cerr<<"-hapsfile"<<"\n\tIf using hap file, use this flag. usuage: -hapsfile <hapfilename.hap>  Must also supply \".samplefile\" and \".mapfile\" file using flags -samplefile and -mapfile"<<endl;
cerr<<"-samplefile"<<"\n\tIf using hapsfile, use this flag for the corresponding sample file. usage: -samplefile <samplefilename.sample>  Must also supply \".hapsfile\" and \".mapfile\" file using flags"<<endl;
//cerr<<"-geneticmapfile"<<"\n\tOutput only reduced elements."<<endl;
//<<'\t'	<< "-geneticmapfile"	 <<'\t'	 << "if using hap file, use this flag to supply the corresponding gen file. usuage: -geneticmapfile <genfile.gen>  Must also supply \".hap\" and \".sample\" file using flags -hapsfile and -sample"<< endl
//<<"\n"<<endl
cerr<<"-window"<<"\n\t-window value[window width to calculate moving averages]"<<endl;
cerr<<"-ibd2"<<"\n\t-ibd2 <threshold value> \tCompute ibd2 and ibd4"<<endl;
cerr<<"-log-file"<<"\n\t-log-file [log file name]"<<endl;



cerr<<"-gap"<<"\n\t[max gap to consolidate two matches]"<<endl;
cerr<<"-pct-err-threshold"<<"\n\t[max percentage of errors in a match after the trim] OR -emp-pie-threshold"<<endl;
cerr<<"-ma-threshold"<<"\n\t[specifies percentile to be drawn from trulyIBD data for MA calculations] OR -empirical-ma-threshold"<<endl;
cerr<<"*"<<"\n\tNote that if both -emp-pie-threshold and empirical-ma-threshold are supplied, then -trueSNP and -trueCM will be ignored"<<endl;
cerr<<"**"<<"\n\t-output-type [ must provide any of these. it can be"<<endl;
cerr<<'\t'<<" MovingAverages  or Error1 or Error2 or Error3 or ErrorRandom1 " << endl;
cerr<<'\t'<<" or ErrorRandom2 or Error3 or ErrorRandom3 or Full "<< endl;
cerr<<'\t'<<" look at the description about how these works in wiki ]"<< endl;
cerr<<"***"<<"\n\t(optional) -holdout-ped [new ped file path] -holdout-map [new map file]"<<endl;

cerr<<"-holdout-threshold [threshold to drop a match with new ped file ]"<<endl;
cerr<<"-holdout-missing [missing value representation in new ped file] "<< endl;

cerr<<'\t'<<"-trueCM [ true match maximum cm length] " << endl;
cerr<<'\t'<<"-trueSNP [ true match SNP length]"<< endl;
cerr<<"****\t"<<"-PIE.dist.length [ can be MOL or any cm distance length "<< endl;
cerr<<'\t'<<" please refer wiki for more details on how to use this option"<< endl;
cerr<<'\t'<<"-count.gap.errors [ TRUE or FALSE to include gap errors in errors count ]"<< endl;
	exit(0);
}

/*void helpShowParameters()
{
	cerr<< "flags:" << endl
			<< " -pedfile [ped file ]" << endl
			<< " -mapfile [map file ]" << endl
			<< " -outfile [ out file ]"<<endl
			<< " -geneticmapfile [ gen file ]"<<endl
			<< " -samplefile [ sample file ]"<<endl
			<< " -hapsfile [ hap file ]"<<endl
			<< '\t' << "-silent" << '\t' << "Suppress all output except for warnings and prompts." << endl
			<< '\t' << "-bin_out" << '\t' << "Output in binary format to save space." << endl
			<< '\t' << "-min_cm_initial" << '\t' << "Minimum length for match to be used for imputation (in cM or MB)." << endl
			<< '\t' << "-err_hom" << '\t' << "Maximum number of mismatching homozygous markers (per slice)." << endl
			<< '\t' << "-err_het" << '\t' << "Maximum number of mismatching heterozygous markers (per slice)." << endl
			<< '\t' << "-from_snp" << '\t' << "Start SNP (rsID)." << endl
			<< '\t' << "-to_snp" << '\t' << "End SNP (rsID)." << endl
			<< '\t' << "-haps" << '\t' << "Print the resolved haplotypes in a seperate HAPS file." << endl
			<< '\t' << "-map" << '\t' << "Genetic distance map." << endl
			<< '\t' << "-bits" << '\t' << "Slice size." << endl
			<< '\t' << "-homoz" << '\t' << "Allow self matches (homozygosity)" << endl
			<< '\t' << "-homoz-only" << '\t' << "Look for autozygous/homozygous segments only, does not detect IBD" << endl
			<< '\t' << "-haploid" << '\t' << "Treat input individual as two fully phased chromosomes with no recombination\n\t\toutput IDs with 0/1 suffix for chromosome destinction" << endl
			<< '\t' << "-h_extend" << '\t' << "Extend from seeds if *haplotypes* match" << endl
			<< '\t' << "-w_extend" << '\t' << "Extend, one marker at a time, beyong the boundaries of a found match" << endl
			<< '\t' << "-reduced" << '\t' << "output only reduced elements" << endl
			<< '\t' << "-no_suffix" << '\t' << "use with -hapoloid to outputthe cell with no .0 or .1 suffix" << endl
			<< '\t' << "--err_w" << '\t' << "use with -wextend to allow error matches like 1,2,3 etc in extending" << endl
			<<'\t'	<< "-hapsfile"	 <<'\t'	 << "if using hap file, use this flag. usuage: -hapsfile <hapfilename.hap>  Must also supply \".sample\" and \".gen\" file using flags -sample and -geneticmapfile"<< endl
			<<'\t'	<< "-sample"	 <<'\t'	 << "if using hap file, use this flag for the corresponding sample file. usuage: -sample <samplefilename.sample>  Must also supply \".hap\" and \".gen\" file using flags -hapsfile and -geneticmapfile"<< endl
			<<'\t'	<< "-geneticmapfile"	 <<'\t'	 << "if using hap file, use this flag to supply the corresponding gen file. usuage: -geneticmapfile <genfile.gen>  Must also supply \".hap\" and \".sample\" file using flags -hapsfile and -sample"<< endl
			<<"\n"<<endl
			<<'\t'<<"-window [window width to calculate moving averages] "<< endl
			<<'\t'<<"-gap [max gap to consolidate two matches]"<< endl
			<<'\t'<<"-pct-err-threshold [max percentage of errors in a match after the trim] OR -emp-pie-threshold" << endl
			<<'\t'<<"-ma-threshold [specifies percentile to be drawn from trulyIBD data for MA calculations] OR -empirical-ma-threshold"<< endl
			<<'\t'<<" Note that if both -emp-pie-threshold and empirical-ma-threshold are supplied, then -trueSNP and -trueCM will be ignored"<< endl
			<<'\t'<<"-output.type [ must provide any of these. it can be "<< endl
			<<'\t'<<" MovingAverages  or Error1 or Error2 or Error3 or ErrorRandom1 " << endl
			<<'\t'<<" or ErrorRandom2 or Error3 or ErrorRandom3 or Full "<< endl
			<<'\t'<<" look at the description about how these works in wiki ]"<< endl
			<<'\t'<<" (optional) -holdout-ped [new ped file path] -holdout-map [new map file] "<< endl
			<<'\t'<<"-holdout-threshold [threshold to drop a match with new ped file ]"<< endl
			<<'\t'<<"-holdout-missing [missing value representation in new ped file] "<< endl
			<<'\t'<<"-log-file [log file name]"<< endl
			<<'\t'<<"-trueCM [ true match maximum cm length] " << endl
			<<'\t'<<"-trueSNP [ true match SNP length]"<< endl
			<<'\t'<<"-PIE.dist.length [ can be MOL or any cm distance length "<< endl
			<<'\t'<<" please refer wiki for more details on how to use this option"<< endl
			<<'\t'<<"-count.gap.errors [ TRUE or FALSE to include gap errors in errors count ]"<< endl
			<<'\t'<<"-ibd2 <threshold value>" <<"\t compute ibd2 and ibd4"<< endl;
	exit(0);
}*/

// end GERMLINE_0001.cpp
