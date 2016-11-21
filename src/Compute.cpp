/*
 * Compute.cpp
 *
 *  Created on: Oct 17, 2016
 *      Author: Piyush
 */
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <stdlib.h>     /* atof */
#include <cstdlib>
#include "Compute.hpp"
#include "Sample.hpp"
#include "Gen.hpp"

/*
 * Linear interpolation
 * http://www.ajdesigner.com/phpinterpolation/linear_interpolation_equation.php
 *
 * */
std::string Compute::samplefiletype = "NULL";
std::vector<size_t> Compute::lenhaps_sample[2];
double linearInterpolation( double x1,  double y1, double x3, double y3,  double x2);
Compute::Compute(std::string hapfilename,std::string samplefilename, std::string genfilename)
{

	haps2vec(hapfilename);
	sample2vec(samplefilename);
	gen2vec(genfilename);
	checkintegrity();	//check if the no of personrs in theHAPS file is equal to the length of the sample file
	createped(hapfilename);
	createmap(samplefilename);

}


void Compute::checkintegrity()
{

lenhaps_sample[0].push_back(HAPS[0].dnaseq.length()/4); // spaces and pairs hence 4
if (lenhaps_sample[0][0] == lenhaps_sample[1][0])
{
	//std::cout<<"SAMPLE AND HAP FILE CONSISTENT"<<std::endl;
}
else
{
	std::cout<<"INCONSISTENT SAMPLE AND HAP FILE... THE NUMBER OF PAIRS IN THE "<<
			   "HAPS FILE IS NOT EQUAL TO THE NUBER OF LINES IN THE SAMPLE FILE... EXITING "<<std::endl;
	exit(0);
}
}



void Compute::createmap(std::string samplefilename)
{
	std::string mapfilename = samplefilename.substr(0,samplefilename.find_last_of("."))+".map";
	std::ofstream fmap(mapfilename.c_str());
	if (!fmap.is_open())
		{
			std::cout<<"ERROR: OUTPUT FILE CANNOT BE OPENED FOR MAP DUMP...EXITING"<<std::endl;
			exit(0);
		}

	long long temp = -1.0;

	for (size_t i = 0; i < HAPS.size();i++)
	{
		fmap<<HAPS[i].intitialIndex<<"\t";
		fmap<<HAPS[i].marker<<"\t";

		for (size_t k = 0 ; k < GEN.size(); k++ )
		{
			temp = std::atoll(HAPS[i].location.c_str()) - std::atoll(GEN[k].position.c_str());

			if ((temp > 0) && ((k+1)==GEN.size()) )
			{
				fmap<<GEN[GEN.size() -1].Genetic_Map<<"\t";
				break;
			}

			else if ( temp < 0)
			{
				if (k > 0)//not the very beginning of the index
				{
					double li = linearInterpolation(std::atof(GEN[k-1].position.c_str()),
										std::atof(GEN[k-1].Genetic_Map.c_str()),
										std::atof(GEN[k].position.c_str()),
										std::atof(GEN[k].Genetic_Map.c_str()),
										std::atof(HAPS[i].location.c_str()));
					fmap<<li<<"\t";
					break;
				}
				else if (k==0)
				{
					fmap<< GEN[0].Genetic_Map<<"\t";
					break;
				}
			}
			else if (temp == 0)
			{
				fmap<<GEN[k].Genetic_Map<<"\t";
				break;
			}
			else if ((temp < 0) && ((k+1)==GEN.size()) )
			{
				fmap<<GEN[GEN.size() -1].Genetic_Map<<"\t";
				break;
			}
		}
		fmap <<HAPS[i].location<<"\n";
		fmap.flush();


	}
fmap.close();
}


double linearInterpolation( double x1,  double y1, double x3, double y3,  double x2)
{


	double y2 = (((x2 - x1)*(y3 - y1))/(x3-x1) ) + y1;
	return y2;

}

void Compute::createped(std::string hapfilename)
{
	std::string pedfilename = hapfilename.substr(0,hapfilename.find_last_of("."))+".ped";
	std::ofstream fped(pedfilename.c_str());
	if (!fped.is_open())
	{
		std::cout<<"ERROR: OUTPUT FILE CANNOT BE OPENED FOR PED DUMP...EXITING"<<std::endl;
		exit(0);
	}

	if (SAMPLE_7_COL.size() > 0		)
	{
		size_t begin = 0;
		size_t end = 4;	//Because " 0 0" is length 4
		size_t len_dna = HAPS[0].dnaseq.length();

		for (size_t i = 0; i < SAMPLE_7_COL.size() ; i++)
		{
			fped<<SAMPLE_7_COL[i].ID_1<<" ";
			fped<<SAMPLE_7_COL[i].ID_2<<" ";
			fped<<SAMPLE_7_COL[i].father<<" ";
			fped<<SAMPLE_7_COL[i].mother<<" ";
			fped<<SAMPLE_7_COL[i].sex<<" ";
			fped<<SAMPLE_7_COL[i].plink_phemo;
			for (size_t k = 0; k < HAPS.size();k++	)
			{
				fped<<HAPS[k].dnaseq.substr(begin,4);
			}
			fped<<"\n";
			fped.flush();
			if ((end + 4) > len_dna )
				break;
			else
				{
					begin = end;	end +=4;
				}
		}
	}

	else if (SAMPLE_3_COL.size() > 0)
	{
		size_t begin = 0;
		size_t end = 4;
		size_t len_dna = HAPS[0].dnaseq.length();

		for (size_t i = 0; i < SAMPLE_3_COL.size() ; i++)
		{
			fped<<SAMPLE_3_COL[i].ID_1<<" ";
			fped<<SAMPLE_3_COL[i].ID_2<<" ";
			fped<<SAMPLE_3_COL[i].father<<" ";
			fped<<SAMPLE_3_COL[i].mother<<" ";
			fped<<SAMPLE_3_COL[i].sex<<" ";
			fped<<SAMPLE_3_COL[i].plink_phemo;


			for (size_t k = 0; k < HAPS.size();k++	)
			{
				fped<<HAPS[k].dnaseq.substr(begin,4);
			}
			fped<<"\n";
			fped.flush();

			if ((end + 4) > len_dna )
				break;
			else
				{
					begin = end;	end +=4;
				}
		}
	}
	else
	{
		std::cout<<"Check the consistency of the .hap file...Program exiting"<<std::endl;
		exit(0);
	}
	fped.close();
}

void Compute::haps2vec(std::string hapfilename)
{
	std::ifstream hapfile (hapfilename.c_str());
	Haps haps;
	std::string line;
	std::string dnaseq;
	std::stringstream ss;

	if (hapfile.is_open())
		{
			while ( std::getline (hapfile,line) )
				{
					ss<<line;
					ss>>haps.intitialIndex;
					ss>>haps.marker;
					ss>>haps.location;
					ss>>dnaseq>>dnaseq;
					getline(ss,dnaseq);
					ss>>haps.dnaseq;
					HAPS.push_back(	Haps(haps.intitialIndex,haps.marker,haps.location,dnaseq));
					ss.str("");
					ss.clear();
				}
			hapfile.close();
		}
	else
		{
			std::cout<<hapfilename<<" : Hap file cannot be opened...Exiting "<<std::endl;
			exit(0);
		}

	//std::cout<<HAPS[0].dnaseq<<std::endl;

/*	for (unsigned int i =0; i <HAPS.size();i++)
	{
		std::cout<<HAPS[i].intitialIndex<<" "<<HAPS[i].marker<<" "<<HAPS[i].location<<HAPS[i].dnaseq<<std::endl;
	}*/
}

void Compute::sample2vec(std::string samplefilename)
{
	std::ifstream samplefile (samplefilename.c_str());//Test_7_Columns.sample
	//std::ifstream samplefile ("./src/Test_7_Columns.sample");
	std::string line;
	std::stringstream ss;
	int whitespaces = 0;
	Sample_3_col sample_3_col; // whitespaces == 2
	Sample_7_col sample_7_col; //  whitespaces == 6
	if (samplefile.is_open())
		{
			getline(samplefile,line);
		}
	else
		{
			std::cout<<samplefilename<<" : Sample file cannot be opened...Exiting "<<std::endl;
			exit(0);
		}


	for (size_t i = 0;i<line.size();i++)
		{
			if (line[i] == ' ')
				whitespaces++;
		}


	if ( whitespaces == 2	)
		{
			samplefiletype = "3_COL";
			std::getline(samplefile, line ); // get rid of 0 0 0 in the second row of the file
			size_t countlines = 0;
			if (samplefile.is_open())
				{

					while(std::getline(samplefile, line ))
						{
							countlines++;	//get the number of liens to get the integrity of the file
							ss<<line;
							ss>>sample_3_col.ID_1;
							ss>>sample_3_col.ID_2;
							SAMPLE_3_COL.push_back(Sample_3_col(sample_3_col.ID_1,sample_3_col.ID_2,"0","0","1","-9"));
							ss.str("");
							ss.clear();
						}
				}
			lenhaps_sample[1].push_back(countlines);
			samplefile.close();
		}

	//std::cout<<SAMPLE_3_COL.size()<<std::endl;
/*	for (int i =0;i<SAMPLE_3_COL.size();i++)
	{
		std::cout<<SAMPLE_3_COL[i].ID_1<<" "<<SAMPLE_3_COL[i].ID_2<<" "<<SAMPLE_3_COL[i].father<<" "<<SAMPLE_3_COL[i].mother<<" "<<SAMPLE_3_COL[i].sex<<" "<<SAMPLE_3_COL[i].plink_phemo<<std::endl;
	}*/
	if ( whitespaces == 6	)
		{
			samplefiletype = "7_COL";
			std::getline(samplefile, line ); // get rid of 0 0 0 in the second row of the file
			if (samplefile.is_open())
				{
					while(std::getline(samplefile, line ))
						{
							ss<<line;
							ss>>sample_7_col.ID_1;
							ss>>sample_7_col.ID_2;
							ss>>sample_7_col.father>>sample_3_col.father;
							ss>>sample_7_col.mother;
							ss>>sample_7_col.sex;
							ss>>sample_7_col.plink_phemo;
							SAMPLE_7_COL.push_back(Sample_7_col(sample_7_col.ID_1,sample_7_col.ID_2,sample_7_col.father, sample_7_col.mother,sample_7_col.sex, sample_7_col.plink_phemo ));
							ss.str("");
							ss.clear();
						}
				}
			samplefile.close();
		}

/*		std::cout<<SAMPLE_7_COL.size()<<std::endl;
		for (int i =0;i<SAMPLE_7_COL.size();i++)
		{
			std::cout<<SAMPLE_7_COL[i].ID_1<<" "<<SAMPLE_7_COL[i].ID_2<<" "<<SAMPLE_7_COL[i].father<<" "<<SAMPLE_7_COL[i].mother<<" "<<SAMPLE_7_COL[i].sex<<" "<<SAMPLE_7_COL[i].plink_phemo<<std::endl;
		}
*/
}



void Compute::gen2vec(std::string genfilename	)
{

	std::ifstream genfile (genfilename.c_str());
	Gen gen;
	std::string line;
	std::string dnaseq;
	std::stringstream ss;
	std::getline (genfile,line);

	if (genfile.is_open())
		{
			while ( std::getline (genfile,line) )
				{
					ss<<line;
					ss>>gen.position;
					ss>>gen.COMBINED_rate;
					ss>>gen.Genetic_Map;
					GEN.push_back(Gen(gen.position, gen.COMBINED_rate, gen.Genetic_Map));
					ss.str("");
					ss.clear();
				}
			genfile.close();
		}
	else
		{
			std::cout<<genfilename<<" : gen file file cannot be opened...Exiting "<<std::endl;
			exit(0);
		}

/*	for (int i =0;i<GEN.size();i++)
	{
		std::cout<<GEN[i].position<<" "<<GEN[i].COMBINED_rate<<" "<<GEN[i].Genetic_Map<<std::endl;
	}*/


}
