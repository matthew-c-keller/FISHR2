/*
 * Compute.hpp
 *
 *  Created on: May 24, 2016
 *      Author: root
 */

#ifndef COMPUTE_HPP_
#define COMPUTE_HPP_
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <boost/tokenizer.hpp>
#include "HandleFlags.hpp"
#include "ReadFiles.hpp"
#include "Ibd.hpp"
#include "Bmid.hpp"
#include "Ped.hpp"

class Compute{

private:
	std::ofstream foutfile;	//	File stream to output
	std::ifstream fpedfile;	//	File to read
	std::ifstream fbmidfile;
	std::ifstream fibdfile;
	unsigned long long start;	//	From ibd file
	unsigned long long stop;	//	From ibd file
	std::vector<std::string> ind1;	//	From ped file
	std::vector<std::string> ind2;	//	From ped file

	std::vector<double> ie;	//contain the elemets for ie
	std::vector<double> trim_ie;	//contain the elemets for ie

	unsigned long long len_left;	//To compute the left of column 4 of map
	unsigned long long len_right;	//To compute the right of column 4 of map
	//unsigned long long len_left_right;

	bool nomatchestrim;
	std::vector <Ibd> IBD;
	std::vector <Bmid> BMID;
	std::vector <Ped> PED;
	bool isreduced;
	bool issilent;


public:
	Compute();
	Compute(HandleFlags,ReadFiles);
	void firstpass(std::string,std::string);
	void calculate(HandleFlags,ReadFiles);
	void compute_ma_ie(unsigned long long,unsigned long long,unsigned long long,int,unsigned long long,std::string,std::string);
	void compute_ma_nm(unsigned long long,unsigned long long,unsigned long long,int,unsigned long long,std::string,std::string);
	void compute_ma_het(unsigned long long,unsigned long long,unsigned long long,int,unsigned long long,std::string );
	void compute_pie_trim2();
	void compute_trim_ie(double trim);
	void compute_trim_ie2 (double trim);
	void compute_left_right(unsigned long long,	unsigned long long );
	void compute_left_right2(unsigned long long,	unsigned long long );



	void convertIBDtovec(std::string);
	void convertBmidtovec();
	void convertPedtovec();
};



#endif /* COMPUTE_HPP_ */
