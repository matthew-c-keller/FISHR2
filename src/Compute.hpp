/*
 * Compute.hpp
 *
 *  Created on: Oct 17, 2016
 *      Author: root
 */

#ifndef COMPUTE_HPP_
#define COMPUTE_HPP_
#include <fstream>
#include "Haps.hpp"
#include "Sample.hpp"
#include "Gen.hpp"

class Compute
{
private:

	static std::vector<size_t> lenhaps_sample[2]; // 0-> haps, 1->sample
	static std::string samplefiletype;

	std::vector<Haps> HAPS;
	std::vector<Gen> GEN;
	std::vector<Sample_3_col> SAMPLE_3_COL;
	std::vector<Sample_7_col> SAMPLE_7_COL;

public:

	Compute(std::string hapfilename,std::string samplefilename, std::string genfilename);
	void haps2vec(std::string hapfilename);
	void sample2vec(std::string samplefilename);
	void gen2vec(std::string genfilename);
	void createped(std::string hapfilename);
	void createmap(std::string samplefilename);
	void checkintegrity();




};


#endif /* COMPUTE_HPP_ */
