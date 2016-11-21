/*
 * Ped.hpp
 *
 *  Created on: Jun 2, 2016
 *      Author: root
 */

#ifndef PED_HPP_
#define PED_HPP_
#include "string"
class Ped
{

public:
	std::string individual;	//1st column of the ped file
	std::string dnasequence;	//AGTC column of the ped file
	Ped();
	Ped(std::string,std::string);
};


#endif /* PED_HPP_ */
