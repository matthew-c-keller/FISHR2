/*
 * Bmid.hpp
 *
 *  Created on: Jun 2, 2016
 *      Author: root
 */

#ifndef BSID_HPP_
#define BSID_HPP_
#include <string>

class Bmid
{
public:

	int fileindex;	//First column of Bmid
	std::string marker;	//Second column of Bmid
	double lencm;	//Third column of Bmid
	unsigned long long location;	//Fourth column of Bmid
	unsigned long long lineno;

	Bmid();
	Bmid(int,std::string,double,unsigned long long,unsigned long long);
};


#endif /* BSID_HPP_ */
