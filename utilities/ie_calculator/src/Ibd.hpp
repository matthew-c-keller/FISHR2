/*
 * Ibd.hpp
 *
 *  Created on: Jun 1, 2016
 *      Author: root
 */

#ifndef IBD_HPP_
#define IBD_HPP_
#include <string>

class Ibd
{
public:
	std::string person1;
	std::string person2;
	unsigned long long begin;
	unsigned long long end;
	unsigned long long dist;
	Ibd();
	Ibd(std::string ,std::string ,unsigned long long ,unsigned long long,unsigned long long);


};


#endif /* IBD_HPP_ */
