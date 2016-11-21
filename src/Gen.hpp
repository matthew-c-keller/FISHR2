/*
 * Gen.hpp
 *
 *  Created on: Oct 17, 2016
 *      Author: root
 */

#ifndef GEN_HPP_
#define GEN_HPP_
#include <string>

class Gen
{
public:
	std::string position;
	std::string COMBINED_rate;
	std::string Genetic_Map;
	Gen();
	Gen(std::string,std::string, std::string);
};



#endif /* GEN_HPP_ */
