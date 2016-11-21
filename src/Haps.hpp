/*
 * Haps.hpp
 *
 *  Created on: Oct 17, 2016
 *      Author: root
 */

#ifndef HAPS_HPP_
#define HAPS_HPP_
#include <string>
#include <vector>

class Haps{

public:

	std::string intitialIndex;
	std::string marker;
	std::string location;
	std::string dnaseq;

	Haps();
	Haps(std::string,std::string,std::string,std::string );

};


#endif /* HAPS_HPP_ */
