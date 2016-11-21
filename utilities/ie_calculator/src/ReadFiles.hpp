/*
 * ReadFiles.hpp
 *
 *  Created on: May 24, 2016
 *      Author: root
 */

#ifndef READFILES_HPP_
#define READFILES_HPP_
#include <string>


class ReadFiles{
private:



public:
	ReadFiles();
	int readpedfile(std::string);
	int readibdfile(std::string);
	int readbmidfile(std::string);
};


#endif /* READFILES_HPP_ */
