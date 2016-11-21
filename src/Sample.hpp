/*
 * Sample.hpp
 *
 *  Created on: Oct 17, 2016
 *      Author: root
 */

#ifndef SAMPLE_HPP_
#define SAMPLE_HPP_
#include <string>
#include <vector>

class Sample_3_col{

public:

	std::string ID_1;
	std::string ID_2;
	std::string father;
	std::string mother;
	std::string sex;
	std::string plink_phemo;
	Sample_3_col();
	Sample_3_col(std::string,std::string,std::string,std::string, std::string,std::string );
};

class Sample_7_col{
	//Don't include the missing column

public:

	std::string ID_1;
	std::string ID_2;
	std::string father;
	std::string mother;
	std::string sex;
	std::string plink_phemo;
	Sample_7_col();
	Sample_7_col(std::string,std::string,std::string,std::string, std::string,std::string );
};



#endif /* SAMPLE_HPP_ */
