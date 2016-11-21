/*
 * Sample.cpp
 *
 *  Created on: Oct 17, 2016
 *      Author: root
 */
#include "Sample.hpp"

Sample_3_col::Sample_3_col()
{

}
Sample_7_col::Sample_7_col()
{

}
Sample_3_col::Sample_3_col(std::string ID_1_, std::string ID_2_, std::string zero1_,std::string zero2_, std::string one_, std::string negnine_)
{
	ID_1 = ID_1_;
	ID_2 = ID_2_;
	father = zero1_;
	mother = zero2_;
	sex = one_;
	plink_phemo = negnine_;
}

Sample_7_col::Sample_7_col(std::string ID_1_, std::string ID_2_, std::string father_,std::string mother_, std::string sex_, std::string plink_phemo_)
{
	ID_1 = ID_1_;
	ID_2 = ID_2_;
	father = father_;
	mother = mother_;
	sex = sex_;
	plink_phemo = plink_phemo_;
}
