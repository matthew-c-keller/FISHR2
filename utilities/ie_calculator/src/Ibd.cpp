/*
 * Ibd.cpp
 *
 *  Created on: Jun 1, 2016
 *      Author: root
 */


#include "Ibd.hpp"

Ibd::Ibd()
	{

	}

Ibd::Ibd(std::string p1,std::string p2,unsigned long long beg, unsigned long long en,unsigned long long d )
	{
		person1= p1;
		person2 = p2;
		begin = beg;
		end = en;
		dist = d;
	}
