/*
 * Bmid.cpp
 *
 *  Created on: Jun 2, 2016
 *      Author: root
 */

#include "Bmid.hpp"


Bmid::Bmid()
{

}

Bmid::Bmid(int fi,std::string m,double le,unsigned long long lo,unsigned long long ln)
{
	fileindex = fi;	//First column of Bmid
	marker = m;	//Second column of Bmid
	lencm = le;	//Third column of Bmid
	location = lo;	//Fourth column of Bmid
	lineno = ln;	//Contain the line number of bmid file
}



