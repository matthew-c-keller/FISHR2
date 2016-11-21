/*
 * HandleFlags.h
 *
 *  Created on: May 24, 2016
 *      Author: root
 */
#ifndef HANDLEFLAGS_H_
#define HANDLEFLAGS_H_

#include <string>

class HandleFlags
{
private:
	static std::string pedfile;
	static std::string bmidfile;
	static std::string ibdfile;
	static std::string outfile;
	static bool reduced;
	static bool silent;
	static double trim;
	static int windowsize;

public:
	HandleFlags();
	int setFlagValues(int argc, char *argv[]);	//	set the values of the private values to the run time arguments passed by the user.
	void displayFlagValues();	//	view flag values
	std::string getpedfilename();	//	view pedfilename
	std::string getbmidfilename();	//	view bmidfilename
	std::string getibdfilename();	//	view ibdfilename
	std::string getoutputfilename();
	double gettrimvalue();			// get the trim flag value
	int getwindowsize();
	bool getisreduced();
	bool getissilent();
};



#endif /* HANDLEFLAGS_H_ */
