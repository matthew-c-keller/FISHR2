#!/bin/bash
make clean
make FISHR2
#ret = !echo $? 
#if [  $ret]
#then
#  echo "FISHR2 compilation error. See README.md"
#else
#	echo"FISHR2 compiled successfully"
#fi


make PARAMETER_FINDER
#ret = !echo $? 

#if [ $ret ]
#then
#  echo "parameter_finder compilation error. See README.md"
#else
#	echo"parameter_finder compiled successfully"
#fi
make IE_CALCULATOR
#ret = !echo $? 
#if [ $ret ]
#then
#  echo "ie_calculator compilation error. See README.md"
#else
#	echo"ie_calculator compiled successfully"
#fi
make GAP