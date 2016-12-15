//#include "Consolidator.hpp"
//#include "ErrorCalculator.cpp"
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <exception>
#include <algorithm> 
#include <assert.h>
#include <sstream>
#include <iomanip>
#include <cfloat>
#include "Consolidator.hpp"
#include "ErrorCalculator.hpp"
#include "BasicDefinitions.h"

//using namespace std;

template <typename T>
std::string NumberToString ( T Number )
{
	std::ostringstream ss;
	ss << Number;
	return ss.str();
}

bool Consolidator::compareFunction(SNP_lrf s1, SNP_lrf s2)
{

	return (s1.start<s2.start);

}

void Consolidator::sortMatches()
{
         for(int i=0;i<person_count;i++)
        {
                for(int j=i;j<person_count;j++)
                {


                        std::sort(m_matches[i][j].begin(),m_matches[i][j].end(),compareFunction);

                }
        }
}

void Consolidator::readMatches(std::string path,int pers_count, ErrorCalculator& eCalculator, int trueSNP, float trueCM, int snipExtend,std::string pedFile )//path->BMATCHFILE file
{


     person_count=pers_count;
     if(pers_count<=0)
     {
       std::cerr<<"wrong BSID file, check it, reading ped file failed"<<std::endl;
       return;
     }
     try
     {
        person_count=pers_count;
        m_matches.resize(pers_count+1);
        m_trueMatches.resize( pers_count + 1 );
        for(int i=0;i<pers_count;++i)
        {
             m_matches[i].resize(pers_count+1);
             m_trueMatches[i].resize( pers_count + 1 );
        }
        unsigned int pid[2];
        unsigned int sid[2];
        unsigned int dif,hom[2];
        std::ifstream file_bmatch(path.c_str(),std::ios::binary);
        if( !file_bmatch )
         {
        	std::cerr<<"unable to open the bmatch file, exiting the program" << std::endl;
        	std::cout<<"we found the bmatch file"<<std::endl;
             exit( -1 );

         }


unsigned long long counttt=0;
        while ( !file_bmatch.eof())
        {
        	counttt++;

        		pid[0] = -1;
                file_bmatch.read( (char*) &pid[0] , sizeof( unsigned int ) );
                if ( pid[0] == -1 ) continue;
                file_bmatch.read( (char*) &pid[1] , sizeof( unsigned int ) );
                file_bmatch.read( (char*) &sid[0] , sizeof( unsigned int ) );
                file_bmatch.read( (char*) &sid[1] , sizeof( unsigned int ) );
                file_bmatch.read( (char*) &dif , sizeof( int ) );
                file_bmatch.read( (char*) &hom[0] , sizeof( bool ) );
                file_bmatch.read( (char*) &hom[1] , sizeof( bool ) );
                 if(pid[0]>=pers_count||pid[1]>=pers_count)
                 {
                	 std::cerr<<"problem with bsid file, check it please"<<std::endl;
                      return;

                 }
                 SNP_lrf snp;


                 snp.start=sid[0];
                 snp.end=sid[1];

//
/*

                 std::cout<<" pid[0]= "<<pid[0];
                 std::cout<<"\tpid[1]= "<<pid[1];
                 std::cout<<"\tsid[0]= "<<sid[0];
                 std::cout<<"\tsid[1]= "<<sid[1];
                 std::cout<<"\tdif= "<<dif;
                 std::cout<<"\thom[0]="<<hom[0];
                 std::cout<<"\thom[1]= "<<hom[1];
                 std::cout<<"\tsnp.start = "<<sid[0];
                 std::cout<<"\tsnp.end=" <<sid[1];
                 std::cout<<std::endl;
//
*/


                 if(pid[0]<=pid[1])
                       m_matches[(pid[0])][(pid[1])].push_back(snp);
                 else  
                       m_matches[(pid[1])][(pid[0])].push_back(snp);
                 if( ( eCalculator.getCMDistance( sid[ 1 ] ) - 
                                eCalculator.getCMDistance( sid[ 0 ] ) ) >= trueCM && 
                                 ( sid[ 1 ] - sid[ 0 ] ) >= trueSNP &&  pid[0] != pid[1] )
                 {
                     if(pid[0]<=pid[1])
                       m_trueMatches[(pid[0])][(pid[1])].push_back(snp);
                     else
                       m_trueMatches[(pid[1])][(pid[0])].push_back(snp);
   
                 }	

        }


        file_bmatch.close();
    }
    catch(std::exception &e)
    {
    	std::cerr<<"Error:"<<e.what()<<std::endl;
       exit( -1 );
    }
        
}
//for reading in a file of user supplied snp weights
void Consolidator::readUserSuppliedSnpWeights( std::string path ){
  float temp;
  std::stringstream ss;
  std::string item,line;
  try
  {
	  std::ifstream file_weights(path.c_str());
    if( !file_weights )
    {
    	std::cerr<< "Weight file cannot be opened, check the path or the file type. Exiting program." <<std::endl;
      exit( -1 );
    }
    while ( getline(file_weights , line) )
    {
      ss.clear(); 
      ss.str( line );
      ss >> temp;
      user_supplied_snp_weights.push_back( temp );
    }
    file_weights.close();
  }
  catch(std::exception &e)
  {
	  std::cerr<<e.what()<<std::endl;
    exit( -1 );
  }
}

void Consolidator::performConsolidation(ErrorCalculator& eCalculator, int gap,int min_snp,float min_cm,int extendSnp)
{
	//cout<<"in consolidator extendsnp"<<extendSnp<<endl;
         int consolidations = 0, removed = 0;
for(int i=0;i<person_count;++i)//for each person
        {
                for(int j=i;j<person_count;++j)//compare with each other person
                {
                        int temp1=-1,temp2=-1;
                         for(int l=0;l<m_matches[i][j].size();++l)//for each match
                         {
                               temp1= m_matches[i][j][l].start;//cout<<"temp1 before = "<<temp1<<endl;
                               temp2= m_matches[i][j][l].end;//cout<<"temp2 before = "<<temp2<<endl;

                                if(temp2==-1||temp1==-1){continue;}

                                for(int k=l+1;k<m_matches[i][j].size();++k) //for each other match
                                {
                                    if((m_matches[i][j][k].start-temp2-1)<=gap)
                                    {
                                        ++consolidations;
                                        temp2=m_matches[i][j][k].end;
                                        m_matches[i][j][l].end=temp2;
                                        m_matches[i][j][k].end=-1;
                                    }
                                    else break;

                               }
     //this may be what is causing our initial drops to never show up...
     
                               if( ( (temp2-temp1)<min_snp) || ( (eCalculator.getCMDistance(temp2)-eCalculator.getCMDistance(temp1))<min_cm) )
                              {
      
                                       ++removed;
                          //             m_matches[i][j][l].end=m_matches[i][j][l].start=-1;
                              }
}
}
}
        /*std::string str = " \n Number of Consolidations: " +
                            NumberToString( consolidations );
        str =  str + " \n Number of matches removed due to initial length: "
                     +  NumberToString( removed );*/
        /*new*/
global_initial = removed;
        consolidated_str = "Number of Consolidations: " + NumberToString( consolidations );
        initial_drop_str = "Number of matches removed due to initial length: " +  NumberToString( removed );
        /*wen*/
       // eCalculator.log( str );
       
}
//overloaded version
void Consolidator::performTrim(ErrorCalculator& e_obj,int window,
                               int ma_snp_ends, float ma_threshold,
                               int min_snp,float min_cm,
                               float per_err_threshold, std::string option,
                               float hThreshold, bool holdOut,float empirical_threshold, float empirical_pie_threshold,int extendSnp)//<piyush> added the param int EXTENDSNP for calculating moving window avg)
{
  int removed1 =0, removed2 = 0, removed3 = 0, removed4 = 0;
  int not_removed = 0;
  int total_count = global_initial;
  bool wrongOption = false;
  float per_err_threshold1;
  if(empirical_pie_threshold >= 0.0){
    per_err_threshold1 = empirical_pie_threshold;
  } else {
    per_err_threshold1 = getPctErrThreshold( per_err_threshold );
  }
  std::stringstream sstr;
  sstr << std::fixed << std::setprecision(10) << per_err_threshold1;
  std::string per_err_value = sstr.str();
  emp_pie_thresh_str = "empirical pie threshold is : " + per_err_value  + " \n";
  float hThreshold1 = 0;
  if( holdOut )
  {
    hThreshold1 = getHoldOutThreshold( hThreshold );
  }
  
  per_err_threshold = per_err_threshold1;
  hThreshold = hThreshold1;

  for(int i=0;i<person_count;i++)
  {
    for(int j=i;j<person_count;j++)
    {
      for(int l=0;l<m_matches[i][j].size();l++)
      {
        total_count++;
        if(m_matches[i][j][l].end==-1)
        {
          continue; 
        }

        int temp1=m_matches[i][j][l].start;
        //cout<<"temp1 start begin= "<<temp1<<endl;

        int temp2=m_matches[i][j][l].end;
        //cout<<"temp2 end begin= "<<temp2<<endl;


if (extendSnp != 0)
{

        /*<piyush1>*/
        if(temp1-extendSnp <0)
                                     {
                                     	temp1=0;
                                     	//cout<<"New value of temp1= "<<temp1<<endl;
                                     }
                                     else
                                     {
                                     	temp1=m_matches[i][j][l].start-extendSnp;
                                     	//cout<<"New value of temp1= "<<temp1<<endl;
                                     }

                                     if(temp2+extendSnp > 4443)// change this constant
                                     {
                                     	temp2=m_matches[i][j][l].end;
                                     	//cout<<"New value of temp2= "<<temp2<<endl;
                                     }
                                     else
                                     {
                                     	temp2=m_matches[i][j][l].end+extendSnp;
                                     	//cout<<"New value of temp2= "<<temp2<<endl;
                                     }


        /*till here*/
                                     /*cout<<"temp1 start after= "<<temp1<<endl;
                                     cout<<"temp2 end after= "<<temp2<<endl;*/


                                     //cout<<"perform trim temp1= "<<temp1<<endl;
                                     //cout<<"perform trim temp2= "<<temp2<<endl;


}

        int pers1 = i, pers2 = j;
        if( option.compare( "ErrorRandom1" ) == 0 || option.compare( "ErrorRandom2" ) == 0 || option.compare( "ErrorRandom3" ) == 0 )
        {
          pers1 = std::rand() % e_obj.getNoOfPersons();
          pers2 = std::rand() % e_obj.getNoOfPersons();
          if( pers1 > pers2 )
          {
            pers1 = pers1 + pers2;
            pers2 = pers1 - pers2;
            pers1 = pers1 - pers2;
          }
        }

        std::vector<std::vector<int> > errors=e_obj.checkErrors(pers1, pers2, temp1, temp2);
        std::vector<int>finalErrors=e_obj.getFinalErrors(errors);//<piyush for errors>

        //cout<<"finalErrors size= "<<finalErrors.size()<<endl;
        /*Inject implied error at start/end of SH here*/
        std::vector<int>::iterator it;
        it = finalErrors.begin(); //go to the start of the vector
        if(finalErrors[0] != 1){
          finalErrors.insert(it,1); //inject an error at position 1, if not already there
        }	
        /*End inject implied error section*/

        std::vector<int>trimPositions;
        std::vector<float>movingAverages;
        float threshold;
        if( (e_obj.isInitialCmDrop(temp1,temp2,min_cm)) || ((temp2-temp1) < min_snp) ){ //initial drop. Don't calculate MA
          trimPositions.push_back(temp1);
          trimPositions.push_back(temp2);
          trimPositions.push_back(1);
        }else{

          movingAverages = e_obj.getMovingAverages(finalErrors,temp1,temp2,window,extendSnp);//<piyush> get moving averages are calculated from this part
          if(empirical_threshold < 0.0){
            threshold = e_obj.getCutoff();
          } else {
            threshold = empirical_threshold;
          }
          trimPositions = e_obj.getTrimPositions(movingAverages,temp1,temp2,threshold,min_cm); 

        }
        //-----------------

        int beforeTrimStart = temp1;	
        int beforeTrimEnd = temp2;
        m_matches[i][j][l].end = temp2 = temp1+trimPositions[1];
        m_matches[i][j][l].start = temp1 = temp1+trimPositions[0];
        int del0 = trimPositions[0];
        int del1 = trimPositions[1];
        float per_err = e_obj.getThreshold(finalErrors,del0,del1,ma_snp_ends);

        //add new weighted option
        /*
         For this new option, we only output SH that are not dropped. So, the output is finalOutput + weighted column.
        */
        if( (option.compare("weightedOutput") == 0) || (option.compare("weightedOutputBP") == 0) ){
          int snp1 = 0, snp2 = 0, hlength = 0;
          float noOfOppHom = 0;
          if( holdOut )
          {
            snp1 = e_obj.getNewSnp( temp1 );
            snp2 = e_obj.getNewSnp( temp2 );
            hlength = snp2 - snp1;
            if( hlength <= 0 )
            {
              hlength = 1;
            }
            noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, m_matches[i][j][l].start, m_matches[i][j][l].end );
          }
          if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){ 
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            removed4++;
            continue;
          }
          if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ //removed2 a tpos.size of 3 indicates trimming due ot cM 
            removed2++;
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            continue;
          }
          if( per_err > per_err_threshold){
            removed1++;
            continue;
          }
          if( holdOut && hThreshold < ( noOfOppHom ) / hlength ){
            removed3++;
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            continue; 
          } //removed3

          not_removed++;
          m_weighted_sh.push_back(Weighted_SH(temp1,temp2,i,j)); //build the vector of SH that passed
          continue;
        }//end weghtedOutput
        /*Add new finalErrorsOutput*/
        if( (option.compare("finalErrorsOutput") == 0) ){
          int snp1 = 0, snp2 = 0, hlength = 0;
          float noOfOppHom = 0;

          if( holdOut )
          {
            snp1 = e_obj.getNewSnp( temp1 );
            snp2 = e_obj.getNewSnp( temp2 );
            hlength = snp2 - snp1;
            if( hlength <= 0 )
            {
              hlength = 1;
            }
            noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, m_matches[i][j][l].start, m_matches[i][j][l].end );
          }


          if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){ 
            std::vector<float>movingAverages;
            temp1 = beforeTrimStart;
            temp2 = beforeTrimEnd;
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            removed4++;
            continue;
          }

          if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ //removed2 a tpos.size of 3 indicates trimming due ot cM
            removed2++;
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            continue;
          }
          if( per_err > per_err_threshold){
            removed1++;
            continue;
          }

          if( holdOut && hThreshold < ( noOfOppHom ) / hlength ){
            removed3++;
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            continue; 
          } //removed3
          not_removed++;
          e_obj.finalErrorsOutput(i,j,temp1,temp2,min_cm,per_err);// <piyush>this is where the final output is written is called at
          continue;
        }//end finalErrorsOutput
        if( (option.compare("FullPlusDropped") == 0) ){
          int snp1 = 0, snp2 = 0, hlength = 0;
          float noOfOppHom = 0;

          if( holdOut )
          {
            snp1 = e_obj.getNewSnp( temp1 );
            snp2 = e_obj.getNewSnp( temp2 );
            hlength = snp2 - snp1;
            if( hlength <= 0 )
            {
              hlength = 1;
            }
            noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, m_matches[i][j][l].start, m_matches[i][j][l].end );
          }


          if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){	
            std::vector<float>movingAverages;
            temp1 = beforeTrimStart;
            temp2 = beforeTrimEnd;
            e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,1);//standardize the error codes
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            removed4++;
            continue;
          }

          if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ //removed2 a tpos.size of 3 indicates trimming due ot cM
            e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,2); 
            removed2++;
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            continue;
          }
          if( per_err > per_err_threshold){
            e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,3);
            removed1++;
            continue;
          }

          if( holdOut && hThreshold < ( noOfOppHom ) / hlength ){
            e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,4);
            removed3++;
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            continue; 
          } //removed3
          not_removed++;
          e_obj.finalOutPut(i,j,temp1,temp2,min_cm);
          continue;
        } //end FullPlusDropped

        //Calculate Error1
        if( (option.compare("Error1") == 0 ) || (option.compare("ErrorRandom1") == 0) || (option.compare("Error") == 0) ){

          if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){ //dropped before trimming
            //don't bother printing out ma for this one. But go back and change it so that it doesn't actually calc it
            std::vector<float>movingAverages;//null	   
            //trying something special in this case. This can be removed once idrops aren't being trimmed
            //test code
            temp1 = beforeTrimStart;
            temp2 = beforeTrimEnd;
            //
            e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,1);
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            removed4++; //seems ok
            continue;
          } 
          if( (( temp2-temp1 ) < min_snp) || ((trimPositions.size() == 3) && (trimPositions[2] == 2) ) ) //dropped after trimming
          {
            e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,2);
            ++removed2;
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            continue;
          }
          if( per_err > per_err_threshold ) //dropped due to pie
          {
            e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,3);
            ++removed1;
            m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
            continue;
          }
          not_removed++;
          e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,0);//no drop
          continue;
        }//end error1

        int snp1 = 0, snp2 = 0, hlength = 0;
        float noOfOppHom = 0;
        if( holdOut )
        {
          snp1 = e_obj.getNewSnp( temp1 );
          snp2 = e_obj.getNewSnp( temp2 );
          hlength = snp2 - snp1;
          if( hlength <= 0 )
          {
            hlength = 1;
          }
          noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, m_matches[i][j][l].start, m_matches[i][j][l].end );
        }
        //update drop order 2/26/14
        if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) )
        {
          ++removed2;
          m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
          continue;
        }

        if( per_err > per_err_threshold )
        {

          ++removed1;
          m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
          continue;
        }
        //probably not removed?
        not_removed++;
        if( option.compare("MovingAverages")==0 ) //make this ma2
        {
          if( holdOut)
          {
            e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp,min_cm,movingAverages,trimPositions,per_err, noOfOppHom, hlength );
          }  
          else
          {
            e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm,movingAverages, trimPositions,per_err );
          }
          continue;
        }

        if(option.compare("Error2")==0 || option.compare( "ErrorRandom2" ) == 0)
        {
          if( holdOut)
          { 
            e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
          }
          else
          {
            e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err);
          }
          continue;
        }
        if ( holdOut && hThreshold < ( noOfOppHom ) / hlength )
        {
          ++removed3;
          m_matches[i][j][l].start= m_matches[i][j][l].end=-1;
          continue;
        }
        if( option.compare("Error3")==0 || option.compare( "ErrorRandom3" ) == 0  )
        {
          e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
        }
      }//l
    }//j
  }//i


  /*ENTERING TESTING AREA DEC 4th 2014*/
  /*****************************
  ******************************/


  /*Now, let's handle weighted output if need be*/
  if( option.compare("weightedOutput") == 0 ){
    float snp_average_count = 0.0;
    int start_position;
    int end_position;
    int genome_length;
    if(isUserSuppliedWeights()){ //the user has supplied their own weights.
      //in this case, the min and max values correspond to the number of lines in the input file,
      //since each line represents a snp. So the min is always 0, and the max is always the number of lines-1.
      start_position = 0;
      end_position = user_supplied_snp_weights.size() - 1;
    }else {
      start_position = find_genome_min();
      end_position = find_genome_max();
    }//end else
    genome_length = (end_position - start_position)+1;
    genome_vector.resize(genome_length,0);
    if(isUserSuppliedWeights()){
      for(int i = 0; i < user_supplied_snp_weights.size(); i++){
        update_genome(i,user_supplied_snp_weights[i]);
      }
    }else{
      /*This next for loop adds one to each snp in a SH. Bypass it if the user gives a files of weights*/
      for(int i = 0; i < m_weighted_sh.size(); i++){
        update_genome(m_weighted_sh[i].snp1, m_weighted_sh[i].snp2);
      }
    }
    //this part is next...will probably need to add stuff to that weighted object...
    snp_average_count = average_snp_count();
    for(int i = 0; i < m_weighted_sh.size(); i++){
      m_weighted_sh[i].snp_weight = update_snp_weight(m_weighted_sh[i].snp1, m_weighted_sh[i].snp2);
    }
    for(int i = 0; i < m_weighted_sh.size(); i++){
      m_weighted_sh[i].final_weight = ( snp_average_count / (m_weighted_sh[i].snp_weight));
      e_obj.weightedOutput(m_weighted_sh[i].per1, m_weighted_sh[i].per2, m_weighted_sh[i].snp1, m_weighted_sh[i].snp2, m_weighted_sh[i].final_weight);
    }
  }

  if (option.compare("weightedOutputBP") == 0){

  //begin new test code section here: Dec 4th 2014
  int genome_length = e_obj.getGenomeBPLength();    
  float adjusted_genome_length = genome_length / 1000.0; //L using kbp for now
  int genome_min = e_obj.getMinimumBP(); std::cout<<"genome_min= "<<genome_min<<std::endl;
  int genome_max = e_obj.getMaximumBP(); std::cout<<"genome_max= "<<genome_max<<std::endl;
  int genome_size_snps = (find_genome_max() - find_genome_min())+1; //used for genome_vector
  float wprime_numerator = 0.0;  //This is Ci / L
  float total_sh_length_sum = 0.0;
  float w2prime_denominator = 0.0;

  genome_vector.resize(genome_size_snps,0); //resize and zero out the genome. shit that needs to be snps.

  //update all of the snp counts in the genome. This looks fine.
  for(int i = 0; i < m_weighted_sh.size(); i++){
    update_genome(m_weighted_sh[i].snp1, m_weighted_sh[i].snp2); 
  }

  //calculate the w' numerator by summing up all of the snp counts and dividing by the genome length.
  //WARNING: This can cause wprime_numerator to overflow. Currently using kbp units to avoid this, but
  //this needs to be addressed.
  for(int i = 0; i < genome_vector.size(); i++){
    wprime_numerator += genome_vector[i] / adjusted_genome_length;
  }
  
  //Calculate w' for each SH.
  for(int i = 0; i < m_weighted_sh.size(); i++){
    float wprime_denominator = 0.0;
    m_weighted_sh[i].mbp_length = (e_obj.getSHBPLength(m_weighted_sh[i].snp1, m_weighted_sh[i].snp2)/1000.0);
    wprime_denominator = get_snps_over_range(m_weighted_sh[i].snp1, m_weighted_sh[i].snp2, m_weighted_sh[i].mbp_length);
    m_weighted_sh[i].wprime = wprime_numerator / wprime_denominator;
  }

  //This is the total length of all SH. This can probably overflow as well...ugh.
  for(int i = 0; i < m_weighted_sh.size(); i++){
    total_sh_length_sum += m_weighted_sh[i].mbp_length;
  }
  
  //Calculate the w2prime denominator - this value is a constant
  for(int i = 0; i < m_weighted_sh.size(); i++){
    float temp = m_weighted_sh[i].mbp_length * m_weighted_sh[i].wprime;
    w2prime_denominator += temp / total_sh_length_sum;
  }

  //Calculate and output w2' for each SH
  for(int i = 0; i < m_weighted_sh.size(); i++){
    m_weighted_sh[i].w2prime = (m_weighted_sh[i].wprime) / w2prime_denominator;
    e_obj.weightedOutput(m_weighted_sh[i].per1, m_weighted_sh[i].per2, m_weighted_sh[i].snp1, m_weighted_sh[i].snp2, m_weighted_sh[i].w2prime);
  }
}

  /*End weighted output*/



  /*END TESTING AREA DEC 4th 2014*/
  /*****************************
  ******************************/
  ma_drop_str = "No of matches removed due to length of trimming by moving averages: " + NumberToString( removed2 );
  pie_drop_str = "No of matches removed due to percentage error: " + NumberToString( removed1 );
  if(holdOut){
  //  str = str+ " \n No of matches removed due hold out ped file checking: "+ NumberToString( removed3 );
  }
  //begin log output
  std::string parameter_string_1 = "\n\n**********Parameters used in program**********\n";
  e_obj.log(parameter_string_1);
  e_obj.log(emp_ma_thresh_str); //keep
  e_obj.log(emp_pie_thresh_str);//keep
  parameter_string_1 = "**********************************************\n\n";
  e_obj.log(parameter_string_1);
  std::string total_count_str = "The total number of SH in the input file was: " + NumberToString(total_count);
  e_obj.log(total_count_str);
  e_obj.log(consolidated_str);
  e_obj.log(initial_drop_str);
  //  e_obj.log(ibg_str);
  e_obj.log(ma_drop_str);
  e_obj.log(pie_drop_str);

  //std::cerr<<IBD2_4_COUNT<<"\t"<<IBD<<std::endl;

  if (!IBD)
  {
  final_sh_str = "Total number of SH that were not dropped is: " + NumberToString(not_removed);
  e_obj.log(final_sh_str);
  }
  else
  {
	  //final_sh_str = "Total number of SH that were not dropped is: " + NumberToString(not_removed+IBD2_4_COUNT);
  }
  TOTAL_IBD1_COUNT = not_removed;
  //e_obj.log(final_sh_str);
}//end performTrim

void Consolidator::performHoldOutTrim( ErrorCalculator& ecal, 
     float threshold, std::string hMiss,
     std::string option )
{
   threshold = getHoldOutThreshold( threshold );
   for(int i=0;i<person_count;i++)
   {
      for(int j=i;j<person_count;j++)
      {
          for(int l=0;l<m_matches[i][j].size();l++)
          {
              if(m_matches[i][j][l].end==-1)
              {
                  continue;
              }
              int temp1=m_matches[i][j][l].start;
              int temp2=m_matches[i][j][l].end;
              int  snp1 = ecal.getNewSnp( temp1 );
              int snp2 = ecal.getNewSnp( temp2 );
              int length = snp2 - snp1;
              float noOfOppHom = ecal.getOppHomThreshold( i, j, m_matches[i][j][l].start, m_matches[i][j][l].end );
              if( length != 0 && threshold < (noOfOppHom) / length )
              {
                 m_matches[i][j][l].start = m_matches[i][j][l].end = -1;
              }
         }
     }
  }
}
void Consolidator::finalOutPut(ErrorCalculator &e,float min_cm, int min_snp )const
{



     for(int i=0;i<person_count;i++)
        {
                for(int j=i;j<person_count;j++)
                {

                         for(int l=0;l<m_matches[i][j].size();l++)
                         { 
                        	 if(m_matches[i][j][l].start==-1||m_matches[i][j][l].end==-1 || ( m_matches[i][j][l].end - m_matches[i][j][l].start ) < min_snp ) continue;
                                e.finalOutPut(i,j,m_matches[i][j][l].start,m_matches[i][j][l].end ,min_cm);                
                         }

}

}

}
void Consolidator::findTrueSimplePctErrors( ErrorCalculator &e_obj, float PIElength, bool holdOut,int window, float ma_threshold, float empirical_ma_threshold )
{
   for(int i=0;i<person_count;i++)
   {
      for(int j=i;j<person_count;j++)
      {
          for(int l=0;l<m_trueMatches[i][j].size();l++)
          {
              if(m_trueMatches[i][j][l].end==-1)
              {
                  continue;
              }
    //-------------------------------------------------------------------------------------------------
	int t1 = m_trueMatches[ i ][ j ][ l ].start +
						  ( m_trueMatches[ i ][ j ][ l ].end -
							m_trueMatches[ i ][ j ][ l ].start ) * 0.25;
	int t2 = m_trueMatches[ i ][ j ][ l ].end -
			   ( m_trueMatches[ i ][ j ][ l ].end -
				m_trueMatches[ i ][ j ][ l ].start ) * 0.25;

	std::vector<std::vector<int> > trueErrors=e_obj.checkErrors( i, j, t1, t2);
	std::vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );
    //This section handles finding the maximum moving averages amongst trulyIBD segments
    std::vector<float> av;
    float current_max;
    if(empirical_ma_threshold < 0.0){
             av = e_obj.getTrueMovingAverages(finalTrueErrors,t1,t2,window);//<true moving avg>
             current_max = av[0];
             for(int q = 1; q < av.size(); q++){
                   if(av[q] > current_max){
                           current_max = av[q];
                   }
              }
             e_obj.addMaxAverage(current_max);
    }
    //------------------------------------------------------------------------------
             int temp1 = m_trueMatches[i][j][l].start;
             int temp2 = m_trueMatches[i][j][l].end;
             float startCM = e_obj.getCMDistance( temp1 );
             float endCM = e_obj.getCMDistance( temp2 );
             float mid1CM = startCM + ( endCM - startCM ) / 2 - PIElength / 2;
             float mid2CM = startCM + ( endCM - startCM ) / 2 + PIElength / 2;
             while( e_obj.getCMDistance( temp1 ) <= mid1CM || e_obj.getCMDistance( temp2 ) >=mid2CM )
             {
                if( e_obj.getCMDistance( temp1 ) <= mid1CM )
                {
                  ++temp1;
                }
                if( e_obj.getCMDistance( temp2 ) >=mid2CM )
                {
                  --temp2;
                }
             }
             std::vector<std::vector<int> > errors=e_obj.checkErrors( i, j, temp1, temp2);

             std::vector<int>finalErrors=e_obj.getFinalErrors( errors );
//                  float per_err = e_obj.getThreshold(finalErrors,temp1, temp2, 0 );
             float per_err = e_obj.getThreshold(finalErrors,temp1,temp2); //overload!
                  m_errors.push_back( per_err );
                 if( holdOut  )
                 {
                        float oppHom =
                                     ( e_obj.getOppHomThreshold( i, j,
                                               temp1,
                                             temp2 ) ) / ( temp2 -temp1 );
                        m_holdOutErrors.push_back( oppHom );
                 }

          }
     }
  }
   
   //this section actually handles the sorting of the max averages, and the setting of the user supplied percentile.
   std::vector<float>maxes;
   float cutoff = empirical_ma_threshold; //assume the user wanted to supply a value. This value will be overwritten shortly if they did not.
   if(empirical_ma_threshold < 0.0){
   maxes = e_obj.getMaxAverages();
   std::sort(maxes.begin(),maxes.end());
   e_obj.setMaxAverage(maxes);
   cutoff = e_obj.getXthPercentile(ma_threshold); //<-make that an actual user input value
   }
   
   e_obj.setCutoff(cutoff);//set the actual threshold to be used when calculating MA in all other SH

   if(empirical_ma_threshold < 0.0){
   //std::string outt = "\n User supplied ma-threshold is: " + NumberToString(ma_threshold);
   //outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
   //e_obj.log(outt);
    ma_thresh_str = "User supplied ma-threshold is: " + NumberToString(ma_threshold);
    emp_ma_thresh_str = "Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
   } else {
   //std::string outt = "\n User supplied empirical-ma-threshold is: " + NumberToString(cutoff);
   //outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
    emp_ma_thresh_str = "Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
   //e_obj.log(outt);
   }
   //----------------------------------------


   std::sort( m_errors.begin(), m_errors.end() );
   std::sort( m_holdOutErrors.begin(), m_holdOutErrors.end() );
   /*std::string str =  " \n No of segments deemed to be IBD for finding empirical error threshold "
                      + NumberToString( m_errors.size() );*/
   ibg_str = "No of segments deemed to be IBD for finding empirical error threshold "
                      + NumberToString( m_errors.size() );
   /*str  = str + " \n No of segments deemed to be IBD for finding empirical error threshold in hold out are "
                      + NumberToString( m_holdOutErrors.size() );*/
        //e_obj.log( str );

}

void Consolidator::findTruePctErrors( ErrorCalculator &e_obj,int ma_snp_ends, bool holdOut,int window,float ma_threshold, float empirical_ma_threshold )
{
   for(int i=0;i<person_count;i++)
   {
      for(int j=i;j<person_count;j++)
      {
          for(int l=0;l<m_trueMatches[i][j].size();l++)
          {
              if(m_trueMatches[i][j][l].end==-1)
              {
                  continue;
              }
     //handle moving averages calculation
     //
     int t1 = m_trueMatches[ i ][ j ][ l ].start +
                          ( m_trueMatches[ i ][ j ][ l ].end -
                            m_trueMatches[ i ][ j ][ l ].start ) * 0.25;
              int t2 = m_trueMatches[ i ][ j ][ l ].end -
                           ( m_trueMatches[ i ][ j ][ l ].end -
                            m_trueMatches[ i ][ j ][ l ].start ) * 0.25;	
     //now we have the positions of the first and last 25% of the truly ibd SH
     //all that's left to do is to pass them into the moving averages function, and obtain the max ma
     //then store that in a vector, sort them, and find the xth percentile of that vector. That will be
     //the ma that we use later
     //for that "finalErrors" parameters, need to get the number of errors along the truly IBD SH first...
              std::vector<std::vector<int> > trueErrors=e_obj.checkErrors( i, j, t1, t2);
              std::vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );

     //handles MA calculations
     std::vector<float> av;
     float current_max;
     if(empirical_ma_threshold < 0.0){
     av = e_obj.getTrueMovingAverages(finalTrueErrors,t1,t2,window);
              current_max = av[0];
              for(int q = 1; q < av.size(); q++){
                   if(av[q] > current_max){
                           current_max = av[q];
                   }
              }
              e_obj.addMaxAverage(current_max);
     }
 	     
     //
              int temp1 = m_trueMatches[ i ][ j ][ l ].start +
                          ( m_trueMatches[ i ][ j ][ l ].end -
                            m_trueMatches[ i ][ j ][ l ].start ) * 0.15; //Should probably stop doing this
              int temp2 = m_trueMatches[ i ][ j ][ l ].end - 
                           ( m_trueMatches[ i ][ j ][ l ].end - 
                            m_trueMatches[ i ][ j ][ l ].start ) * 0.15;
              int start =0, end =0, fend = ( temp2 -temp1 )  ;
               
 //since we are using MOL at this point, this will pick out a random SH from the set of non-truly IBD SH 
 //and use that length to define the region over which we find PIE. Unless you are changing something with MOL,
 //don't ever read this next block
                  int randPers1, randPers2, pos;
                  randPers1 = std::rand() % person_count;
                  randPers2 = std::rand() % person_count;
                  if( randPers1 > randPers2 )
                  {
                     randPers1 = randPers1 + randPers2;
                     randPers2 = randPers1 - randPers2;
                     randPers1 = randPers1 - randPers2;
                  }
                  while( m_matches[ randPers1 ][ randPers2 ].size() <= 0 )
                  {
                    randPers1 = std::rand() % person_count;
                    randPers2 = std::rand() % person_count;
                   if( randPers1 > randPers2 )
                   {
                      randPers1 = randPers1 + randPers2;
                      randPers2 = randPers1 - randPers2;
                      randPers1 = randPers1 - randPers2;
                   }

                  }
                  pos = std::rand() % m_matches[ randPers1 ][ randPers2 ].size();
                  int len = m_matches[ randPers1 ][ randPers2 ][ pos ].end 
                            - m_matches[ randPers1 ][ randPers2 ][ pos ].start;
                  if( len >= fend || len <= 0)
                  {
                      continue;
                  } 
                  //temp1 = temp1;//Why?
                  temp2 = temp1 + len;
 //end crazy MOL stuff
                  std::vector<std::vector<int> > errors=e_obj.checkErrors( i, j, temp1, temp2);

                  std::vector<int>finalErrors=e_obj.getFinalErrors( errors );
//                  float per_err = e_obj.getThreshold(finalErrors,temp1,temp2,ma_snp_ends );
 float per_err = e_obj.getThreshold(finalErrors,temp1,temp2);//overload
                  m_errors.push_back( per_err );
                 if( holdOut  )
                 {
                        float oppHom = ( e_obj.getOppHomThreshold( i, j, temp1, temp2 ) ) / ( temp2 -temp1 );
                        m_holdOutErrors.push_back( oppHom );
                 }
          }
       }  
   }
   std::vector<float>maxes;
   float cutoff = empirical_ma_threshold;
   if(empirical_ma_threshold < 0.0){
   maxes = e_obj.getMaxAverages();
   std::sort(maxes.begin(),maxes.end());
   e_obj.setMaxAverage(maxes);
   cutoff = e_obj.getXthPercentile(ma_threshold); 
   }
   e_obj.setCutoff(cutoff);//set the actual threshold to be used when calculating MA in all other SH
   //
   std::sort( m_errors.begin(), m_errors.end() );
   std::sort( m_holdOutErrors.begin(), m_holdOutErrors.end() );
   std::string str =  " \n No of elements in error check are: "
                      + NumberToString( m_errors.size() );
   str  = str + " \n No of elements in hold  error check are: "
                      + NumberToString( m_holdOutErrors.size() );

        e_obj.log( str );

}


float Consolidator::getHoldOutThreshold( float threshold )
{
 int pos = ( m_holdOutErrors.size() ) * threshold;
 if( pos >= m_holdOutErrors.size() || pos < 0 )
 {
	 std::cerr<< "wrong threshold value. it should be between 0-99.99";
    exit( -1 );
 }
 return m_holdOutErrors[ pos ];  
}
float Consolidator::getPctErrThreshold( float threshold)
{
   int pos = ( m_errors.size()  ) * threshold;
 if( pos >= m_errors.size() || pos < 0 )
 {
	 std::cerr<< "wrong threshold value. it should be between 0-99.99"
        << "errors list size is " << m_errors.size()<<std::endl ;
    exit( -1 );
 }
 return m_errors[ pos ];
}
/*Weighted output functions*/
void Consolidator::update_genome(int snp1,int snp2){
  int s1,s2;
  s1 = snp1;
  s2 = snp2;
  for(int i = s1; i <= s2; i++){
    genome_vector[i] += 1;
  }
}

void Consolidator::update_genome(int snp, float weight){
  if( (snp < 0) || (snp > genome_vector.size()) ){
	  std::cerr << "Error in updating genome. SNP: " << snp << " is outside of genome length." << std::endl;
    exit(-1);
  }
  genome_vector[snp] = weight;
}

void Consolidator::print_genome(){
	std::cout << "Printing Geneome..." << std::endl;
  for(int i = 0; i < genome_vector.size(); i++){
	  std::cout << genome_vector[i] << ", ";
  }
}
float Consolidator::update_snp_weight(int snp1,int snp2){
  int length = ((snp2 - snp1) + 1); //counts are inclusive, so [10,20] has 11 elements, not 10.
  float sum = 0.0;
  float weight = 0.0;
  for(int i = snp1; i <= snp2; i++){
    sum += genome_vector[i];
  }
  weight = sum / length;
  return weight;
}
int Consolidator::find_genome_min(){
  int min = m_weighted_sh[0].snp1;
  for(int i = 1; i < m_weighted_sh.size(); i++){
    if(m_weighted_sh[i].snp1 < min){
      min = m_weighted_sh[i].snp1;
    }
  }
  return min;
}
int Consolidator::find_genome_max(){
  int max = m_weighted_sh[0].snp2;
  for(int i = 1; i < m_weighted_sh.size(); i++){
    if(m_weighted_sh[i].snp2 > max){
      max = m_weighted_sh[i].snp2;
    }
  }
  return max;
}
float Consolidator::average_snp_count(){
  float sum = 0.0;
  float avg = 0.0;
  int zeros = 0;
  for(int i = 0; i < genome_vector.size(); i++){
    if(genome_vector[i] == 0){
      zeros += 1;
    } else {
      sum += genome_vector[i];
    }
  }
  avg = sum / ( genome_vector.size() - zeros);
  return avg;
}

bool Consolidator::isUserSuppliedWeights(){
  return user_supplied_snp_weights.size();
}
/**/
float Consolidator::get_snps_over_range(int snp1, int snp2, float weight){
  int length = ((snp2 - snp1) + 1); //counts are inclusive, so [10,20] has 11 elements, not 10.
  float sum = 0.0;
  for(int i = snp1; i <= snp2; i++){
    sum += genome_vector[i] / weight;
  }
  return sum;
}
