//#include "ErrorCalculator.hpp"
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <exception>
#include <sstream>
#include <fstream>
#include <time.h>
#include <numeric>
#include <algorithm>
#include <assert.h>
#include <cctype>
#include <cstring>
#include <boost/tokenizer.hpp>
#include "ErrorCalculator.hpp"
#include "ErrorFinderManager.hpp"
#include "Consolidator.hpp"
#include "GERMLINE.h"
#include "Ped.hpp"


//using namespace std;
void ErrorCalculator::createLogFile( std::string path )
{
     path  = path + std::string( ".log" );
     m_logger.open( path.c_str(), std::ofstream::out );
     if( !m_logger )
     {
       std::cerr<<"unable to create log file, default out file is FISH.log"<< std::endl;
        m_logger.open( "FISH.log", std::ofstream::out  );
     }
}
void ErrorCalculator::log( std::string& str )
{
    m_logger<< str<<std::endl;
}
void ErrorCalculator::readBmidFile(std::string path)
{
  std::stringstream ss;
  std::string item,line;
     try
     {
       std::ifstream file_bmid(path.c_str());
       if( !file_bmid )
       {
         std::cerr<< "bmid file cannot be oppened. exiting program" <<std::endl;
            exit( -1 );
       }
       Marker cur_marker;
        while ( getline(file_bmid , line) )
        {
                ss.clear(); ss.str( line );

                        ss >> cur_marker.chr >> cur_marker.rsid >> cur_marker.cm_distance >> cur_marker.bp_distance;

                        marker_id.push_back( cur_marker );
       }
      file_bmid.close();
     }
     catch(std::exception &e)
     {
       std::cerr<<e.what()<<std::endl;
      exit( -1 );
     }




}
void ErrorCalculator::changeMapFile( std::string path )
{
     mapper.resize( marker_id.size() );
     std::stringstream ss;
     std::string item,line;
     try
     {
       std::ifstream file_bmid(path.c_str());
       if( !file_bmid )
       {
         std::cerr<< " hold out map file cannot be oppened. exiting program" <<std::endl;
            exit( -1 );
       }

       Marker cur_marker;
       int oPos = 0, nPos = 0;
        while ( getline(file_bmid , line) )
        {
             ss.clear(); ss.str( line );
             ss >> cur_marker.chr >> cur_marker.rsid >> cur_marker.cm_distance >> cur_marker.bp_distance;
             hMarker_id.push_back( cur_marker );
             while( oPos < marker_id.size() &&
                    marker_id[oPos].bp_distance <= hMarker_id[nPos].bp_distance)
       {
                  mapper[ oPos++ ] = nPos;
             }
             ++nPos;
       }
       --nPos;
       while( oPos < marker_id.size() )
       {
           mapper[ oPos++ ] = nPos;
       }
      file_bmid.close();
      if( oPos != marker_id.size() )
      {
          std::cerr<< "the new map file is missing some snps from the old snps which is invalid, match funtion may not work correctly" <<std::endl;

      }
     }
     catch(std::exception &e)
     {
       std::cerr<<e.what()<<std::endl;
      exit( -1 );
     }

}
void ErrorCalculator::readBsidFile(std::string path)
{
   try
   {
     std::ifstream file_bsid(path.c_str() );
       if( !file_bsid )
       {
         std::cerr<< "bsid file cannot be oppened. exiting program" <<std::endl;
            exit( -1 );
       }

       std::stringstream ss;
       std::string item,line;
       while( getline(file_bsid , line) )
        {
                ss.clear();ss.str(line);
                while(getline(ss,item,' '))
                {
                        sample_id.push_back( item );
                        pers_count++;
                }
        }
        pers_count/=2;
        file_bsid.close();


  }
  catch(std::exception &e)
  {
    std::cerr<<e.what()<<std::endl;
     exit( -1 );
  }


}
void ErrorCalculator::readPedFile(std::string path, std::string missing)
{
         if(pers_count<=0)
         {
           std::cerr<<"read bsid file first"<<std::endl;
             return;
         }
         ped_file.clear();
         ped_file.resize(pers_count);
         std::vector < char > locator;
         int start=0, j=0;
         bool warn = false;
         try
         {
           std::ifstream file_ped(path.c_str());
                if( !file_ped )
                {
                  std::cerr<< "original ped file cannot be opened. exiting program" <<std::endl;
                     exit( -1 );
                }

                std::string item,line;
                while(getline(file_ped,line))
                {
                  std::vector<std::string> v;
                  std::string item;
                  std::stringstream ss1(line);
          while(getline(ss1,item,' '))
          {
            v.push_back(item);
          }

                        if(start == 0)
          {
              for(int t = 6; ( t + 1 ) < v.size(); t+=2 )
                {
                   locator.push_back( v[t][0] );
                }
              ++start;
          }
                         for(int t = 6, i = 0; ( t + 1 ) < v.size() && i < locator.size(); t+=2, ++i )
              {
                           SNPs_lrf s;
                           if( missing.compare( v[t] )==0 ||missing.compare( v[t+1] )==0 )
                {
                  warn = true;
                  s.SNP1 = 0;
                  s.SNP2 = 1;
                }
              s.SNP1 = (locator[ i ] == v[t][0] ) ? false : true;
              s.SNP2 = ( locator[ i ] == v[t+1][0] ) ? false : true;
              ped_file[j].push_back(s);
              }
                        j++;
                }

         if(warn)
               {
           std::cerr<<"warning: In ped file "<<path<< ": missing values of type: " << missing << " characters found" << std::endl;
               }
          }
         catch(std::exception &e)
         {
           std::cerr<<e.what()<<std::endl;
                 exit( -1 );
         }
}

void ErrorCalculator::convertPedtovec(std::string path)
{
  std::string line;
  std::string individual;
  std::string dnaseq;
  std::stringstream ss;
  Ped populate;

  std::ifstream file_ped(path.c_str());
  if (!file_ped.is_open() )
  {
    std::cerr<<"Problem opening Ped file to load..Exiting"<<std::endl;
    exit(0);
  }

  while(getline(file_ped,line))
    {
      ss<<line;
      ss>>populate.individual>>populate.individual;
      individual = populate.individual;
      ss>>populate.individual>>populate.individual>>populate.individual>>populate.individual;
      getline(ss,dnaseq,'\n');
      dnaseq.erase(std::remove_if(dnaseq.begin(), dnaseq.end(), ::isspace), dnaseq.end());
      //std::cout<<dnaseq<<std::endl;
      PED.push_back(Ped(individual,dnaseq));
      ss.str("");
      ss.clear();
    }
  file_ped.close();
}




//<piyush-- calculate how many number of AGCT seq ar there in the ped file first line>
/*
int returnACGTcountPED(string pedFile)
{
  cout<<"in here"<<endl;
  return 55;

}
*/
void ErrorCalculator::readHPedFile(std::string path, std::string missing)
{
         if(pers_count<=0)
         {
           std::cerr<<"read bsid file first"<<std::endl;
             return;
         }
         hped_file.clear();
         hped_file.resize(pers_count);
         std::vector < char > locator;
         int start=0, j=0;
         bool warn = false;
         try
         {
           std::ifstream file_ped(path.c_str());
                if( !file_ped )
                {
                  std::cerr<< "hold out file cannot be oppened. exiting program" <<std::endl;
                    exit( -1 );
                 }

                std::string item,line;
                while(getline(file_ped,line))
                {
                  std::vector<std::string> v;
                  std::string item;
                  std::stringstream ss1(line);
                        while(getline(ss1,item,' '))
                        {
                                v.push_back(item);
                        }
                        std::vector< std::string >::iterator it;
                       it = find( sample_id.begin(), sample_id.end(), v[0] );

                        if( it >= sample_id.end() )
                        {
                          std::string str = "This person is missing in original "
                                 " ped file: " + v[0];
                             log( str );
                             continue;
                        }
                        j = ( it - sample_id.begin() ) /2 ;
                        if(start == 0)
      {
                              for(int t = 6; ( t + 1 ) < v.size(); t+=2 )
                              {
                                   locator.push_back( v[t][0] );
                              }
                              ++start;
      }
                         for(int t = 6, i = 0; ( t + 1 ) < v.size() && i < locator.size(); t+=2, ++i )
                         {
                           SNPs_lrf s;
                                if( missing.compare( v[t] )==0 ||missing.compare( v[t+1] )==0 )
        {
                                        warn = true;
          s.SNP1 = 0;
                                        s.SNP2 = 1;
        }

                                s.SNP1 = (locator[ i ] == v[t][0] ) ? false : true;
                                s.SNP2 = ( locator[ i ] == v[t+1][0] ) ? false : true;
                                hped_file[j].push_back(s);
                        }
                      //  j++;
                }

         if(warn)
               {
           std::cerr<<"warning: In ped file "<<path<< ": missing values of type: " << missing << " characters found" << std::endl;
               }
          }
         catch(std::exception &e)
         {
           std::cerr<<e.what()<<std::endl;
                 exit( -1 );
         }
}
//New function for signaling whether or not an SH is an initial drop
bool ErrorCalculator::isInitialCmDrop(int snp1, int snp2, float minLength){
  if( ( (marker_id[snp2].cm_distance - marker_id[snp1].cm_distance) < minLength )  ){
    return true;
  } else {
    return false;
  }
}


//New algorithm to start at middle, and step outwards for trim positions:
//Nate 2/4/2014
  std::vector<int> ErrorCalculator::getTrimPositions(std::vector<float>averages,int snp1,int snp2,float threshold,float minLength){
  std::vector<int> trimPositions;
  int length = snp2 - snp1;
  int midpoint = (int)((length/2.0)+0.5);
  int distance_from_midpoint = 30; //this value by default. Will change if SH.length < 61
  int start = midpoint-distance_from_midpoint, end = midpoint+distance_from_midpoint;

  //debug items
  std::vector<int> starts;
  std::vector<int> ends;
  starts.push_back(start);
  ends.push_back(end);

  //perform a quick initial check on the SH. If it is below the initial SH cM length, then let it be known
  if( ( (marker_id[snp2].cm_distance - marker_id[snp1].cm_distance) < minLength )  ){

    trimPositions.push_back(snp1);
    trimPositions.push_back(snp2);
    trimPositions.push_back(1); //same as error drop code
    return trimPositions;
  }
  if(length <= 61) {
  //In this case, we can't move +- 30 on either side, so we need to determine this size dynamically, based on the length of the SH.
    distance_from_midpoint = (int)( (0.28 * length) + 0.5 );
    trimPositions.push_back(0);
    trimPositions.push_back(length-1);
    return trimPositions;
  }
  for(int i = (midpoint-distance_from_midpoint); i >= 0; i--){ //find the start of the SH. Originally had a value of 30 hardcoded.
    if(averages[i]>threshold){
      start = i+1;
      starts.push_back(start);
      break; //once you've found a start point, get out of there.
    }
    start--;
    starts.push_back(start);
  }
  for(int i = midpoint+distance_from_midpoint; i < length; i++){ //find the end point
    if(averages[i]>threshold){
      end = i-1;
      ends.push_back(end);
      break;
    }
    end++;
    ends.push_back(end);
  }
  //at this point, you have your start and end positions, so check the cM distance, if it passes, then push it. Otherwise, it will have to be
  //thrown out, so push something very small to fail the SNP threshold
  //cout << "In TRIM, the distance is : " << ((marker_id[(snp1+end)].cm_distance) - (marker_id[(snp1+start)].cm_distance)) << " and the minlength is : " << minLength << endl
  if(start == -1) start = 0;

  if( ((marker_id[(snp1+end)].cm_distance) - (marker_id[(snp1+start)].cm_distance)) >= minLength ){
          trimPositions.push_back(start);
    trimPositions.push_back(end);
  }else{ //in this case, we have failed the cM length constraint, so drop it
    trimPositions.push_back(start);
    trimPositions.push_back(end);
    trimPositions.push_back(2);//eo output code
  }
  return trimPositions;
}

  std::vector<int> ErrorCalculator::getTrimPositions(std::vector<float>averages,int snp1,int snp2, float threshold_start, float threshold_end,float minLength,int ma_err_ends)
{
    std::vector<int> trimPositions;
        int length=snp2-snp1;
        int start=0,end=length;
        if(ma_err_ends>=length)ma_err_ends=length;
        for(int i=0;i<ma_err_ends;++i)
        {
                   if(averages[i]>threshold_start)
                   {
                          start=i;
                   }
        }
        for(int i=length-1;i>(length-ma_err_ends);--i)
        {
                   if(averages[i]>threshold_end)
                   {
                          end=i;
                   }
        }
                        trimPositions.push_back(start);
                        trimPositions.push_back(end);
          return trimPositions;
}

//getTrueMovingAverages() is a slightly differnt version of getMovingAverages, and it is intended to be
//used only when calculating the TrueIBD segments. The main difference comes in where it starts calculating
  std::vector<float>  ErrorCalculator::getTrueMovingAverages(std::vector<int> errors,int snp1,int snp2,int width){
    std::vector<float> averages;
    std::vector<int> e_Places;
  int length = snp2-snp1;
  //move in one window length on both ends
  int start = snp1;// + (width/2);
  int end = snp2;// - (width/2);
  int previous = 0;
  int present;
  averages.resize((length-width));
//  averages.resize(length);
  e_Places.resize(length+(width/2));

  for(int i = 0; i < errors.size(); i++){
    e_Places[errors[i]] = 1;
  }

  for(int i = 0; i<(width); i++){
    averages[0] += e_Places[i];
  }
  present = width;
  for(int i = 1; ((i < averages.size())) ; i++){
    averages[i] = averages[i-1] + (e_Places[present]) - (e_Places[previous]);
    present++;
    previous++;
  }

  for(int i = 0; i < averages.size(); i++){
    averages[i] /=width;
  }
  return averages;
}
//======================
//------------------with new, backwards algorithm
  std::vector<float> ErrorCalculator::getTrueMovingAverages2(std::vector<int> errors,int snp1, int snp2, int width){
    std::vector<float> averages;
    std::vector<int> e_Places;
        int length = snp2-snp1;
  int start = snp1;
  int end = snp2;
  int previous = 0;
  int present;
  averages.resize((length-width));
  e_Places.resize(length+(width/2));
  for(int i = 0; i < errors.size(); i++){
                e_Places[errors[i]] = 1;
        }
        for(int i = 0; i<(width); i++){
                averages[0] += e_Places[i];
        }
        present = width;
  for(int i = 1; ((i < averages.size())) ; i++){
                averages[i] = averages[i-1] + (e_Places[present]) - (e_Places[previous]);
                present++;
                previous++;
        }
        for(int i = 0; i < averages.size(); i++){
                averages[i] /=width;
        }
        return averages;
}

/* MOVING AVERAGEPIYUSH
vector<float>  ErrorCalculator::getMovingAverages(vector<int> errors,int snp1,int snp2,int width, int extendSNP)//snp1,snp2=start,end;<piyush> added the param int EXTENDSNP for calculating moving window avg
{
  if(extendSNP%2!=0)//if inputted extended snp is even, divide the window equally
  {
    cout<<"extendSNP=  "<<extendSNP<<endl;
    cout<<"Value is not odd, program is exiting"<<endl;
    exit(0);
  //
  }

  if(extendSNP==0)
  {
    extendSNP= width; // window equals to the extended window
  }
  vector<float> averages;
  vector<int> e_Places;
  int beg=0;
  int end=0;
  int length= snp2-snp1;


  averages.resize(length+1);//0-254,270



  end=errors.size()-1;//<piyush, a1>'th number 1st refers to the first element //255
  //cout<<"errors[end] ="<<errors[end]<<endl;//255

  e_Places.resize(errors[end]+1);
  //cout<<"errors[errors.size()-1] = "<<errors[errors.size()-1]<<endl; // <piyush, a1>maximum error size index errors.size()-1

//e_Places[errors[errors.size()-1]]=1;//<piyush, a1>just marking the end of the segment as 1,can be omitted
  //not really necessary just for say that the last value of error[smthing] is an error and hence 1
  //can be averted.//maybe keep it we'll see
   for(int i=0;i<errors.size();i++)
     {
     e_Places[errors[i]]=1;
     //cout<<"e_Places[errors[i]] = "<<e_Places[errors[i]]<<" errors[i] = "<<errors[i]<<endl;
     }

//cout<<"e_Places.size()= "<<e_Places.size()<<endl;

//cout<<"averages.resize(length) ="<<averages.size()<<endl;//255
int cBeg=snp1;
int cEnd=snp2;
int iter=0;
float avg=0;
int pBeg=0;
int pEnd=0;
pBeg=cBeg - ( (extendSNP/2)-1);
pEnd=cEnd + (extendSNP/2);
//cout<<"cBeg= "<<cBeg<<"cEnd= "<<cEnd<<endl;
//cout<<"pBeg= "<<pBeg<<"pEnd= "<<pEnd<<endl;
int incr=0;//errors[beg];
int avgC=length;


while(iter<=avgC)
  {
  int i;
  //cout<<"snp1= "<<snp1<<"snp2= "<<snp2<<endl;

      if(extendSNP%2==0)
      {
         for ( i=incr-(extendSNP/2)+1 ;i<=incr+(extendSNP/2);i++)
         {
            avg=i;
          //  cout<<"i= "<<i<<endl;
           if(i<0 ||i>errors[end])
           {
            // cout<<"errors[end]"<<errors[end]<<endl;
             avg=0;
             averages[iter]=averages[iter]+avg;
            //  cout<<"averages[iter]= "<<averages[iter]<<"iter= "<<iter<<endl;

           }
           else
           {
           averages[iter]=averages[iter]+e_Places[i];
        //  cout<<"averages[iter]= "<<averages[iter]<<"iter= "<<iter<<endl;
           }
         }

      }
      else
      {
        cout<<"Extenderror exception" <<endl;
        exit(0);
      }
      //cout<<"i"<<i-1<<" averages["<<iter<<"]= "<< averages[iter]<<" iter= " <<iter<<endl;

  incr++;
  iter++;
  snp1++;
  }
//cout<<"Size of average="<<averages.size();


for (int i=0;i<averages.size();i++){
  //cout<<"average b "<<averages[i]<<endl;
  averages[i]/=extendSNP;
  //cout<<"average a "<<averages[i]<<endl;
  //cout<<"i="<<i<<endl;
}
//cout<<"Size of average="<<averages.size()<<endl;



return averages;
}

*/



/*


        for(int i=0;i<presentWidth;i++)
        {
                averages[0]+=e_Places[i];//e[0]+...+e[4] if the win size = 10 (i.e. for the first half)

        }
  int lastPos = averages.size() - 1;
  e_Places[lastPos] = 1; //account for error at the end of the SH...why this isn't already showing up is unknown
  //we will need to calculate the last position here as well, and the summations will all change
  for(int i = lastPos; i > (lastPos - (width/2)); i--){
    averages[(averages.size()-1)] += e_Places[i];
  }
  present -= 1;
  int midpoint = (int)((length/2.0)+0.5);
  //forwards
  for(int i = 1; i < midpoint; i++){
    if(present<width){
      averages[i] = averages[i-1] + e_Places[++present];
    }else{
      averages[i] = averages[i-1] + (e_Places[present] - e_Places[previous]);
      previous++;
      present++;
    }
  }
  //backwards
  int backwardPresent = lastPos - (width/2);
  if(backwardPresent <= 0){
    cerr << backwardPresent << endl;
  }
  int backwardPrevious = lastPos;
  int avLast = averages.size() - 2;
  for(int i = avLast; i >= midpoint; i--){
    if(backwardPresent > (lastPos - width)){
      averages[i] = averages[i+1] + e_Places[backwardPresent--];//try both ++pre and post++
    }else{
      averages[i] = averages[i+1] + (e_Places[backwardPresent] - e_Places[backwardPrevious]);
      backwardPresent--;
      backwardPrevious--;
    }
  }
  //this can all stay the same
        for(int i=0;i<length;i++)
        {
                if(averages[i]<0)averages[i]=0;
                if(i<width/2)
                {
                        averages[i]/=(width/2+i);
                }

                else if((length-i-1)<width/2)
                {
                        averages[i]/=(width/2+(length-i-1));

                }

                else
                {
                        averages[i]/=width;
                }

        }
        return averages;
}
*/

//------------------
//======================
//adding reference as of 4/25..

  std::vector<float>  ErrorCalculator::getMovingAverages(std::vector<int> errors,int snp1,int snp2,int width, int extendSNP)//snp1,snp2=start,end;<piyush> added the param int EXTENDSNP for calculating moving window avg
{
  //cout<<"getMovingAverages snp1= "<<snp1<<endl;
  //cout<<"getMovingAverages snp2 = "<<snp2<<endl;
  //cout<<"getMovingAverages errors size ="<<errors.size()<<endl;
  std::vector<float> averages;
      int length= snp2-snp1;
      std::vector<int> e_Places;
    averages.resize(length);
    e_Places.resize(length+(width/2));
      e_Places[(e_Places.size()-1)]=1;//this is injecting an error at the last possible position. We may ignore this from now on.

    //fill the IE array(errors) with the error locations listed in e_Places[]
    for(int i=0;i<errors.size();i++)
    {
          //this is to adjust for off-by-one indexing I believe
         if( (errors[i] - 1) < 0){

           }else{
              e_Places[(errors[i]-1)]=1;
           }
    }//end for

    int presentWidth = width/2;
    int previous=0,present=presentWidth;
    for(int i=0;i<presentWidth;i++)
    {
            averages[0]+=e_Places[i];

    }
int lastPos = averages.size() - 1;
e_Places[lastPos] = 1; //account for error at the end of the SH...why this isn't already showing up is unknown
//we will need to calculate the last position here as well, and the summations will all change
for(int i = lastPos; i > (lastPos - (width/2)); i--){
  averages[(averages.size()-1)] += e_Places[i];
}
present -= 1;
int midpoint = (int)((length/2.0)+0.5);
//forwards
for(int i = 1; i < midpoint; i++){
  if(present<width){
    averages[i] = averages[i-1] + e_Places[++present];
  }else{
    averages[i] = averages[i-1] + (e_Places[present] - e_Places[previous]);
    previous++;
    present++;
  }
}
//backwards
int backwardPresent = lastPos - (width/2);
if(backwardPresent <= 0){
  std::cerr << backwardPresent << std::endl;
}
int backwardPrevious = lastPos;
int avLast = averages.size() - 2;
for(int i = avLast; i >= midpoint; i--){
  if(backwardPresent > (lastPos - width)){
    averages[i] = averages[i+1] + e_Places[backwardPresent--];//try both ++pre and post++
  }else{
    averages[i] = averages[i+1] + (e_Places[backwardPresent] - e_Places[backwardPrevious]);
    backwardPresent--;
    backwardPrevious--;
  }
}
//this can all stay the same
    for(int i=0;i<length;i++)
    {
            if(averages[i]<0)averages[i]=0;
            if(i<width/2)
            {
                    averages[i]/=(width/2+i);
            }

            else if((length-i-1)<width/2)
            {
                    averages[i]/=(width/2+(length-i-1));

            }

            else
            {
                    averages[i]/=width;
            }

    }
    return averages;
}

std::vector<int>ErrorCalculator::getFinalErrors(std::vector<std::vector< int> > errors)const
{
    std::vector<int>finalPositions;
        int p1=0,p2=0,p3=0,p4=0;
        int len1=0;
        int maxlen=0,prev=0,present=0;
        for(int i=0;i<4;i++) //little unclear as to what this is doing. Perhaps this is always snp2-(snp1+1), which would be the entire length and hence the longest?
        {
                if(errors[i].size()<2)return errors[i];
        }
        int first_pos=getMax(errors[0][p1],errors[1][p2],errors[2][p3],errors[3][p4]); //get the max distance of the first set, this is the longest match and hence is an...error?

        finalPositions.push_back(first_pos); //here we push back an actual error position. Things may be as simple as always setting finalPositions[0] = 0, just to ensure that there is always an
              //error at position 0 in the SH. Similarily, we could add an error to the end of the SH by setting finalPositions[finalPositions.size()-1] = SH.size() as well
   while(p1<errors[0].size()||p2<errors[1].size()||p3<errors[2].size()||p4<errors[3].size())
   {
        maxlen=0;
        len1=0;
        int epos=finalPositions[finalPositions.size()-1];
        while(p1<errors[0].size()&&errors[0][p1]<=epos)
        {
                p1++;

        }
        if(p1<errors[0].size())
        {
                len1=errors[0][p1]-epos;

                if(maxlen<len1)
                {
                        maxlen=len1;
                        present=errors[0][p1];

                }
        }
  while(p2<errors[1].size()&&errors[1][p2]<=epos)
        {
                p2++;

        }
        if(p2<errors[1].size())
        {
                len1=errors[1][p2]-epos;
                if(maxlen<len1)
                {
                        maxlen=len1;


                        present=errors[1][p2];

                }
        }
        while(p3<errors[2].size()&&errors[2][p3]<=epos)
        {
                p3++;

        }
        if(p3<errors[2].size())
        {
                len1=errors[2][p3]-epos;
                if(maxlen<len1)
                {
                        maxlen=len1;


                        present=errors[2][p3];

                }
        }

        while(p4<errors[3].size()&&errors[3][p4]<=epos)
        {
                p4++;

        }

   if(p4<errors[3].size())
        {
                len1=errors[3][p4]-epos;
                if(maxlen<len1)
                {
                        maxlen=len1;
                        present=errors[3][p4];

                }
        }
        if(present>finalPositions[finalPositions.size()-1])
                finalPositions.push_back(present);

  }

  return finalPositions;
}
std::vector<std::vector<int> >  ErrorCalculator::checkErrors(int pers1,int pers2,int snp1,int snp2)const
{
  std::vector<std::vector<int > > errors;
  std::vector<int> finalErrors;
        errors.resize(4);

        if(pers1>=ped_file.size()||pers2>=ped_file.size())
        {

          std::cerr<<"error occured here pers1 as "<<pers1<<" pers2 as "<<pers2<<std::endl;
                exit(-1);

        }

        else if(snp1>=ped_file[pers1].size()||snp2>=ped_file[pers1].size())
        {
          std::cerr<< "snp1 = " << snp1 << " snp2 = "<< snp2<<std::endl;
          std::cerr<<"something went wrong with bmid file or ped file"<<std::endl;
                exit(-1);
        }

        for(int i=snp1;i<snp2;i++)
        {
            if( countGapError == false )
            {
                if(ped_file[pers1][i].SNP1==ped_file[pers1][i].SNP2&&ped_file[pers2][i].SNP1==ped_file[pers2][i].SNP2&&ped_file[pers1][i].SNP1!=ped_file[pers2][i].SNP1)

                     continue;
           }
                if(ped_file[pers1][i].SNP1!=ped_file[pers2][i].SNP1)
                {
                        errors[0].push_back(i-snp1);
                }

                if(ped_file[pers1][i].SNP1!=ped_file[pers2][i].SNP2)
                {
                        errors[1].push_back(i-snp1);
                }
                if(ped_file[pers1][i].SNP2!=ped_file[pers2][i].SNP1)
                {
                        errors[2].push_back(i-snp1);
                }

                if(ped_file[pers1][i].SNP2!=ped_file[pers2][i].SNP2)
                {
                        errors[3].push_back(i-snp1);
                }
   }
                errors[0].push_back(snp2-snp1+1);
                errors[1].push_back(snp2-snp1+1);
                errors[2].push_back(snp2-snp1+1);
                errors[3].push_back(snp2-snp1+1);
        


        return errors;
}
float ErrorCalculator::getOppHomThreshold( int pers1, int pers2, int snp1Old, int snp2Old ) const
{

  if( hped_file[pers1].size() ==0 || hped_file[pers2].size() ==0 )
  {
      return 0;
  }
  if( snp1Old >= mapper.size() || snp2Old >= mapper.size()  )
  {
    std::cerr<< " In hold out phase wrong snps entered: "
            << snp1Old << " or " <<snp2Old <<std::endl;
        exit( -1 );
  }
  int snp1 = mapper[ snp1Old ];
  int snp2 = mapper[ snp2Old ];
  if(pers1>=hped_file.size()||pers2>=hped_file.size())
  {
    std::cerr<<"error occured here pers1 as "<<pers1<<" pers2 as "<<pers2<<std::endl;
     exit(-1);
  }
  else if(snp1>=hped_file[pers1].size()||snp2>=hped_file[pers1].size()

            || snp1>=hped_file[pers2].size()||snp2>=hped_file[pers2].size() )
  {
    std::cerr<<"something went wrong with bmid file or hped file"<< snp1
      << " " << snp2 << " "<< hped_file[pers1].size()  <<std::endl;
     exit(-1);
  }
  float OHCount = 0;
  for(int i = snp1;i < snp2; ++i)
  {
     if( hped_file[pers1][i].SNP1 == hped_file[pers1][i].SNP2 && hped_file[pers2][i].SNP1 == hped_file[pers2][i].SNP2 && hped_file[pers1][i].SNP1 != hped_file[pers2][i].SNP1 )
     {
        ++OHCount;
     }
  }
  return OHCount;
}
int ErrorCalculator::getMax(int &a,int &b, int &c, int &d)const
{
        int max1, max;
        if(a>b)max1=a;
        else max1=b;
        if(c>d)max=c;
        else max=d;
        if(max1>max) max=max1;
        return max;

}
void ErrorCalculator::finalOutPut(int pers1,int pers2,int snp1,int snp2, float min_cm)
{

        if((marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<min_cm)
        {
                return;
        }
 
        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
          std::cerr<<"There is an error with the ped and bsid files, please check that they are correct."<<std::endl;
               return;

        }

             if(snp1>=marker_id.size()||snp2>=marker_id.size())
             {
               std::cerr<<"something went wrong with the bmid file while outputing. The snp1 is: "<<snp1<<" snp2 is "<<snp2<<" marker id size is "<<marker_id.size()<<std::endl;
                     return;

             }
             //Addressing the use of IID instead of FID by stepping up one index in the master sample_id vector. This has been tested, and will be
             //applied to all outputs. It goes from sample_id[pers1*2] --> sample_id[(pers1*2)+1]

             if (!IBD)
         {
          std::cout<<sample_id[(pers1*2)+1]<<"\t"
          <<sample_id[(pers2*2)+1]<<"\t"
          <<marker_id[snp1].bp_distance<<"\t"
          <<marker_id[snp2].bp_distance<<"\t"
          <<(snp2-snp1)<<"\t"
          <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<std::endl;
         }

             else
         {
          std::cout<<sample_id[(pers1*2)+1]<<"\t"
          <<sample_id[(pers2*2)+1]<<"\t"
          <<marker_id[snp1].bp_distance<<"\t"
          <<marker_id[snp2].bp_distance<<"\t"
          <<(snp2-snp1)<<"\t"
          <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"IBD1"<<std::endl;


          /*long SH_L = marker_id[snp1].bp_distance;  //snp1 relative distance i.e. line number of bmid
          long SH_R = marker_id[snp2].bp_distance;  //snp2 ^^*/
          std::vector<bool> status;
          status = ma_het(pers1, pers2, snp1, snp2,  min_cm);

/*          if (status[0] && status[1]) // status 1 stands for ibd4
            {

              std::cout<<"IBD4"<<std::endl;
              return;
              //check for ibd4
            }
          if (status[0])
            {
              std::cout<<sample_id[(pers1*2)+1]<<"\t"
              <<sample_id[(pers2*2)+1]<<"\t"
              <<marker_id[snp1].bp_distance<<"\t"
              <<marker_id[snp2].bp_distance<<"\t"
              <<(snp2-snp1)<<"\t"
              <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"IBD2"<<std::endl;
                return;
            }*/


          if ((status[0] == false)&&(status[1] == false)  )
            {
              status.clear();
              bool status;
              status = ma_het_nm(pers1, pers2, snp1, snp2,  min_cm);
              /*if (status)
              {
                std::cout<<sample_id[(pers1*2)+1]<<"\t"
                <<sample_id[(pers2*2)+1]<<"\t"
                <<marker_id[snp1].bp_distance<<"\t"
                <<marker_id[snp2].bp_distance<<"\t"
                <<(snp2-snp1)<<"\t"
                <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"IBD2"<<std::endl;
              }*/

            }
         }
}


bool ErrorCalculator::ma_het_nm(int pers1,int pers2,int snp1,int snp2, float min_c)
{
  std::string pers1dnasnp1 = seperateSNP(PED[pers1].dnasequence, 0, PED[pers1].dnasequence.length(),0);
  std::string pers1dnasnp2 = seperateSNP(PED[pers1].dnasequence, 0, PED[pers1].dnasequence.length(),1);
  std::string pers2dnasnp1 = seperateSNP(PED[pers2].dnasequence, 0, PED[pers2].dnasequence.length(),0);
  std::string pers2dnasnp2 = seperateSNP(PED[pers2].dnasequence, 0, PED[pers2].dnasequence.length(),1);
  std::vector <double> indices = returnhighest_nm(pers1dnasnp1, pers1dnasnp2, pers2dnasnp1, pers2dnasnp2, snp1, snp2);  //index[2] has the highest cm dist value

  if (indices[2]< min_c )
  {
    //do nothing
    //return false;
  }
  else
  {
    std::cout<<sample_id[(pers1*2)+1]<<"\t"
    <<sample_id[(pers2*2)+1]<<"\t"
    <<marker_id[indices[0]].bp_distance<<"\t"
    <<marker_id[indices[1]].bp_distance<<"\t"
    <<(indices[1]-indices[0])<<"\t"
    <<indices[2]<<"\t"<<"IBD2"<<std::endl;

    //print it
    return true;
  }




}

std::vector<bool> ErrorCalculator::ma_het(int pers1,int pers2,int snp1,int snp2, float min_cm)  //No NM required. Just send (11,21), (11,22),(12,21),(12,22) to het and test for max of them > threashold
{

  std::string pers1dnasnp1 = seperateSNP(PED[pers1].dnasequence, 0, PED[pers1].dnasequence.length(),0);
  std::string pers1dnasnp2 = seperateSNP(PED[pers1].dnasequence, 0, PED[pers1].dnasequence.length(),1);
  std::string pers2dnasnp1 = seperateSNP(PED[pers2].dnasequence, 0, PED[pers2].dnasequence.length(),0);
  std::string pers2dnasnp2 = seperateSNP(PED[pers2].dnasequence, 0, PED[pers2].dnasequence.length(),1);

  std::vector <double> indices1 = returnhighest_het(pers1dnasnp1,pers1dnasnp2,snp1,snp2) ;//snpindexBegin, snpindexEnd, maxcount index 0,1,2 respectively
  std::vector <double> indices2 = returnhighest_het(pers2dnasnp1,pers2dnasnp2,snp1,snp2); // ^^

  std::vector<bool> status; // 2 values. 0 contains information about,  if ibd2 <withoutNM>, is true or false and  1 contains information about ibd4 is true or false
  //std::vector<double> status; // 2 values. 0 contains information about,  if ibd2 <withoutNM>, is true or false and  1 contains information about ibd4 is true or false

  //std::cerr<<indices1[2]<<std::endl;

  if (max(indices1[2],indices2[2]) >= IBD_THRESHOLD )
    {
      status.push_back(true);
      if ( indices1[2]>=indices2[2]    )
        {
          std::cout<<sample_id[(pers1*2)+1]<<"\t"
          <<sample_id[(pers2*2)+1]<<"\t"
          <<marker_id[indices1[0]].bp_distance<<"\t"
          <<marker_id[indices1[1]].bp_distance<<"\t"
          <<(indices1[1]-indices1[0])<<"\t"
          <<indices1[2]<<"\t"<<"IBD2"<<std::endl;
        }
      else
        {
          std::cout<<sample_id[(pers1*2)+1]<<"\t"
          <<sample_id[(pers2*2)+1]<<"\t"
          <<marker_id[indices2[0]].bp_distance<<"\t"
          <<marker_id[indices2[1]].bp_distance<<"\t"
          <<(indices2[1]-indices2[0])<<"\t"
          <<indices2[2]<<"\t"<<"IBD2"<<std::endl;
        }



      if (min(indices1[2],indices2[2]) >= IBD_THRESHOLD )
        {
          status.push_back(true);
          if ( indices1[2]>=indices2[2]    )
            {
              std::cout<<sample_id[(pers1*2)+1]<<"\t"
              <<sample_id[(pers2*2)+1]<<"\t"
              <<marker_id[indices2[0]].bp_distance<<"\t"
              <<marker_id[indices2[1]].bp_distance<<"\t"
              <<(indices2[1]-indices2[0])<<"\t"
              <<indices2[2]<<"\t"<<"IBD4"<<std::endl;

              std::cout<<sample_id[(pers1*2)+1]<<"\t"
              <<sample_id[(pers2*2)+1]<<"\t"
              <<marker_id[indices2[0]].bp_distance<<"\t"
              <<marker_id[indices2[1]].bp_distance<<"\t"
              <<(indices2[1]-indices2[0])<<"\t"
              <<indices2[2]<<"\t"<<"IBD4"<<std::endl;
            }
          else
            {
              std::cout<<sample_id[(pers1*2)+1]<<"\t"
              <<sample_id[(pers2*2)+1]<<"\t"
              <<marker_id[indices1[0]].bp_distance<<"\t"
              <<marker_id[indices1[1]].bp_distance<<"\t"
              <<(indices1[1]-indices1[0])<<"\t"
              <<indices1[2]<<"\t"<<"IBD4"<<std::endl;

              std::cout<<sample_id[(pers1*2)+1]<<"\t"
              <<sample_id[(pers2*2)+1]<<"\t"
              <<marker_id[indices1[0]].bp_distance<<"\t"
              <<marker_id[indices1[1]].bp_distance<<"\t"
              <<(indices1[1]-indices1[0])<<"\t"
              <<indices1[2]<<"\t"<<"IBD4"<<std::endl;
            }

        }
      else
        {
          status.push_back(false);
        }
    }
  else // if ibd2 <withoutNM>, if it ain't ibd2<withoutNM>, then it ain't ibd4
    {
      status.push_back(false);
      status.push_back(false);
    }

  return status;
}

std::string ErrorCalculator::seperateSNP(std::string dnasequence,int zero,int length,int snp1_or_snp2 ) //snp1_or_snp2 --->0 is snp1, 1 is snp2
{
  std::string concat = "";
  if (snp1_or_snp2 == 0)
    {
      for (int i =  0; i < length ; i=i+2 )
        concat+=dnasequence[i];
    }
  else
    {
      for (int i = 1; i < length  ; i=i+2 )
        concat+=dnasequence[i];
    }
  return concat;
}

std::vector<double> ErrorCalculator::returnhighest_het(std::string s1,std::string s2,int startIndex, int endIndex)
{
  double maxcount = 0.0 ;
  double rcount = 0.0;

  int runningcount = 0;
  int snpindexBegin = 0;
  int snpindexBegin_temp = 0;
  int snpindexEnd = 0;
  int snpindexEnd_temp = 0;
  std::vector <double> indices;//snpindexBegin, snpindexEnd, maxcount

  for (int i = startIndex ; i <=endIndex; i++)
  {
    if ((s1[i]!=s2[i] ) )
      {
        if (rcount > maxcount)
        {
          snpindexBegin = snpindexBegin_temp;
          snpindexEnd = snpindexEnd_temp;
          maxcount = rcount;
        }
          runningcount = 0;
          rcount=0.0;
      }
    else if(s1[i]==s2[i])
      {
        if (runningcount == 0)
        {
          snpindexBegin_temp = i;
        }
        else
        {
          snpindexEnd_temp = i;
          rcount =  (marker_id[snpindexEnd_temp].cm_distance - marker_id[snpindexBegin_temp].cm_distance );
        }
        runningcount++;


/*           std::cerr<<marker_id.size()<<std::endl;
           for (int i = 0;i<marker_id.size();i++)
           {
             std::cerr<<marker_id[i].chr<<"\t"<<marker_id[i].rsid<<"\t"<<marker_id[i].cm_distance<<"\t"<<marker_id[i].bp_distance<<endl;
           }*/
      }
    if (  !((i+1)<=endIndex)  )
      {
        if (rcount > maxcount)
          {
            snpindexBegin = snpindexBegin_temp;
            snpindexEnd = snpindexEnd_temp;
            maxcount = rcount;
          }
      }
  }
  indices.push_back(snpindexBegin);
  indices.push_back(snpindexEnd);
  indices.push_back(maxcount);
  return indices;

  //std::cout<<snpindexBegin<<"\t"<<snpindexEnd<<"\t"<<maxcount<<std::endl;
}
std::vector<double> ErrorCalculator::returnhighest_nm(std::string s11,std::string s12,std::string s21,std::string s22,int startIndex, int endIndex)
{
  double maxcount = 0.0;
  double rcount = 0.0;
  int runningcount = 0;
  int snpindexBegin = 0;
  int snpindexBegin_temp = 0;
  int snpindexEnd = 0;
  int snpindexEnd_temp = 0;
  std::string s1,s2;
  std::vector <double> indices;//snpindexBegin, snpindexEnd, maxcount


  for (int i = startIndex ; i <=endIndex; i++)
  {
    s1 = s11[i]+s12[i];
    s2 = s21[i]+s22[i];
    std::sort(s1.begin(), s1.end());
    std::sort(s2.begin(), s2.end());


    if ((s1!=s2 ) )
      {
        if (rcount > maxcount)
        {
          snpindexBegin = snpindexBegin_temp;
          snpindexEnd = snpindexEnd_temp;
          maxcount = rcount;
        }
          runningcount = 0;
          rcount=0.0;
      }
    else if(s1==s2)
      {
        if (runningcount == 0)
        {
          snpindexBegin_temp = i;
        }
        else
        {
          snpindexEnd_temp = i;
          rcount =  (marker_id[snpindexEnd_temp].cm_distance - marker_id[snpindexBegin_temp].cm_distance );
        }
        runningcount++;
      }
    if (  !((i+1)<=endIndex)  )
      {
        if (rcount > maxcount)
          {
            snpindexBegin = snpindexBegin_temp;
            snpindexEnd = snpindexEnd_temp;
            maxcount = rcount;
          }
      }
  }

  //std::cout<<snpindexBegin<<"\t"<<snpindexEnd<<"\t"<<maxcount<<std::endl;
  indices.push_back(snpindexBegin);
  indices.push_back(snpindexEnd);
  indices.push_back(maxcount);
  return indices;
}
/*std::vector<int> returnhighest(std::string s1,std::string s2,int startIndex, int endIndex)
{
  int maxcount = 0 ;
  int runningcount = 0;
  int snpindexBegin = 0;
  int snpindexBegin_temp = 0;
  int snpindexEnd = 0;
  int snpindexEnd_temp = 0;
  std::vector <int> indices;//snpindexBegin, snpindexEnd, maxcount

  for (int i = startIndex ; i <=endIndex; i++)
  {
    if ((s1[i]!=s2[i] ) )
      {
        if (runningcount > maxcount)
        {
          snpindexBegin = snpindexBegin_temp;
          snpindexEnd = snpindexEnd_temp;
          maxcount = runningcount;
        }

          runningcount = 0;
      }
    else if(s1[i]==s2[i])
      {
        if (runningcount == 0)
        {
          snpindexBegin_temp = i;
        }
        else
        {
          snpindexEnd_temp = i;
        }
        runningcount++;
      }
    if (  !((i+1)<=endIndex)  )
      {
        if (runningcount > maxcount)
          {
            snpindexBegin = snpindexBegin_temp;
            snpindexEnd = snpindexEnd_temp;
            maxcount = runningcount;
          }
      }
  }
  indices.push_back(snpindexBegin);
  indices.push_back(snpindexEnd);
  indices.push_back(maxcount);
  return indices;

  //std::cout<<snpindexBegin<<"\t"<<snpindexEnd<<"\t"<<maxcount<<std::endl;
}*/









void ErrorCalculator::finalErrorsOutput(int pers1, int pers2, int snp1, int snp2, float min_cm, float per_err){
  if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
  {
    std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
    return;
  }

  if(snp1>=marker_id.size()||snp2>=marker_id.size())
  {
    std::cerr<<"some thing went wrong with bmid file while outputing. The snp1 is: "<<snp1<<" snp2 is "<<snp2<<" marker id size is "<<marker_id.size()<<std::endl;
    return;
  }
  std::cout<<sample_id[(pers1*2)+1]<<"\t"
  <<sample_id[(pers2*2)+1]<<"\t"
  <<marker_id[snp1].bp_distance<<"\t"
  <<marker_id[snp2].bp_distance<<"\t"
  <<(snp2-snp1)<<"\t"
  <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"
  << per_err << std::endl;
}

void ErrorCalculator::weightedOutput(int pers1, int pers2, int snp1, int snp2, float weight){
  std::cout<<sample_id[(pers1*2)+1]<<"\t"
  <<sample_id[(pers2*2)+1]<<"\t"
  <<marker_id[snp1].bp_distance<<"\t"
  <<marker_id[snp2].bp_distance<<"\t"
  <<(snp2-snp1)<<"\t"
  <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)
  <<"\t"<<weight<<std::endl;
}

/*template <class T>
void ErrorCalculator::fullPlusDroppedOutput( int pers1,int pers2,int snp1,int snp2,int min_snp, float min_cm, std::vector<T> positions,float pct_err,int reason ){
        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
          std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
    <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<positions.size();++i)
        {
          std::cout<<positions[i]<<"/";
        }
        std::cout<<"\t"<< pct_err;
        std::cout << "\t"<< reason << std::endl;
}
//new error1 output
template <class T>
void ErrorCalculator::middleOutPut(int pers1, int pers2, int snp1, int snp2, int min_snp, float min_cm, std::vector<T> positions, float pct_err, int start, int end){
  if( (snp1==-1) || (snp2==-1) ) return;
  if( (sample_id.size() <= (pers1*2)) || (sample_id.size() <= (pers2*2)) ){
    std::cerr << "Something went wrong with the bsid file while outputting data." << std::endl;
    return;
  }
  std::cout << sample_id[(pers1*2)+1] << "\t" << sample_id[(pers2*2)+1] << "\t" << marker_id[snp1].bp_distance << "\t" << marker_id[snp2].bp_distance << "\t"
       << start << "/" << end << "\t" << (snp2-snp1) << "\t" << (marker_id[snp2].cm_distance - marker_id[snp1].cm_distance) << "\t" << "/";
  for(int i = 0; i < positions.size(); i++){
    std::cout << positions[i] << "/";
  }
  std::cout << "\t" << pct_err << std::endl;
}
//error1 output
template <class T>
void ErrorCalculator::middleOutPut( int pers1,int pers2,int snp1,int snp2, int min_snp, float min_cm,std::vector<T> positions, float pct_err )
{
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }

        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
          std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
    <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<positions.size();++i)
        {
          std::cout<<positions[i]<<"/";
        }
        std::cout<<"\t"<< pct_err <<std::endl;

}
template < class T >
void ErrorCalculator::middleOutPut( int pers1,int pers2,int snp1,int snp2,int min_snp, float min_cm, std::vector<T> positions,float pct_err,std::string reason){
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }
        if((marker_id[snp2].cm_distance - marker_id[snp1].cm_distance)<min_cm || (snp2-snp1) < min_snp )
        {
             return;
        }

        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
          std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
    <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<positions.size();++i)
        {
          std::cout<<positions[i]<<"/";
        }
        std::cout<<"\t"<< pct_err;
        std::cout << "\t"<< reason << std::endl;
}


template < class T >
void ErrorCalculator::middleOutPut( int pers1,int pers2,
                                  int snp1,int snp2, int min_snp,
                                  float min_cm,std::vector< T > positions,
                                  std::vector< int > trims, float pct_err )
{
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }
        if((marker_id[snp2].cm_distance - marker_id[snp1].cm_distance)<min_cm || (snp2-snp1) < min_snp )
        {
          std::cout << "Testing a theory. " << std::endl;
             return;
        }

        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
          std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
                <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<trims.size();++i)
        {
          std::cout<<trims[i]<<"/";
        }
        std::cout<<"\t";
        for(int i=0;i<positions.size();++i)
        {
          std::cout<<positions[i]<<"/";
        }
        std::cout<<"\t"<< pct_err <<std::endl;

}
//this one is for ma1, the one above is I think for ma2
template < class T >
void ErrorCalculator::middleOutPut( int pers1,int pers2,
                                  int snp1,int snp2, int min_snp,
                                  float min_cm,std::vector< T > positions,
                                  std::vector< int > trims, float pct_err,int start,int end )
{
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }
        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
          std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
                <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"
         << pct_err << "\t" << start <<"/" << end << "/";
        for(int i=0;i<trims.size();++i)
        {
          std::cout<<trims[i]<<"/";
        }
        std::cout<<"\t";
        for(int i=0;i<positions.size();++i)
        {
          std::cout<<positions[i]<<"/";
        }

}*/
//==================================================
//---------------------------------
//This will be used for the error output option.
//The inputs to this function are as follows:
//pers1 and pers1 are the ids of the two people
//snp1 and snp2 are the locations of the snps on the chromosome(I think...)
//min_snp is the minimum length in snp for a sh to be dropped(possibly depricated
//min_cm is the minimum length in centiMorgans for a sh to be dropped(possibly depricated)
//positions is a vector of the moving averages at each point
//errors is a vector of where each IE is located.
//pct_err is PIE value for this sh
//start/end Trim is the snp where the SH starts/stops AFTER trimming
//start/end is the snp where the SH starts/stops BEFORE trimming
//reason is an int[0,3] that gives additional info about the SH.
//0 => not dropped
//1 => dropped due to length before trimming
//2 => dropped due to length after trimming
//3 => dropped due to PIE
/*template < class T > void ErrorCalculator::errorOutput(int pers1, int pers2, int snp1, int snp2, int min_snp, float min_cm, std::vector< T > positions,std::vector< int > errors, float pct_err, int startTrim,int endTrim, int start, int end, int reason)
{
  if( (snp1 == -1) || (snp2 == -1) ) return;
  if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
    std::cerr<<"something went wrong with bsid file while outputing"<<std::endl;
               return;
        }
  std::cout << sample_id[(pers1*2)+1] << "\t" << sample_id[(pers2*2)+1] << "\t" << marker_id[snp1].bp_distance << "\t" << marker_id[snp2].bp_distance << "\t"
             << (snp2-snp1) << "\t" << (marker_id[snp2].cm_distance-marker_id[snp1].cm_distance) << "\t" << pct_err << "\t" << start << "/" << end << "\t" << startTrim <<"/" << endTrim << "\t" << reason
       << "\t" << "/";
  for(int i=0;i<errors.size();++i)
        {
    std::cout << errors[i] << "/";
        }
  std::cout << "\t";
        for(int i=0;i<positions.size();++i)
        {
          std::cout << positions[i] << "/";
        }
        std::cout << std::endl;
}
//-----------------------------------------
//==================================================
template < class T >
void ErrorCalculator::middleHoldOutPut( int pers1,
                                    int pers2, int snp1, int snp2,
                                    int min_snp, float min_cm,
                                    std::vector< T > positions, std::vector< int > trims,
                                    float pct_err, int numberofOppHoms,
                                    int length )
{
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }
        if((marker_id[snp2].cm_distance - marker_id[snp1].cm_distance)<min_cm || (snp2-snp1) < min_snp )
        {
             return;
        }

        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
          std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
                <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<trims.size();++i)
        {
          std::cout<<trims[i]<<"/";
        }
        std::cout<<"\t";
        for(int i=0;i<positions.size();++i)
        {
          std::cout<<positions[i]<<"/";
        }
        std::cout<<"\t"<< pct_err << "\t"<< numberofOppHoms << "\t" << length <<std::endl;

}*/
/***********************************************************************************
 *Input: a vector<int> of finalErrors(error positions, indexed from 1), start and end snps snp1 and snp2. ma_err_ends is no longer used
 *Output: A float returning the PIE threshold for this SH
 *Description: This version contains the new method for calculating the numerator of PIE. The denominator is simply the length of the post-trimmed SH, which in this
  case is snp2-snp1. The numerator is calculated by summing up all of the implied errors from the finalErrors vector that fall within the interval (snp1,snp2). Some
  care is taken to adjust the finalErrors vector values by 1, so that they index correctly. This indexing difference is simply a result of things have been handled in this
  system historically, so just roll with it.
 ***********************************************************************************/
float ErrorCalculator::getThreshold(std::vector<int> finalErrors,int snp1,int snp2, int ma_err_ends)
{
           int length = snp2-snp1;
       float sum=0;
           for(int i=0;i<finalErrors.size();i++)
           {
    finalErrors[i] = finalErrors[i] - 1; //this is to adjust to array indexing
    if( (finalErrors[i] > snp1) && (finalErrors[i] < (snp2)))
                {
                     ++sum;
                }
           }
           return sum/length;
}
//this overloaded version is called when dealing with TrulyIBD segments
/*******************************************************************************
 *Input: a vector<int> of finalErrors, start and ends snps snp1 and snp2
 *Output: a float returning the PIE threshold for this SH
 *Description: Same as the above version of getThreshold, except for the fact that trulyIBD segments are not yet trimmed when this function is called
 *******************************************************************************/
float ErrorCalculator::getThreshold(std::vector<int> finalErrors, int snp1, int snp2){
  float sum = 0.0;
  int length = snp2 - snp1;
  for(int i = 0; i < finalErrors.size(); i++){
    if( (finalErrors[i] > 0) && (finalErrors[i] < length)){
      ++sum;
    }
  }
  return sum / length;
}

float ErrorCalculator::getCMDistance(int position)
{
      if(position>=marker_id.size())
             {
        std::cerr<<"wrong cm distance reuqested"<<std::endl;
                       exit( -1 );
             }
       return marker_id[position].cm_distance;

}
int ErrorCalculator::getNewSnp( int snp )
{
   if( snp < mapper.size() )
   {
       return mapper[ snp ];
   }
   return -1;
}
void ErrorCalculator::countGapErrors( bool countErrors )
{
   countGapError = countErrors;
}
/******************************************************************************************************************
Input: float x, the percentile that we are interested in returning
Output: float that is the value at the xth percentile
Description: This is used with the max moving averages of our truly IBD segments.
Once we have gone through every trulyIBD segment and stored its maximum moving average value into the sorted
vector max_averages in the calling function, we then call this function with argument x, which simply indexes into
max_averages and returns the value at that percentile. For example, suppose that max_averages is:
{0,0.1,0.2,0.3,0.5,0.6,0.7,0.8,0.9,1.0}. If we let x = 0.7, then we are requesting the 70th percentile
********************************************************************************************************************/
float ErrorCalculator::getXthPercentile(float x){
  assert((x >= 0.0) && (x<= 1.0));//avoid buffer overruns
  int index = (int) ((x * max_averages.size())+0.5); //round up for now
  index = index - 1; //array indexing
  //just to be safe, make sure that this index is still in range
  assert((index >= 0) && (index < max_averages.size()));
  return max_averages[index];
}
