#include "ErrorCalculator.hpp"
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <exception>
#include <sstream>
#include <fstream>
#include <time.h>
#include <assert.h>
using namespace std;
void ErrorCalculator::createLogFile( std::string path )
{
     path  = path + std::string( ".log" );
     m_logger.open( path.c_str(), std::ofstream::out );
     if( !m_logger )
     {
        cerr<<"unable to create log file, default out file is FISH.log"<< endl;
        m_logger.open( "FISH.log", std::ofstream::out  ); 
     }
}
void ErrorCalculator::log( std::string& str )
{
    m_logger<< str<<endl;
}
void ErrorCalculator::readBmidFile(string path)
{


     stringstream ss;
     string item,line;
     try
     {
       ifstream file_bmid(path.c_str());
       if( !file_bmid )
       {
            cerr<< "bmid file cannot be oppened. exiting program" <<endl; 
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
     catch(exception &e)
     {
      cerr<<e.what()<<endl;
      exit( -1 );
     }

}
void ErrorCalculator::changeMapFile( string path )
{
     mapper.resize( marker_id.size() );
     stringstream ss;
     string item,line;
     try
     {
       ifstream file_bmid(path.c_str());
       if( !file_bmid )
       {
            cerr<< " hold out map file cannot be oppened. exiting program" <<endl;
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
          std::cerr<< "the new map file is missing some snps from the old snps which is invalid, match funtion may not work correctly" <<endl;

      }
     }
     catch(exception &e)
     {
      cerr<<e.what()<<endl;
      exit( -1 );
     }

}
void ErrorCalculator::readBsidFile(string path)
{
   try
   {
       ifstream file_bsid(path.c_str() );
       if( !file_bsid )
       {
            cerr<< "bsid file cannot be oppened. exiting program" <<endl;
            exit( -1 );
       }

       stringstream ss;
       string item,line;
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
  catch(exception &e)
  {
     cerr<<e.what()<<endl;
     exit( -1 );
  }


}
void ErrorCalculator::readPedFile(string path, string missing)
{
         if(pers_count<=0)
         {
             cerr<<"read bsid file first"<<endl;
             return;
         }
         ped_file.clear();
         ped_file.resize(pers_count);
         vector < char > locator;
         int start=0, j=0;
         bool warn = false;
         try
         {
                ifstream file_ped(path.c_str());
                if( !file_ped )
                {
                    cerr<< "original ped file cannot be oppened. exiting program" <<endl;
                     exit( -1 );
                }

                string item,line;
                while(getline(file_ped,line))
                {
                        vector<string> v;
                        string item;
                        stringstream ss1(line);
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
                                SNPs s;
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
                    cerr<<"warning: In ped file "<<path<< ": missing values of type: " << missing << " characters found" << endl; 
               }
          }
         catch(exception &e)
         {
                 cerr<<e.what()<<endl;
                 exit( -1 );
         }
}

void ErrorCalculator::readHPedFile(string path, string missing)
{
         if(pers_count<=0)
         {
             cerr<<"read bsid file first"<<endl;
             return;
         }
         hped_file.clear();
         hped_file.resize(pers_count);
         vector < char > locator;
         int start=0, j=0;
         bool warn = false;
         try
         {
                ifstream file_ped(path.c_str());
                if( !file_ped )
                {
                    cerr<< "hold out file cannot be oppened. exiting program" <<endl;
                    exit( -1 );
                 }

                string item,line;
                while(getline(file_ped,line))
                {
                        vector<string> v;
                        string item;
                        stringstream ss1(line);
                        while(getline(ss1,item,' '))
                        {
                                v.push_back(item);
                        }
                        vector< string >::iterator it;
                       it = find( sample_id.begin(), sample_id.end(), v[0] );

                        if( it >= sample_id.end() )
                        {
                             string str = "This person is missing in original " 
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
                                SNPs s;
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
                    cerr<<"warning: In ped file "<<path<< ": missing values of type: " << missing << " characters found" << endl; 
               }
          }
         catch(exception &e)
         {
                 cerr<<e.what()<<endl;
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
//new algorithm for GL program only. Used to calculate a new cM distance
float ErrorCalculator::newCmLength(int snp1, int snp2, float cut_length){
  return ( (marker_id[snp2].cm_distance - marker_id[snp1].cm_distance) - cut_length);
}

float ErrorCalculator::adjustCmLength(int snp, float val){
  return ((marker_id[snp].cm_distance) + val);
}
//this algorithm is used by the GL program to find the snp closest to a given cM value.
//start_snp is the old snp location, cm_length is the new cM value, end_snp is an upper bound for the 
//search space, and direction is in [0,1], 0 indicated that we are moving up the SH (ascending),
//1 indicates that we are moving down the SH (descending).
int ErrorCalculator::snpFinder(float cm_length,int start_snp, int end_snp, int direction){
  if(start_snp == end_snp){
    cerr << "Error: When performing cM trim, the start snp was greater than or equal to the end_snp. Try a smaller cM length." << endl;
    exit(1);
  }
  if(direction == 0){ 
    while(start_snp < end_snp){
      if(start_snp < 1){
        start_snp++;
        continue;
      }
      if(cm_length > marker_id[start_snp].cm_distance){
        start_snp++;
        continue;
      } else {
        if( (cm_length - marker_id[start_snp - 1].cm_distance) >= (marker_id[start_snp].cm_distance - cm_length)){
          return start_snp;
        } else {
          return start_snp - 1;
        }
      }
    }//end while
  }else{ //going backwards
    while(end_snp > start_snp){
      if(end_snp > marker_id.size()){
        end_snp--;
        continue;
      }
      if(cm_length < marker_id[end_snp].cm_distance){
        end_snp--;
        continue;
      } else {
        if( (cm_length - marker_id[end_snp].cm_distance) >= (marker_id[end_snp+1].cm_distance - cm_length) ){
          return end_snp+1;
        } else {
          return end_snp;
        }
      }
    }//end while
  }
}//end snpFinder

//New algorithm to start at middle, and step outwards for trim positions:
//Nate 2/4/2014
vector<int> ErrorCalculator::getTrimPositions(std::vector<float>averages,int snp1,int snp2,float threshold,float minLength){
	vector<int> trimPositions;
	int length = snp2 - snp1; 
	int midpoint = (int)((length/2.0)+0.5);
	int distance_from_midpoint = 30; //this value by default. Will change if SH.length < 61
	int start = midpoint-distance_from_midpoint, end = midpoint+distance_from_midpoint;

	//debug items
	vector<int> starts;
	vector<int> ends;
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
		//For now, move out 20% on either side...this should probably be standardized at some point.
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
		//just trying this out as a theory for now. Will need a better solution in the future, even if this idea does work.
		trimPositions.push_back(start);
		trimPositions.push_back(end);
		trimPositions.push_back(2);//eo output code
	}
	return trimPositions;
}

vector<int> ErrorCalculator::getTrimPositions(std::vector<float>averages,int snp1,int snp2, float threshold_start, float threshold_end,float minLength,int ma_err_ends)
{
        vector<int> trimPositions;
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
vector<float>  ErrorCalculator::getTrueMovingAverages(vector<int> errors,int snp1,int snp2,int width){
	vector<float> averages;
	vector<int> e_Places;
	int length = snp2-snp1;
	//move in one window length on both ends
	int start = snp1;// + (width/2);
	int end = snp2;// - (width/2);
	int previous = 0;
	int present;
	averages.resize((length-width));
//	averages.resize(length);
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
vector<float> ErrorCalculator::getTrueMovingAverages2(vector<int> errors,int snp1, int snp2, int width){
	vector<float> averages;
        vector<int> e_Places;
        int length = snp2-snp1;
	int start = snp1;
	int end = snp2;
	int previous = 0;
	int present;
	averages.resize((length-width));
	e_Places.resize(length+(width/2));//not sure on this one
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
//------------------
//======================
//adding reference as of 4/25..
vector<float>  ErrorCalculator::getMovingAverages(vector<int> errors,int snp1,int snp2,int width)
{
        vector<float> averages;
	int length= snp2-snp1;
        vector<int> e_Places;
        averages.resize(length);

        e_Places.resize(length+(width/2));
        //fill the IE array(errors) with the error locations listed in e_Places[]
        for(int i=0;i<errors.size();i++)
        {
              
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
vector<int>ErrorCalculator::getFinalErrors(vector<vector< int> > errors)const
{
        vector<int>finalPositions;
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
vector<vector<int> >  ErrorCalculator::checkErrors(int pers1,int pers2,int snp1,int snp2)const
{
        vector<vector<int > > errors;
        vector<int> finalErrors;
        errors.resize(4);

        if(pers1>=ped_file.size()||pers2>=ped_file.size())
        {

                cerr<<"error occured here pers1 as "<<pers1<<" pers2 as "<<pers2<<endl;
                exit(-1);

        }

        else if(snp1>=ped_file[pers1].size()||snp2>=ped_file[pers1].size())
        {
                cerr<< "snp1 = " << snp1 << " snp2 = "<< snp2<<endl;
                cerr<<"something went wrong with bmid file or ped file"<<endl;
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
        cerr<< " In hold out phase wrong snps entered: " 
            << snp1Old << " or " <<snp2Old <<endl;
        exit( -1 );
  }
  int snp1 = mapper[ snp1Old ];
  int snp2 = mapper[ snp2Old ];
  if(pers1>=hped_file.size()||pers2>=hped_file.size())
  {
     cerr<<"error occured here pers1 as "<<pers1<<" pers2 as "<<pers2<<endl;
     exit(-1);
  }
  else if(snp1>=hped_file[pers1].size()||snp2>=hped_file[pers1].size() 

            || snp1>=hped_file[pers2].size()||snp2>=hped_file[pers2].size() )
  {
     cerr<<"something went wrong with bmid file or hped file"<< snp1
      << " " << snp2 << " "<< hped_file[pers1].size()  <<endl;
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
               cerr<<"some thing went wrong with bsid file while outputing"<<endl;
               return;

        }
	
             if(snp1>=marker_id.size()||snp2>=marker_id.size())
             {
                     cerr<<"some thing went wrong with bmid file while outputing. The snp1 is: "<<snp1<<" snp2 is "<<snp2<<" marker id size is "<<marker_id.size()<<endl;
                     return;

             }
                cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
             <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<endl;
        

}
template <class T>
void ErrorCalculator::fullPlusDroppedOutput( int pers1,int pers2,int snp1,int snp2,int min_snp, float min_cm, vector<T> positions,float pct_err,int reason ){
        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
               cerr<<"some thing went wrong with bsid file while outputing"<<endl;
               return;
        }

   cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
    <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<positions.size();++i)
        {
          cout<<positions[i]<<"/";
        }
        cout<<"\t"<< pct_err;
        cout << "\t"<< reason << endl;
}
//new error1 output
template <class T>
void ErrorCalculator::middleOutPut(int pers1, int pers2, int snp1, int snp2, int min_snp, float min_cm, vector<T> positions, float pct_err, int start, int end){
	if( (snp1==-1) || (snp2==-1) ) return;
	if( (sample_id.size() <= (pers1*2)) || (sample_id.size() <= (pers2*2)) ){
		cerr << "Something went wrong with the bsid file while outputting data." << endl;
		return;
	}
	cout << sample_id[(pers1*2)+1] << "\t" << sample_id[(pers2*2)+1] << "\t" << marker_id[snp1].bp_distance << "\t" << marker_id[snp2].bp_distance << "\t"
	     << start << "/" << end << "\t" << (snp2-snp1) << "\t" << (marker_id[snp2].cm_distance - marker_id[snp1].cm_distance) << "\t" << "/";
	for(int i = 0; i < positions.size(); i++){
		cout << positions[i] << "/";
	}
	cout << "\t" << pct_err << endl;
}
//error1 output
template <class T>
void ErrorCalculator::middleOutPut( int pers1,int pers2,int snp1,int snp2, int min_snp, float min_cm,vector<T> positions, float pct_err )
{
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }

        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
               cerr<<"some thing went wrong with bsid file while outputing"<<endl;
               return;
        }

	 cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
		<<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<positions.size();++i)
        {
        	cout<<positions[i]<<"/";     
        }
        cout<<"\t"<< pct_err <<endl;

}
template < class T >
void ErrorCalculator::middleOutPut( int pers1,int pers2,int snp1,int snp2,int min_snp, float min_cm, vector<T> positions,float pct_err,string reason){
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
               cerr<<"some thing went wrong with bsid file while outputing"<<endl;
               return;
        }

   cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
    <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<positions.size();++i)
        {
          cout<<positions[i]<<"/";     
        }
        cout<<"\t"<< pct_err;
        cout << "\t"<< reason << endl;
}


template < class T >
void ErrorCalculator::middleOutPut( int pers1,int pers2,
                                  int snp1,int snp2, int min_snp, 
                                  float min_cm,vector< T > positions, 
                                  vector< int > trims, float pct_err )
{
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }
        if((marker_id[snp2].cm_distance - marker_id[snp1].cm_distance)<min_cm || (snp2-snp1) < min_snp )
        {
	     cout << "Testing a theory. " << endl;
             return;
        }

        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
               cerr<<"some thing went wrong with bsid file while outputing"<<endl;
               return;
        }

         cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
                <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<trims.size();++i)
        {
                cout<<trims[i]<<"/";
        }
        cout<<"\t";
        for(int i=0;i<positions.size();++i)
        {
                cout<<positions[i]<<"/";
        }
        cout<<"\t"<< pct_err <<endl;

}
//this one is for ma1, the one above is I think for ma2
template < class T >
void ErrorCalculator::middleOutPut( int pers1,int pers2,
                                  int snp1,int snp2, int min_snp,
                                  float min_cm,vector< T > positions,
                                  vector< int > trims, float pct_err,int start,int end )
{
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }

        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
               cerr<<"some thing went wrong with bsid file while outputing"<<endl;
               return;
        }

         cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
                <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"
	       << pct_err << "\t" << start <<"/" << end << "/";
        for(int i=0;i<trims.size();++i)
        {
                cout<<trims[i]<<"/";
        }
        cout<<"\t";
        for(int i=0;i<positions.size();++i)
        {
                cout<<positions[i]<<"/";
        }

}
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
template < class T > void ErrorCalculator::errorOutput(int pers1, int pers2, int snp1, int snp2, int min_snp, float min_cm, std::vector< T > positions,std::vector< int > errors, float pct_err, int startTrim,int endTrim, int start, int end, int reason)
{
	if( (snp1 == -1) || (snp2 == -1) ) return;
	if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
               cerr<<"something went wrong with bsid file while outputing"<<endl;
               return;
        }
	cout << sample_id[(pers1*2)+1] << "\t" << sample_id[(pers2*2)+1] << "\t" << marker_id[snp1].bp_distance << "\t" << marker_id[snp2].bp_distance << "\t"
             << (snp2-snp1) << "\t" << (marker_id[snp2].cm_distance-marker_id[snp1].cm_distance) << "\t" << pct_err << "\t" << start << "/" << end << "\t" << startTrim <<"/" << endTrim << "\t" << reason 
	     << "\t" << "/";
	for(int i=0;i<errors.size();++i)
        {
                cout << errors[i] << "/";
        }
        cout << "\t";
        for(int i=0;i<positions.size();++i)
        {
                cout << positions[i] << "/";
        }
	cout << endl;
}
//-----------------------------------------
//==================================================
//overloaded version for new GL prog 8/8/2014
template < class T > void ErrorCalculator::errorOutput(int pers1, int pers2, int snp1, int snp2, int min_snp, float min_cm, std::vector< T > positions,std::vector< int > errors, float pct_err, int start,int end,float PIE, float MAMAX, float randomPIE, float RANDOMmaMAX )
{
  if( (snp1 == -1) || (snp2 == -1) ) return;
  if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
               cerr<<"something went wrong with bsid file while outputing"<<endl;
               return;
        }
  cout << sample_id[(pers1*2)+1] << "\t" << sample_id[(pers2*2)+1] << "\t" << marker_id[snp1].bp_distance << "\t" << marker_id[snp2].bp_distance << "\t"
             << (snp2-snp1) << "\t" << (marker_id[snp2].cm_distance-marker_id[snp1].cm_distance) << "\t" << start << "/" << end << "\t" << PIE << "\t" << MAMAX << "\t" << randomPIE << "\t" << RANDOMmaMAX 
       << "\t" << "/";
	
	if(errors.size() >= 1){
  		for(int i=0;i<(errors.size() - 1);++i)
        	{
                	cout << errors[i] << "/";
        	}
	}
        cout << "\t";
        for(int i=0;i<positions.size();++i)
        {
                cout << positions[i] << "/";
        }
  cout << endl;
}
//-----------------------------------------
//==================================================
template < class T >
void ErrorCalculator::middleHoldOutPut( int pers1, 
                                    int pers2, int snp1, int snp2, 
                                    int min_snp, float min_cm,
                                    vector< T > positions, vector< int > trims, 
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
               cerr<<"some thing went wrong with bsid file while outputing"<<endl;
               return;
        }

         cout<<sample_id[(pers1*2)+1]<<"\t"
             <<sample_id[(pers2*2)+1]<<"\t"
                <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<trims.size();++i)
        {
                cout<<trims[i]<<"/";
        }
        cout<<"\t";
        for(int i=0;i<positions.size();++i)
        {
                cout<<positions[i]<<"/";
        }
        cout<<"\t"<< pct_err << "\t"<< numberofOppHoms << "\t" << length <<endl;

}
/***********************************************************************************
 *Input: a vector<int> of finalErrors(error positions, indexed from 1), start and end snps snp1 and snp2. ma_err_ends is no longer used
 *Output: A float returning the PIE threshold for this SH
 *Description: This version contains the new method for calculating the numerator of PIE. The denominator is simply the length of the post-trimmed SH, which in this
  case is snp2-snp1. The numerator is calculated by summing up all of the implied errors from the finalErrors vector that fall within the interval (snp1,snp2). Some
  care is taken to adjust the finalErrors vector values by 1, so that they index correctly. This indexing difference is simply a result of things have been handled in this
  system historically, so just roll with it.
 ***********************************************************************************/
float ErrorCalculator::getThreshold(vector<int> finalErrors,int snp1,int snp2, int ma_err_ends)
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
float ErrorCalculator::getThreshold(vector<int> finalErrors, int snp1, int snp2){
	float sum = 0.0;
	int length = snp2 - snp1;
	for(int i = 0; i < finalErrors.size(); i++){
		if( (finalErrors[i] > 0) && (finalErrors[i] < length)){
			++sum;
		}
	}
	return sum / length; //for some reason, the numerator is one too high in the GL program. 8/12/14. as of 8/14/14, that was causing -pie
}

float ErrorCalculator::getCMDistance(int position)
{
      if(position>=marker_id.size())
             {
                       cerr<<"wrong cm distance reuqested"<<endl;
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
