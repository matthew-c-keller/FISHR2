#include "Consolidator.hpp"
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
using namespace std;


template <typename T>
  string NumberToString ( T Number )
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }

bool Consolidator::compareFunction(SNP s1, SNP s2)
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

void Consolidator::readMatches(string path,int pers_count, ErrorCalculator& eCalculator, int trueSNP, float trueCM )
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
        ifstream file_bmatch(path.c_str(),ios::binary);
        if( !file_bmatch )
         {
             cerr<<"unable to open the bmatch file, exiting the program" << endl;
             exit( -1 );

         }
   
        while ( !file_bmatch.eof())
        {
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
                      cerr<<"problem with bsid file, check it please"<<endl;
                      return;

                 }
                 SNP snp;
                 snp.start=sid[0];
                 snp.end=sid[1];

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
    catch(exception &e)
    {
       cerr<<"Error:"<<e.what()<<endl;
       exit( -1 );
    }
        
}
void Consolidator::performConsolidation(ErrorCalculator& eCalculator, int gap,int min_snp,float min_cm)
{
         int consolidations = 0, removed = 0;
	 for(int i=0;i<person_count;++i)//for each person
        {
                for(int j=i;j<person_count;++j)//compare with each other person
                {
                        int temp1=-1,temp2=-1;
                         for(int l=0;l<m_matches[i][j].size();++l)//for each match
                         {

                               temp1= m_matches[i][j][l].start;

                                temp2= m_matches[i][j][l].end;

        
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
			      
                               if( ( (temp2-temp1)<min_snp) || ( (eCalculator.getCMDistance(temp2)-eCalculator.getCMDistance(temp1))<min_cm) )
                              {
				       
                                       ++removed;
                              }
			}
		}
	}
        /*new*/
	      global_initial = removed;
        consolidated_str = "Number of Consolidations: " + NumberToString( consolidations );
        initial_drop_str = "Number of matches removed due to initial length: " +  NumberToString( removed );
        /*wen*/
       
       
}
//overloaded version
void Consolidator::performTrim(ErrorCalculator& e_obj,int window,
                               int ma_snp_ends, float ma_threshold,
                               int min_snp,float min_cm,
                               float per_err_threshold, string option,
                               float hThreshold, bool holdOut,float empirical_threshold, float empirical_pie_threshold,float cut_value,float cm_cut_value)
{
  int removed1 =0, removed2 = 0, removed3 = 0, removed4 = 0;
  int not_removed = 0;
  int total_count = global_initial;
  bool wrongOption = false;
  if((cut_value >= 0.0) && (cm_cut_value >= 0.0) ){
    cerr << "ERROR: You have specified both a cut value percentage and a cM cut value, but only one is allowed. Please try again." << endl;
    exit(1);
  }
  //8/14/14
  //Adding header for output
  cout << "id1" << "\t" << "id2" << "\t" << "marker_id1" << "\t" << "marker_id2" << "\t" << "snp_length" << "\t" << "cm_length" << "\t" << "start/end" << "\t" << "pie" << "\t" << "ma_max" << "\t" << "random_pie" << "\t" << "random_ma_max" << "\t" << "errors" << "\t" << "moving_averages" << endl;
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
        if(i == j){
          continue; //ignore same person matches
        }

        int pers1 = i, pers2 = j;
        //cerr << "For this snp, before applying the cut argument, the length is: " << ((m_matches[i][j][l].end - m_matches[i][j][l].start)) << " SNPs" << endl;
        if((m_matches[i][j][l].end - m_matches[i][j][l].start) < min_snp){
          //reduced argument forces an initial drop due to SNPS
          continue;
        }
        //reduced argument check for mincm
        if(e_obj.isInitialCmDrop(m_matches[i][j][l].start,m_matches[i][j][l].end,min_cm)){
          continue;
        }
        
        int temp1,temp2;
        if(cut_value >= 0.0){
        //remove first and last 25% by default, this can also be a user-specified value
        float new_cut_value = (1 - cut_value) * 0.5; //cut it in half for each side
        temp1 = m_matches[ i ][ j ][ l ].start +
        (int)(( m_matches[ i ][ j ][ l ].end -
        m_matches[ i ][ j ][ l ].start ) * new_cut_value);
        temp2 = m_matches[ i ][ j ][ l ].end -
        (int)(( m_matches[ i ][ j ][ l ].end -
        m_matches[ i ][ j ][ l ].start ) * new_cut_value);
        } else {
          //check that cm arg does not exceed the SH's cM length
        float new_cm_length = e_obj.newCmLength(m_matches[i][j][l].start, m_matches[i][j][l].end,cm_cut_value); //change that cut_value parameter.
        new_cm_length = (new_cm_length / 2.0);
        if(new_cm_length <= 0.0){
          cerr << "Error: The cM length that you specified was larger than one or more of the SH lengths." << endl;
          exit(1);
        }
        float new_cm_start = e_obj.adjustCmLength(m_matches[i][j][l].start, new_cm_length);
        float new_cm_end =  e_obj.adjustCmLength(m_matches[i][j][l].end, ((-1.0) * new_cm_length) );
        if(new_cm_start >= new_cm_end){
          cerr << "ERROR: start value in cM greater than end value. Exiting..." << endl;
          exit(1);
        }
        temp1 = e_obj.snpFinder(new_cm_start,m_matches[i][j][l].start,m_matches[i][j][l].end,0);
        temp2 = e_obj.snpFinder(new_cm_end,m_matches[i][j][l].start,m_matches[i][j][l].end,1);
        }
        //we need to calculate for both this person and a random person.
        int randpers1;
        int randpers2;





        //for now, let's enable randomness by default
        randpers1 = std::rand() % e_obj.getNoOfPersons();
        randpers2 = std::rand() % e_obj.getNoOfPersons();
        if( randpers1 > randpers2 )
        {
          randpers1 = randpers1 + randpers2;
          randpers2 = randpers1 - randpers2;
          randpers1 = randpers1 - randpers2;
        }

        //we now have pers1,2 and randpers1,2

        vector<vector<int> > errors=e_obj.checkErrors(pers1, pers2, temp1, temp2);
        vector<int>finalErrors=e_obj.getFinalErrors(errors);
        vector<vector<int> > random_errors = e_obj.checkErrors(randpers1,randpers2,temp1,temp2);
        vector<int>random_final_errors = e_obj.getFinalErrors(random_errors);
        vector<float>movingAverages;
        vector<float>random_moving_averages;
        float threshold;
  
        movingAverages = e_obj.getMovingAverages(finalErrors,temp1,temp2,window);
        random_moving_averages = e_obj.getMovingAverages(random_final_errors,temp1,temp2,window);
        if(empirical_threshold < 0.0){
         threshold = e_obj.getCutoff(); //this is the empirical average threshold for moving averages
        } else {
          threshold = empirical_threshold;
        }
        
        float per_err = e_obj.getThreshold(finalErrors,temp1,temp2);
        float random_per_err = e_obj.getThreshold(random_final_errors,temp1,temp2);
        //find maximum of moving averages
        int max_pos = 0;
        int rmax_pos = 0;
        float max_ma = 0.0;
        float rmax_ma = 0.0;

        for(int z = 1; z < movingAverages.size();z++){
          if(movingAverages[z] > movingAverages[max_pos]){
            max_pos = z;
          }
        }
        max_ma = movingAverages[max_pos];
        for(int c = 1; c < random_moving_averages.size(); c++){
          if(random_moving_averages[c] > random_moving_averages[rmax_pos]){
            rmax_pos = c;
          }
        }
        rmax_ma = random_moving_averages[rmax_pos];

        if( (option.compare("Error1") == 0 ) || (option.compare("ErrorRandom1") == 0) || (option.compare("Error") == 0) ){
          not_removed++;
          e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,per_err,max_ma,random_per_err,rmax_ma);
          continue;
        }//end error1
      }//l
    }//j
  }//i

  ma_drop_str = "No of matches removed due to length of trimming by moving averages: " + NumberToString( removed2 );
  pie_drop_str = "No of matches removed due to percentage error: " + NumberToString( removed1 );

}//endperftrim()

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

        /*What do these two functions do? Is this necessary for being able to find errors, or 
        is it only useful for MA calculations?*/
        vector<vector<int> > trueErrors=e_obj.checkErrors( i, j, t1, t2);
        vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );
        /*x*/

        //This section handles finding the maximum moving averages amongst trulyIBD segments
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

        /*Here they are again.
        */
        vector<vector<int> > errors=e_obj.checkErrors( i, j, temp1, temp2);

        vector<int>finalErrors=e_obj.getFinalErrors( errors );
        //                  float per_err = e_obj.getThreshold(finalErrors,temp1, temp2, 0 );
        float per_err = e_obj.getThreshold(finalErrors,temp1,temp2); //overload!
        m_errors.push_back( per_err );
        /*x*/


      }//end for(l)
    }//end for(j)
  }//end for(i)

    //this section actually handles the sorting of the max averages, and the setting of the user supplied percentile.
    vector<float>maxes;
    float cutoff = empirical_ma_threshold; //assume the user wanted to supply a value. This value will be overwritten shortly if they did not.
    if(empirical_ma_threshold < 0.0){
      maxes = e_obj.getMaxAverages();
      std::sort(maxes.begin(),maxes.end());
      e_obj.setMaxAverage(maxes);
      cutoff = e_obj.getXthPercentile(ma_threshold); //<-make that an actual user input value
    }

    e_obj.setCutoff(cutoff);//set the actual threshold to be used when calculating MA in all other SH

    if(empirical_ma_threshold < 0.0){
      ma_thresh_str = "User supplied ma-threshold is: " + NumberToString(ma_threshold);
      emp_ma_thresh_str = "Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
    } else {
      emp_ma_thresh_str = "Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
  }
  //----------------------------------------
  std::sort( m_errors.begin(), m_errors.end() );
  std::sort( m_holdOutErrors.begin(), m_holdOutErrors.end() );
  ibg_str = "No of segments deemed to be IBD for finding empirical error threshold "
  + NumberToString( m_errors.size() );
}//end ftspe



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
	      vector<vector<int> > trueErrors=e_obj.checkErrors( i, j, t1, t2);
              vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );

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
                  temp1 = temp1;
                  temp2 = temp1 + len;
		  //end crazy MOL stuff
                  vector<vector<int> > errors=e_obj.checkErrors( i, j, temp1, temp2);

                  vector<int>finalErrors=e_obj.getFinalErrors( errors );
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
   vector<float>maxes;
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
    cerr<< "wrong threshold value. it should be between 0-99.99";
    exit( -1 );
 }
 return m_holdOutErrors[ pos ];  
}
float Consolidator::getPctErrThreshold( float threshold)
{
   int pos = ( m_errors.size()  ) * threshold;
 if( pos >= m_errors.size() || pos < 0 )
 {
    cerr<< "wrong threshold value. it should be between 0-99.99"
        << "errors list size is " << m_errors.size()<<endl ;
    exit( -1 );
 }
 return m_errors[ pos ];
}

