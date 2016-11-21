#include <fstream>
#include <sstream>
#include<string>
#include<vector>
class SNPs;
class Marker;
class ErrorCalculator
{
         public:
         ErrorCalculator():pers_count(0){}
   
         void createLogFile( std::string path  );
       
         void log( std::string& log );

         void readBmidFile(std::string path);
 
         void readBsidFile(std::string path);

         void readPedFile(std::string path, std::string missing);
         
         void finalOutPut( int pers1,int pers2,
                          int snp1,int snp2,float min_cm );
         template < class T >
          void fullPlusDroppedOutput( int pers1,int pers2,int snp1,int snp2,int min_snp, float min_cm, std::vector<T> positions,float pct_err,int reason );
	
	 template < class T >
	 void middleOutPut( int pers1,int pers2, int snp1, int snp2, int min_snp, float min_cm, std::vector< T > positions, float pct_err, int start, int end);

	 template < class T >
         void middleOutPut( int pers1,int pers2,
                            int snp1,int snp2, int min_snp, 
                            float min_cm, std::vector<T >positions, 
                            float pct_err );
        template < class T >
         void middleOutPut( int pers1,int pers2,int snp1,int snp2, 
                            int min_snp, float min_cm, 
                            std::vector< T >positions, 
                            std::vector< int > trims, float pct_err );
        template < class T >
         void middleOutPut( int pers1,int pers2,int snp1,int snp2,
                            int min_snp, float min_cm,std::vector< T >positions, float pct_err, std::string reason);

	template < class T >
         void middleOutPut( int pers1,int pers2,int snp1,int snp2,
                            int min_snp, float min_cm,
                            std::vector< T >positions,
                            std::vector< int > trims, float pct_err,int start, int end );


        template < class T >
         void middleHoldOutPut( int pers1,int pers2,int snp1,int snp2, 
                            int min_snp, float min_cm, 
                            std::vector< T >positions,
                            std::vector< int > trims,
                            float pct_err, int nSNPS, int length );
	//this will be used for the new "erro" output option
	template < class T > void errorOutput(int p1, int p2, int snp1, int snp2, int min_snp, float min_cm, std::vector< T >positions,std::vector< int > trims, float pct_err, int startTrim, 
					      int endTrim, int start, int end, int reason);
    template < class T > void errorOutput(int p1, int p2, int snp1, int snp2, int min_snp, float min_cm, std::vector< T >positions,std::vector< int > trims, float pct_err, int startTrim, 
                          int endTrim, float pie, float maMAX, float randomPie, float randomMAmax);
         
         std::vector<int>getFinalErrors(
                            std::vector<std::vector< int> > errors )const;       
         std::vector<std::vector <int> >checkErrors(
                                         int pers1,int pers2,
                                         int snp1,int snp2)const;

         int getMax(int &a,int &b,int &c,int &d)const;

          std::vector<int>getTrimPositions(
                                       std::vector<float> averages,
                                       int snp1,int snp2,float threshold,
                                       float minLength );
         
         std::vector<int>getTrimPositions(
                                std::vector<float> averages,
                                int snp1,int snp2, float threshold_start, 
                                float threshold_end,float minLength,
                                int ma_snp_ends);
        

         std::vector<float>getMovingAverages(
                               std::vector<int> errors,
                               int snp1,int snp2,
                               int width);

	 std::vector<float>getTrueMovingAverages(std::vector<int> errors,int snp1,int snp2,int width);
      
	 std::vector<float>getTrueMovingAverages2(std::vector<int> errors,int snp1,int snp2, int width);
 
         float getThreshold(
                          std::vector<int>trimPositions,
                          std::vector<int>finalErrors );

	float getThreshold(std::vector<int> finalErrors,
                               int snp1,int snp2, int ma_err_ends ); 
	float getThreshold(std::vector<int> finalErrors, int snp1, int snp2);

         int getNoOfPersons(){return pers_count;}

         virtual ~ErrorCalculator(){}; 

         float getCMDistance(int position);

	 void changeMapFile( std::string mapfile);
            
         void readHPedFile( std::string pedfile, std::string missing );

         float getOppHomThreshold( int pers1, int pers2, int snp1, int snp2 ) const;
         int getNewSnp ( int );
  
         void countGapErrors( bool );

         std::vector<float>getMaxAverages(){return max_averages;}//new
	 void addMaxAverage(float av){max_averages.push_back(av);}//new
         void setMaxAverage(std::vector<float> av){
		max_averages = av;
	 }
         float getXthPercentile(float x);//new
         void setCutoff(float c){cutoff = c;}
	 float getCutoff(){return cutoff;}
	 bool isInitialCmDrop(int snp1, int snp2, float minLength);
     float newCmLength(int snp1, int snp2, float cut_length);
     float adjustCmLength(int snp, float val);
     int snpFinder(float cm_length,int start_snp, int end_snp, int direction);
         private:
         int pers_count;
         std::vector<std::vector<SNPs > >ped_file;
         std::vector<std::vector<SNPs > >hped_file;
	 std::vector< int > mapper;
         std::vector<Marker>marker_id;
	 std::vector<Marker>hMarker_id;
         std::vector<std::string>sample_id;
         std::vector<float>max_averages; //new nate 2/3/2014
         std::ofstream m_logger;
         bool countGapError;
	 float cutoff;
};    

class SNPs
{
       public:
       bool SNP1,SNP2;

};
class Marker
{
	public:
        short   chr;
        std::string  rsid;
        float   cm_distance;
        long    bp_distance;
};
