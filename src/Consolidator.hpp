#ifndef CONSOLIDATOR_H
#define CONSOLIDATOR_H
//#include "ErrorCalculator.cpp"
#include "ErrorCalculator.hpp"

class SNP_lrf
{
 public:
 int start,end;

};
class Weighted_SH{
public:
    int snp1, snp2;
    int per1,per2;
    float final_weight;
    float snp_weight;
    float mbp_length;
    float wprime;
    float w2prime;
    Weighted_SH(int s1,int s2,int p1,int p2){
        snp1 = s1;
        snp2 = s2;
        per1 = p1;
        per2 = p2;
    }
};
class Consolidator
{
        public:
	//overloaded version to support new -ma-threshold argument and to depricate -ma-err-threshold-start and end
    //update as of 7/30/14 to incorporate new emp-pie-thresh
	void performTrim(ErrorCalculator& e_obj,int window,int ma_snp_ends,float ma_threshold,int min_snp,float min_cm,float per_err_threshold,std::string option,float hThreshold,bool holdout,float empirical_threshold, float empirical_pie_threshold,int extendSNP);//<piyush> added the param int EXTENDSNP for calculating moving window avg
 
       
	void performConsolidation(ErrorCalculator& ecal, 
              int gap, int min_snp, float min_cm,int extendSnp);

        void performHoldOutTrim( ErrorCalculator& ecal, 
             float threshold, std::string hMiss, std::string );      
 
        void readMatches( std::string path, int pers_count, 
                     ErrorCalculator& eCalculator, int trueSNP, float trueCM, int snipExtend,std::string pedFile  ); //<piyush 1> added int snipExtend for extending snip
        void readUserSuppliedSnpWeights( std::string path );

        void finalOutPut(ErrorCalculator &e,float min_cm, int min_snp )const;
        void findTruePctErrors( ErrorCalculator &e,
                                int ma_snp_ends, bool holdOut,int window,float ma_threshold, float empirical_ma_threshold );//overlaoding to add window argument and threshold as well
        void findTrueSimplePctErrors( ErrorCalculator &e,
                                float PIElength, bool holdOut, int window, float ma_threshold, float empirical_ma_threshold );
        float getPctErrThreshold( float threshold);
        float getHoldOutThreshold( float threshold );
        /*Weighted output functions*/
        void update_genome(int snp1,int snp2);
        void update_genome(int snp, float weight);
        void print_genome();
        bool isUserSuppliedWeights();
        float update_snp_weight(int snp1,int snp2);
        int find_genome_min();
        int find_genome_max();
        float average_snp_count();
        /**/
        float get_snps_over_range(int snp1, int snp2, float weight);

	//new methods to handle calculations of moving averages -nate 2/11/2014
	std::vector < std::vector < std::vector < SNP_lrf > > > getTrueMatches(){
		return m_trueMatches;
	}

	void setTrueMatches(std::vector < std::vector < std::vector < SNP_lrf > > > x){
		m_trueMatches = x;
	}
        private:
       int person_count;

        std::vector < std::vector < std::vector < SNP_lrf > > > m_matches;
        std::vector < std::vector < std::vector < SNP_lrf > > > m_trueMatches;
        std::vector< float > m_errors;
        std::vector< float > m_holdOutErrors;
        //new, for weighting algo
        std::vector < Weighted_SH > m_weighted_sh;
        //this should actually hold the ints for the genome array in the weighted algorithm
        std::vector <int> genome_vector;
        //this vector is if a user supplies a file of snp weights
        std::vector<float> user_supplied_snp_weights;
        void sortMatches();

        static bool compareFunction(SNP_lrf s1, SNP_lrf s2);
	std::string consolidated_str;
	std::string initial_drop_str;
	std::string emp_pie_thresh_str;
	std::string ma_drop_str;
	std::string pie_drop_str;
	std::string emp_ma_thresh_str;
	std::string ibg_str;
	std::string final_sh_str;
	std::string ma_thresh_str;
	int global_initial;
};
#endif // ndef CONSOLIDATOR_H
