#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
bool REDUCE=false, CM_FLAG=false, SNP_FLAG=false, ALLOW_ERR=false; 
int MIN_SNP=0, all_err=1,count_err[100000]={0};
float MIN_CM=0.0;
bool found=false,looks=false;
using namespace std;

struct Marker
{
	short	chr;
	string	rsid;
	float	cm_distance;
	long	bp_distance;
};

float get_distance( Marker& m_left , Marker& m_right , bool& genetic )
{
	if ( m_left.cm_distance == -1 || m_right.cm_distance == -1 )
	{
		genetic = false;
		return m_right.bp_distance - m_left.bp_distance;
	} else
	{
		genetic = true;
		return m_right.cm_distance - m_left.cm_distance;
	}
}

int main (int argc, char* argv[])
{
	if(argc < 4|| argc>8){
			cerr << "Usage: " << argv[0] << " [BMATCH FILE] [BSID FILE] [BMID FILE] or [BMATCH FILE] [BSID FILE] [BMID FILE] -reduced [min_snp] [min_cm] [all_err per snips 100,1000 etc] " << endl;
			return 0;
	}


		if(strcmp(argv[4],"-reduced")==0)
		{
			 REDUCE=true;

			if(argc==8){
				 CM_FLAG=true;
					 SNP_FLAG=true;
					MIN_SNP=atoi(argv[5]);
					MIN_CM=atof(argv[6]);
					all_err=atoi(argv[7]);
					}
			
				
		}

	 else
		 cerr << "Usage: " << argv[0] << " [BMATCH FILE] [BSID FILE] [BMID FILE] or [BMATCH FILE] [BSID FILE] [BMID FILE] -reduced [max_snp] [max_cm] [no. of snips to allow 1 error]" << endl;
	

	string line, discard;
	ifstream file_bmatch( argv[1] , ios::binary );
	ifstream file_bsid( argv[2] );
	ifstream file_bmid( argv[3] );
	ifstream fout;
	fout.open("output");
	if(!file_bmatch || !file_bsid || !file_bmid ) { cerr << "file could not be opened" << endl; return 0; }

	stringstream ss;
	
	// load samples
	vector< string > sample_id;
	while( getline(file_bsid , line) )
	{
		sample_id.push_back( line );
	}
	file_bsid.close();
	
	// load markers
	vector< Marker > marker_id;
	Marker cur_marker;
	while ( getline(file_bmid , line) )
	{
		ss.clear(); ss.str( line );

			ss >> cur_marker.chr >> cur_marker.rsid >> cur_marker.cm_distance >> cur_marker.bp_distance;

		marker_id.push_back( cur_marker );
	}
	file_bmid.close();
	
	// load matches
	unsigned int pid[2];
	unsigned int sid[2],min1;
	int dif;
	float min2;
	int tempsid1=-1,tempsid2=-1,tempsid3=-1,tempsid4=-1;
	string tempper1="",tempper2="",tempper3="",tempper4="";
	float temp_cm=0;
	bool hom[2] , genetic;
	while ( !file_bmatch.eof() )
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

		min1=sid[1]-sid[0]+1;
		min2=get_distance( marker_id[sid[0]] , marker_id[sid[1]] , genetic ); 

		if(SNP_FLAG&&CM_FLAG)if(min1>MIN_SNP||min2>MIN_CM) continue;
		string sub1=sample_id[ pid[0]],sub2;
		int n1=sub1.length(), n2=sub1.find(" "),n3=sub1.find(".0 "),n4=sub1.find(".1 ");
		if(n2>2)if(n3>0||n4>0) sub1=sub1.substr(n2,n1-n2-2);
			else sub1=sub1.substr(n2,n1-n2);

	
		sub2=sample_id[pid[1]];
		n1=sub2.length(); n2=sub2.find(" ");n3=sub2.find(".0 ");n4=sub2.find(".1 ");
                if(n2>2)if(n3>0||n4>0) sub2=sub2.substr(n2,n1-n2-2);
                        else sub2=sub2.substr(n2,n1-n2);
	
			if(tempper1.compare(sub1)==0&&tempper2.compare(sub2)==0&&sub1!=sub2&&sid[0]<tempsid2&&sid[0]/(all_err)<=count_err[pid[0]])
		{
		
			count_err[pid[0]]++;
			tempsid2=sid[1];
			
			 continue;found=true;looks=true;
		
		}
		if(tempsid1!=-1)
		{	
			min1=tempsid2-tempsid1+1;
			 min2=get_distance( marker_id[tempsid1] , marker_id[tempsid2] , genetic );
                cout
                        <<tempper1 << '\t';
                cout
                        << tempper2 << '\t';
                if(!REDUCE)
                cout
                        << marker_id[tempsid1 ].chr << '\t';
                cout
                        << marker_id[ tempsid1 ].bp_distance << ' '
                        << marker_id[ tempsid2 ].bp_distance << ' ';

        if(!REDUCE){    cout
                        << marker_id[ tempsid1 ].rsid << '\t'
                        << marker_id[ tempsid2 ].rsid << '\t';
                }
                        cout
                        <<min1 << '\t'
                        << setiosflags(ios::fixed) << setprecision(2)
                        <<min2 << '\t';

                if(!REDUCE){
                if ( genetic ) cout << "cM\t"; else cout << "MB\t";
                        cout << dif << '\t';
                        if ( hom[0] ) cout << "1\t"; else cout << "0\t";
                        if ( hom[1] ) cout << "1\t"; else cout << "0\t";
                }
		found=false;
		cout<< endl;	
		}

	/*	*/
		
		tempper1=sub1;
                tempper2=sub2;
                tempsid1=sid[0];
                tempsid2=sid[1];
	//	cout << endl;
	}
	if(tempsid1!=-1)
                {
                        min1=tempsid2-tempsid1+1;
                         min2=get_distance( marker_id[tempsid1] , marker_id[tempsid2] , genetic );
                cout
                        <<tempper1 << '\t';
                cout
                        << tempper2 << '\t';
                if(!REDUCE)
                cout
                        << marker_id[tempsid1 ].chr << '\t';
                cout
                        << marker_id[ tempsid1 ].bp_distance << ' '
                        << marker_id[ tempsid2 ].bp_distance << ' ';

        if(!REDUCE){    cout
                        << marker_id[ tempsid1 ].rsid << '\t'
                        << marker_id[ tempsid2 ].rsid << '\t';
                }
                        cout
                        <<min1 << '\t'
                        << setiosflags(ios::fixed) << setprecision(2)
                        <<min2 << '\t';

                if(!REDUCE){
                if ( genetic ) cout << "cM\t"; else cout << "MB\t";
                        cout << dif << '\t';
                        if ( hom[0] ) cout << "1\t"; else cout << "0\t";
                        if ( hom[1] ) cout << "1\t"; else cout << "0\t";
                }
                found=false;
                cout<< endl;
                }


	fout.close();
		file_bmatch.close();
	



}
