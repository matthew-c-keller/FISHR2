

#include <iostream>
//#include <list>
#include <vector>
#include <string.h>
#include <sstream>
#include <fstream>
#include <time.h>
#include <cstdlib>
//#include <dirent.h>
#include <cmath>
//#include <algorithm> // for min(1,2)
#include <random> //for normal rv
//#include <iomanip>      // std::setprecision

//#include "CommFunc.h"


using namespace std;


void ras_show_help(void);
void option_proc(int option_num, char* option_str[]);
void option_do(void);
void ras_read_infile(string GCfile);
long int ras_FileLineNumber(string file_name);
long int ras_FileColNumber(string file_name);
int ras_read_pheno(string file_name);
int ras_read_ped(string file_name);
int ras_read_map(string file_name);
int ras_comp_allele_freq(void);



// haps to ped file
int ras_read_shapeit_haps(string file_name);
int ras_read_shapeit_sample(string file_name);
int ras_convert_shapeit_haps_to_ped(void);

vector<string> split(string str,string sep);




// convert .haps to .ped
bool flag_haps2ped=false;
bool flag_code01=false;
vector< vector<int> > matrix;
long int n_col=0;
vector< vector<string> > matrix_shapeit_sample;
vector< vector<string> > matrix_shapeit_allele;


// help
bool flag_Help=false;

// version
string GLOB_VER=" 0.2 "; //length should be 5
bool flag_version=false;
void ras_show_version(void);


// input
bool flag_InputFile=false;
string Inputfile="gap";
bool flag_InputFile2=false;
string Inputfile2="gap";

int nskip=0;
string sep=" ";
long int n_SNP = 0;
long int n_IND = 0;
long int n_SNP2 = 0;
long int n_IND2 = 0;
double GLOB_MAF=0;
double GLOB_THR=0;
long int GLOB_NROW = 0;
long int GLOB_NCOL = 0;

//ped file
vector< vector<string> > matrix_ped;
vector< vector<string> > matrix_ped_Allele;
vector< vector<long> > matrix_ped_Allele_frq;

//map file
vector< vector<string> > matrix_map_chr_rs;
vector< double > matrix_map_cm;
vector< long > matrix_map_bp;

//dat file
//vector< vector<string> > matrix_dat_str;
vector< vector<double> > matrix_dat;
vector< vector<double> > matrix_dat_testset;
int ras_read_dat_int(string file_name);
int ras_read_dat_int2(string file_name);
//int ras_read_dat_str(string file_name);


//other
bool flag_dirGC=false;
string dirGC="";
//output
string outfile="out";
// phenotypes
bool flag_PhenoFile=false;
string PhenoFile="pheno";
vector< vector<double> > matrix_pheno;
long int n_Pheno = 0;

//string logfile="";




// DataMan
bool flag_DataMan=false;
int DataMan_OE_col=0; //0=deactive, 1=odd, 2=even
int DataMan_OE_row=0; //0=deactive, 1=odd, 2=even
bool flag_transpose=false;







class RasDataMan{
public:
    vector< vector<string> > matrix_dat_str;
    long int nrow=0;
    long int ncol=0;
    int DataMan_OE_col=0; //0=deactive, 1=odd, 2=even
    int DataMan_OE_row=0; //0=deactive, 1=odd, 2=even
    string infile="";
    string outfile="";
    bool flag_transpose=false;
    
public:
    int main_DataMan(void);
private:
    int ras_selcols(void);
    int ras_selrows(void);
    int ras_read_dat_str(string file_name);
    int write_transpose(void);
};





//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

int main(int argc, char *argv[])
{
    cout<<endl;
    cout<<"********************************************************************************"<<endl;
    cout<<"*                        GAP (Genetic Analysis Program)                        *"<<endl;
    cout<<"*                                 Version " << GLOB_VER << "                                *"<<endl;
    cout<<"*                (C) 2014 Rasool Tahmasbi and Matthew C. Keller                *"<<endl;
    cout<<"*                    Institute for Behavioral Genetics (IBG)                   *"<<endl;
    cout<<"*                       University of Colorado at Boulder                      *"<<endl;
    cout<<"********************************************************************************"<<endl;
    cout<<endl;

    long int time_used=0, start=time(NULL);
    time_t curr=time(0);
    cout << "Analysis started: " << ctime(&curr) << endl;

    try{
        option_proc(argc, argv);
        option_do();
    }
    catch (string& s){
        cout << "caught it: \"" << s << "\"" << endl;
    }

    cout << endl << endl;
    curr=time(0);
    cout << "Analysis finished: " << ctime(&curr);
    time_used=time(NULL)-start;
    cout<<"Computational time: "<<time_used/3600<<":"<<(time_used%3600)/60<<":"<<time_used%60<<endl<<endl;

    return 0;
}


void option_proc(int option_num, char* option_str[])
{
    int i=0;
    int argc=option_num;
    vector<char *> argv(option_num+2);
    for(i=0; i<option_num; i++) argv[i]=option_str[i];
    argv[option_num]=(char *) "temp";
    argv[option_num+1]=(char *) "temp";
    
    cout << "Command and Options:" << endl;
    for (i=1; i<argc; i++){
        //=====================================================================
        // general
        if(strcmp(argv[i],"--file")==0){
            flag_InputFile=true;
            Inputfile=(string)argv[++i];
            cout<<"  --file "<<argv[i]<<endl;
        }
        if(strcmp(argv[i],"--file2")==0){
            flag_InputFile2=true;
            Inputfile2=(string)argv[++i];
            cout<<"  --file2 "<<argv[i]<<endl;
        }
        if(strcmp(argv[i],"--nskip")==0){
            nskip=atoi(argv[++i]);
            cout<<"  --nskip "<<argv[i]<<endl;
        }
        if(strcmp(argv[i],"--sep")==0){
            sep=(string) argv[++i];
            if(sep[0]=='\\' && sep[1]=='t') sep="\t";
            cout<<"  --sep "<<argv[i]<<endl;
        }
        if(strcmp(argv[i],"--out")==0){
            outfile=argv[++i];
            cout<<"  --out "<< argv[i] <<endl;
        }
        if(strcmp(argv[i],"--pheno")==0){
            flag_PhenoFile=true;
            PhenoFile=argv[++i];
            cout<<"  --pheno " << argv[i] <<endl;
        }
        if(strcmp(argv[i],"--maf")==0){
            GLOB_MAF=atof(argv[++i]);
            cout << "  --maf " << argv[i] <<endl;
        }
        if(strcmp(argv[i],"--thr")==0){
            GLOB_THR=atof(argv[++i]);
            cout << "  --thr " << argv[i] <<endl;
        }
        
        //=====================================================================
        // -haps2ped
        if(strcmp(argv[i],"--haps2ped")==0){
            flag_haps2ped=true;
            cout<<"  --haps2ped (Command)"<<endl;
        }   
        if(strcmp(argv[i],"--code01")==0){
            flag_code01=true;
            cout<<"  --code01"<<endl;
        }   
        //=====================================================================
        // -help
        if(strcmp(argv[i],"--help")==0 || strcmp(argv[i],"-h")==0 || strcmp(argv[i],"?")==0){
            flag_Help=true;
            cout<<"  " << argv[i] << " (Command)" <<endl;
        }
        //=====================================================================
        // -version
        if(strcmp(argv[i],"--vesion")==0 || strcmp(argv[i],"-v")==0){
            flag_version=true;
            cout<<"  " << argv[i] << " (Command)" <<endl;
        }
        //=====================================================================
        // --dataman
        if(strcmp(argv[i],"--dataman")==0){
            flag_DataMan=true;
            cout<<"  --dataman (Command)"<<endl;
        }
        if(strcmp(argv[i],"--selcol")==0){
            cout << "  " << argv[i];
            DataMan_OE_col=atoi(argv[++i]);
            cout << " " << argv[i] <<endl;
        }
        if(strcmp(argv[i],"--selrow")==0){
            cout << "  " << argv[i];
            DataMan_OE_row=atoi(argv[++i]);
            cout << " " << argv[i] <<endl;
        }
        if(strcmp(argv[i],"--transpose")==0){
            cout << "  " << argv[i] << endl;
            flag_transpose=true;
        }
    }
    cout << endl;
}




void option_do(void)
{
    // --------------------------------------------------------------------
    // help
    if (flag_Help){
        ras_show_help();
    }
    // --------------------------------------------------------------------
    // version
    if (flag_version){
        ras_show_version();
    }
    // --------------------------------------------------------------------
    // To convert .haps to .ped: '-haps2ped -file [file_name] -out [file_out]
    else if (flag_haps2ped){
        if (flag_InputFile){
            ras_read_shapeit_haps(Inputfile+".haps");
            ras_read_shapeit_sample(Inputfile+".sample");
            ras_convert_shapeit_haps_to_ped();
        }
        else{
            throw("Error: No input file.");
        }
    }
    
    // --------------------------------------------------------------------
    // RasDataMan
    else if (flag_DataMan){
        //if (flag_InputFile){
        //    ras_read_dat_str(Inputfile);
        //}
        //else{
        //    throw("Error: No input file.");
        //}
        
        
        
        RasDataMan ras_Dataman;
        //ras_Dataman.matrix_dat_str=matrix_dat_str;
        //ras_Dataman.nrow=GLOB_NROW;
        //ras_Dataman.ncol=GLOB_NCOL;
        ras_Dataman.DataMan_OE_col=DataMan_OE_col;
        ras_Dataman.DataMan_OE_row=DataMan_OE_row;
        ras_Dataman.infile=Inputfile;
        ras_Dataman.outfile=outfile;
        ras_Dataman.flag_transpose=flag_transpose;
        // run main func
        ras_Dataman.main_DataMan();
    }
    
    
    // --------------------------------------------------------------------
    // no command
    else{
        cout<<"No COMMAND is used!"<<endl << endl;
        ras_show_help();
    }
}


void ras_read_infile(string GCfile)
{
    long int number_of_lines = 0;
    ifstream ifile(GCfile.c_str());
    if(!ifile) throw("Error: can not open the file ["+GCfile+"] to read.");
    cout<<"Reading GenCall file from ["+GCfile+"]."<<endl;
    string line;
    
    while (getline(ifile, line)){
        number_of_lines++;
    }
    cout<< "Lines: " << number_of_lines <<endl;
    
    ifile.clear();
    ifile.close();
}



vector<string> split(string str,string sep)
{
    char* cstr=const_cast<char*>(str.c_str());
    char* current;
    vector<string> arr;
    current=strtok(cstr,sep.c_str());
    while(current!=NULL){
        arr.push_back(current);
        current=strtok(NULL,sep.c_str());
    }
    return arr;
}



long int ras_FileLineNumber(string file_name)
{
    long int number_of_lines = 0;
    ifstream ifile(file_name.c_str());
    if(!ifile) throw("Error: can not open the file ["+file_name+"] to read.");

    string line;
    while (getline(ifile, line)){
        number_of_lines++;
    }    
    ifile.clear();
    ifile.close();
    return number_of_lines;
}



long int ras_FileColNumber(string file_name)
{
    long int number_of_cols = 0;
    ifstream ifile(file_name.c_str());
    if(!ifile) throw("Error: can not open the file ["+file_name+"] to read.");

    vector<string> st1;
        
    string line;
    if (getline(ifile, line)){
        st1=split(line,sep);
        number_of_cols=st1.size();
    }    
    ifile.clear();
    ifile.close();
    return number_of_cols;
}


int ras_read_shapeit_haps(string file_name)
{
    ifstream ifile(file_name.c_str());
    if(!ifile) throw("Error: can not open the file ["+file_name+"] to read.");

    n_SNP=ras_FileLineNumber(file_name);
    n_col=ras_FileColNumber(file_name);
    
    cout<< "Loading [" << file_name << "] with " << n_SNP << " rows and " << n_col << " cols ..." << endl;
    cout<< "  Number of SNPs: " << n_SNP << endl;
    cout<< "  Number of inds: " << (n_col-5)/2 << endl;
    
    
    matrix.resize(n_SNP, vector<int>(n_col-5, 0) ); //
    matrix_shapeit_allele.resize( n_SNP , vector<string>(2, "") );
    
    vector<string> st1;

    string line;
    long int j=0,i=0;
    while (getline(ifile, line)){ // for each SNP
        st1=split(line,sep);
        matrix_shapeit_allele[j][0]=st1[3]; // First allele
        matrix_shapeit_allele[j][1]=st1[4]; // Second allele
        for(i=5; i<st1.size(); i++) // 2*inds
        {
            matrix[j][i-5]=atoi(st1[i].c_str());
        }
        j++;
        cout << "\r  " << (j*100/n_SNP) << "% completed." << flush;
    }    
    cout << endl;
    ifile.clear();
    ifile.close();
    
    return 0;
}



int ras_read_shapeit_sample(string file_name)
{
    ifstream ifile(file_name.c_str());
    if(!ifile) throw("Error: can not open the file ["+file_name+"] to read.");

    long int n_r=ras_FileLineNumber(file_name);
    long int n_c=ras_FileColNumber(file_name);
    cout<< "Loading [" << file_name << "] with " << n_r << " rows and " << n_c << " cols ..." << endl;
    cout<< "  Number of inds: " << (n_r-2) << endl;

    matrix_shapeit_sample.resize(n_r-2 , vector<string>(n_c, "") );

    vector<string> st1;

    string line;
    long int j=0,i=0;
    getline(ifile, line); //discard first line
    getline(ifile, line); //discard second line
        
    while (getline(ifile, line)){
        st1=split(line,sep);
        for(i=0; i<st1.size(); i++)
        {
            matrix_shapeit_sample[j][i]=st1[i];
        }
        j++;
        cout << "\r  " << (j*100/(n_r-2)) << "% completed." << flush;
    }    
    cout << endl;
    ifile.clear();
    ifile.close();
    
    return 0;
}



int ras_convert_shapeit_haps_to_ped(void)
{
    cout << "Writing ["+outfile+".ped] ..." << endl;
    ofstream myfile((outfile + ".ped").c_str());
    
    long int i,j;
    long int p=(n_col-5)/2; // p=number of inds
    if (flag_code01){
        for(j=0;j<p;j++)
        {
            //myfile << matrix_shapeit_sample[j][0] << sep << matrix_shapeit_sample[j][1] << sep;
            //myfile << matrix_shapeit_sample[j][3] << sep << matrix_shapeit_sample[j][4] << sep;
            //myfile << matrix_shapeit_sample[j][5] << sep << matrix_shapeit_sample[j][6] << sep;
            myfile << matrix_shapeit_sample[j][0] << sep << matrix_shapeit_sample[j][1] << sep; //FID, IID
            myfile << "0" << sep << "0" << sep; // PID, MID
            myfile << "0" << sep << "-9"; // sex, phenotype
            for(i=0;i<n_SNP;i++)
            {
                myfile << sep << matrix[i][2*j] << sep << matrix[i][2*j+1];
            }
            myfile << endl;
            cout << "\r  " << (j*100/(p-1)) << "% completed." << flush;
        }
    }
    else{
        for(j=0;j<p;j++)
        {
            myfile << matrix_shapeit_sample[j][0] << sep << matrix_shapeit_sample[j][1] << sep;
            myfile << matrix_shapeit_sample[j][3] << sep << matrix_shapeit_sample[j][4] << sep;
            myfile << matrix_shapeit_sample[j][5] << sep << matrix_shapeit_sample[j][6];
            for(i=0;i<n_SNP;i++)
            {
                int al1=matrix[i][2*j];
                int al2=matrix[i][2*j+1];
                myfile << sep << matrix_shapeit_allele[i][al1]  << sep << matrix_shapeit_allele[i][al2];
            }
            myfile << endl;
            cout << "\r  " << (j*100/(p-1)) << "% completed." << flush;
        }
    }
    cout << endl;
    myfile.close();
    
    matrix.clear();
    matrix_shapeit_sample.clear();
    matrix_shapeit_allele.clear();
    
    return 0;
}

//////////////
void ras_show_version(void)
{
    cout << "Version is: \"" << GLOB_VER << "\"" << endl;
}


//////////////
void ras_show_help(void)
{
    cout<<"================================================================================"<<endl;
    //cout<<"Notation: The {--par} means that the parameter '--par' is optional."<<endl;
    //cout<<endl;
    cout<<"List of commands:"<<endl;
    cout<<"      '--help' or '-h'"<<endl;
    cout<<"      '--haps2ped'"<<endl;
    cout<<"      '--dataman'"<<endl;
    cout<<endl;
    
    cout<<"----------------------------------------"<<endl;
    cout<<"* To get help:"<<endl;
    cout<<"      'gap --help', 'gap -h' or simply 'gap'."<<endl;
    cout<<endl;
    
    cout<<"----------------------------------------"<<endl;
    cout<<"* To convert .haps to .ped:"<<endl;
    cout<<"      'gap --haps2ped --file [file_name] --out [file_out] --code01'. "<<endl;
    cout<<endl;
    cout<<"  --haps2ped:"<<endl;
    cout<<"        This is the Command."<<endl;
    cout<<"  --file:"<<endl;
    cout<<"        The input files are [file_name.haps] and [file_name.sample], so use"<<endl;
    cout<<"        [file_name] without extension."<<endl;
    cout<<"  --out (optional):"<<endl;
    cout<<"        The default output file name is [out.ped]."<<endl;
    cout<<"  --code01 (optional):"<<endl;
    cout<<"        The optional parameter '--code01' recodes the alleles to 0 or 1."<<endl;
    cout<<endl;
    
    cout<<"----------------------------------------"<<endl;
    cout<<"* Manage data:"<<endl;
    cout<<"      'gap --dataman --file [file_name] --selcol 1 --selrow 1 --out [file_out]'. "<<endl;
    cout<<"      'gap --dataman --file [file_name] --transpose --out [file_out]'. "<<endl;
    cout<<endl;
    cout<<"  --dataman:"<<endl;
    cout<<"        This is the Command."<<endl;
    cout<<"  --file:"<<endl;
    cout<<"        The input file is [file_name.ext] with extension."<<endl;
    cout<<"  --selcol (optional):"<<endl;
    cout<<"        Select odd (1) or even columns (2) (Default is 1)."<<endl;
    cout<<"  --selrow (optional):"<<endl;
    cout<<"        Select odd (1) or even rows (2) (Default is 1)."<<endl;
    cout<<"  --out (optional):"<<endl;
    cout<<"        The default output file name."<<endl;
    cout<<"  --transpose (optional):"<<endl;
    cout<<"        Transpose the input file and save it to the output."<<endl;
    cout<<endl;

    cout<<"================================================================================"<<endl;
}




int ras_read_map(string file_name)
{
    ifstream ifile(file_name.c_str());
    if(!ifile) throw("Error: can not open the file ["+file_name+"] to read.");

    long int n_r=ras_FileLineNumber(file_name);
    long int n_c=ras_FileColNumber(file_name);
    cout<< "Loading [" << file_name << "] with " << n_r << " rows and " << n_c << " cols ..." << endl;
    cout<< "  Number of SNPs: " << n_r << endl;

    matrix_map_chr_rs.resize(n_r,vector<string>(2, ""));
    matrix_map_cm.resize(n_r,0);
    matrix_map_bp.resize(n_r,0);

    
    vector<string> st1;

    string line;
    long int j=0,i=0;

    while (getline(ifile, line)){
        st1=split(line,sep);
        matrix_map_chr_rs[j][0]=st1[0];
        matrix_map_chr_rs[j][1]=st1[1];
        matrix_map_cm[j]=atof(st1[2].c_str());
        matrix_map_bp[j]=atol(st1[3].c_str());
        j++;
        cout << "\r  " << (j*100/n_r) << "% completed." << flush;
    }    
    cout << endl;
    ifile.clear();
    ifile.close();
    return 0;
}


int ras_read_ped(string file_name)
{
    ifstream ifile(file_name.c_str());
    if(!ifile) throw("Error: can not open the file ["+file_name+"] to read.");

    long int n_r=ras_FileLineNumber(file_name);
    long int n_c=ras_FileColNumber(file_name);
    cout<< "Loading [" << file_name << "] with " << n_r << " rows and " << n_c << " cols ..." << endl;
    n_IND=n_r;
    n_SNP=(n_c-6)/2;
    cout<< "  Number of inds: " << n_IND << endl;
    cout<< "  Number of SNPs: " << n_SNP << endl;

    
    matrix_ped.resize(n_IND, vector<string>(n_c, "") ); //FID IID PID MID SEX PHENO A1 A2 ...
    matrix_ped_Allele.resize(n_SNP, vector<string>(2, ""));
    matrix_ped_Allele_frq.resize(n_SNP, vector<long int>(2, 0));

    
    vector<string> st1;

    string line;
    long int j=0,i=0;

    while (getline(ifile, line)){
        st1=split(line,sep);
        for(i=0; i<st1.size(); i++)
        {
            matrix_ped[j][i]=st1[i];
        }
        j++;
        cout << "\r  " << (j*100/n_r) << "% completed." << flush;
    }    
    cout << endl;
    ifile.clear();
    ifile.close();
    
    //compute the allele freq.
    ras_comp_allele_freq();
    return 0;

}


int ras_comp_allele_freq(void)
{
    long int j=0,i=0,k=0;
    cout<< "Computing allele frequency ..." <<endl;

    // compute MAF
    for(i=0; i<n_SNP; i++){
        long int al1=0,al2=0;
        int co=0;
        for (j=0; j<n_IND; j++){
            for (k=0; k<=1; k++){
                string s1=matrix_ped[j][6+2*i+k];
                if (s1!="0"){
                    if (co==0){
                        matrix_ped_Allele[i][0]=s1;
                        al1++;
                        co++;
                    }
                    else if (co==1){
                        if (matrix_ped_Allele[i][0]==s1){
                            al1++;
                        }
                        else{                    
                            matrix_ped_Allele[i][1]=s1;
                            al2++;
                            co++;
                        }                
                    }
                    else{ //co==2
                        if (matrix_ped_Allele[i][0]==s1){
                            al1++;
                        }
                        else if (matrix_ped_Allele[i][1]==s1){
                            al2++;
                        }
                        else{
                            cout << "Error: more than two allele for SNP."<<endl;
                            throw("Error: more than two allele for SNP.");
                        }
                    }
                }
            }
        }
        matrix_ped_Allele_frq[i][0]=al1;
        matrix_ped_Allele_frq[i][1]=al2;
        cout << "\r  " << (i*100/(n_SNP-1)) << "% completed." << flush;
    }
    cout << endl;

    //test allele.frq
    //ofstream myfile("allel.frq");
    //for(i=0; i<n_SNP; i++){
    //    myfile << matrix_ped_Allele[i][0] << sep << matrix_ped_Allele[i][1] << endl;
    //}
    //myfile.close();
    
    return 0;
}



int ras_read_pheno(string file_name)
{
    ifstream ifile(file_name.c_str());
    if(!ifile) throw("Error: can not open the file ["+file_name+"] to read.");

    long int n_r=ras_FileLineNumber(file_name);
    long int n_c=ras_FileColNumber(file_name);
    cout<< "Loading [" << file_name << "] with " << n_r << " rows and " << n_c << " cols ..." << endl;
    cout<< "  Number of inds: " << n_r << endl;
    n_Pheno=n_c-2;
    cout<< "  Number of different phenotypes: " << n_Pheno << endl;

    matrix_pheno.resize(n_r , vector<double>(n_Pheno, -9) );

    vector<string> st1;

    string line;
    long int j=0,i=0;

    while (getline(ifile, line)){
        st1=split(line,sep);
        for(i=2; i<st1.size(); i++)
        {
            matrix_pheno[j][i-2]=atof(st1[i].c_str());
        }
        j++;
        cout << "\r  " << (j*100/n_r) << "% completed." << flush;
    }
    cout << endl;
    ifile.clear();
    ifile.close();
    
    return 0;
}





///////////////////////////
int ras_read_dat_int(string file_name)
{
    ifstream ifile(file_name.c_str());
    if(!ifile) throw("Error: can not open the file ["+file_name+"] to read.");
    
    long int n_r=ras_FileLineNumber(file_name); //nomber of rows
    long int n_c=ras_FileColNumber(file_name); //nomber of cols
    n_IND=n_r;
    n_SNP=n_c;
    
    
    cout<< "Allocating memory ..." << flush;
    matrix_dat.resize(n_r, vector<double>(n_c, 0) );
    cout<< ", Done." << endl;

    cout<< "Loading [" << file_name << "] with " << n_r << " rows and " << n_c << " cols ..." << endl;

    
    vector<string> st1;
    
    string line;
    long int j=0,i=0;
    
    while (getline(ifile, line)){
        st1=split(line,sep);
        for(i=0; i<st1.size(); i++)
        {
            matrix_dat[j][i]=atof(st1[i].c_str());
        }
        j++;
        cout << "\r  " << (j*100/n_r) << "% completed." << flush;
    }
    cout << endl;
    ifile.clear();
    ifile.close();
    
    return 0;
}


///////////////////////////
int ras_read_dat_int2(string file_name)
{
    ifstream ifile(file_name.c_str());
    if(!ifile) throw("Error: can not open the file ["+file_name+"] to read.");
    
    long int n_r=ras_FileLineNumber(file_name); //nomber of rows
    long int n_c=ras_FileColNumber(file_name); //nomber of cols
    n_IND2=n_r;
    n_SNP2=n_c;
    
    
    cout<< "Allocating memory ..." << flush;
    matrix_dat_testset.resize(n_r, vector<double>(n_c, 0) );
    cout<< ", Done." << endl;

    cout<< "Loading [" << file_name << "] with " << n_r << " rows and " << n_c << " cols ..." << endl;
    
    
    vector<string> st1;
    
    string line;
    long int j=0,i=0;
    
    while (getline(ifile, line)){
        st1=split(line,sep);
        for(i=0; i<st1.size(); i++)
        {
            matrix_dat_testset[j][i]=atof(st1[i].c_str());
        }
        j++;
        cout << "\r  " << (j*100/n_r) << "% completed." << flush;
    }
    cout << endl;
    ifile.clear();
    ifile.close();

    return 0;
}






//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

int RasDataMan::main_DataMan(void)
{
    if (flag_transpose){
        ras_read_dat_str(infile);
        write_transpose();
    }
    
    else
    {
        if (DataMan_OE_col>0 && DataMan_OE_row>0)
        {
            return 0;
        }
        else if (DataMan_OE_col>0)
        {
            ras_selcols();
            return 0;
        }
        else if (DataMan_OE_row>0)
        {
            return 0;
        }
        else
        {
            return 0;
        }
    }
    return 0;

}


int RasDataMan::ras_selcols(void)
{
    ifstream ifile(infile.c_str());
    if(!ifile) throw("Error: can not open the file ["+infile+"] to read.");
    
    ofstream ofile(outfile.c_str());
    
    long int n_r=ras_FileLineNumber(infile); //nomber of rows
    long int n_c=ras_FileColNumber(infile); //nomber of cols
    
    cout<< "Loading [" << infile << "] with " << n_r << " rows and " << n_c << " cols and saving to [" << outfile << "] ..." << endl;
    
    
    vector<string> st1;
    
    string line;
    long int j=0,i=0;
    
    long int st=0;
    if (DataMan_OE_col==1) st=1; //odd
    if (DataMan_OE_col==2) st=0; //even
    while (getline(ifile, line))
    {
        st1=split(line,sep);
        for(i=st; i<st1.size(); i=i+2)
        {
            ofile << st1[i] << sep;
        }
        ofile << endl;
        j++;
        cout << "\r  " << (j*100/n_r) << "% completed." << flush;
    }
    cout << endl;
    ifile.clear();
    ifile.close();
    ofile.close();
    return 0;

}


int RasDataMan::write_transpose(void)
{
    cout << "Writing the transposed to ["+outfile+"] ..." << endl;
    ofstream myfile(outfile.c_str());
    
    long int i,j;

    for (j=0; j<ncol; j++){ //first col, because it should be transposed
        for(i=0; i<nrow; i++){
            myfile << matrix_dat_str[i][j] << sep;
        }
        myfile << endl;
        cout << "\r  " << (j*100/(ncol-1)) << "% completed." << flush;
    }
    cout << endl;
    
    myfile.close();
    return 0;
}


int RasDataMan::ras_read_dat_str(string file_name)
{
    ifstream ifile(file_name.c_str());
    if(!ifile) throw("Error: can not open the file ["+file_name+"] to read.");
    
    long int n_r=ras_FileLineNumber(file_name); //nomber of rows
    long int n_c=ras_FileColNumber(file_name); //nomber of cols
    nrow=n_r;
    ncol=n_c;
    
    
    cout<< "Allocating memory ..." << flush;
    matrix_dat_str.resize(n_r, vector<string>(n_c, "") );
    cout<< ", Done." << endl;
    
    cout<< "Loading [" << file_name << "] with " << n_r << " rows and " << n_c << " cols ..." << endl;
    
    vector<string> st1;
    
    string line;
    long int j=0,i=0;
    
    while (getline(ifile, line)){
        st1=split(line,sep);
        for(i=0; i<st1.size(); i++)
        {
            matrix_dat_str[j][i]=st1[i];
        }
        j++;
        cout << "\r  " << (j*100/n_r) << "% completed." << flush;
    }
    cout << endl;
    ifile.clear();
    ifile.close();
    
    return 0;
}








