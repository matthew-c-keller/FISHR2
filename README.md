# FISHR2
FISHR2 is the newest version of the FISHR (Find IBD Shared Haplotypes Rapidly) program and supplants the old version, FISHR. FISHR2 detects IBD segments between individuals from genomewide SNP data. It is implemented in C++ and designed to work on a UNIX/LINUX environment. Major changes from V1:

1) Can be used directly on SHAPEIT formatted (HAPS/SAMPLE) data in addition to the "phased PED" file format used by GERMLINE and FISHR v.1. 

2) Running FISHR2 is now a single-step process. The user runs FISHR2 directly on the phased data rather than running GERMLINE2 first on the phased data and then FISHR on the GERMLINE output.

3) IBD2 and IBD4  shared segments (where 2 or 4 IBD segments exist at the same location between two individuals) can now be detected by FISHR2.



# Using FISHR2
We highly recommend looking at a complete worked through example, which we supply in the following R script: FISHR2.complete.example.R

The FISHR2 manual is available on the main repository (FISHR2.manual.pdf).

The original publication on FISHR v.1 describes the basic algorithm FISHR2 also uses and compares its performance to several leading IBD detection algorithms. This manuscript (doi:10.1038/ejhg.2017.6) is currently in press and is or will be available at http://matthewckeller.com/html/publications.html once published. It can also be obtained as an advanced online publication at http://www.nature.com/ejhg/journal/vaop/ncurrent/pdf/ejhg20176a.pdf

Use "./FISHR2 -help" command to view the list of flags and their usage. 



# To Download:
- 1. In your terminal write  
"git clone https://github.com/matthew-c-keller/FISHR2.git" to download the repository

- 2. Other way is to go to this link
"https://github.com/matthew-c-keller/FISHR2"
Click clone or download -> Download zip -> Extract zip file.

- 3. If you want to install git use:
	- For ubuntu based system:
	sudo apt-get update
	sudo apt-get install git
	- For Mac users, visit the link below and follow the instructions:
	http://mac.github.com.
	- For other operating systems, visit the following link and follow instructions:
	https://help.github.com/articles/set-up-git/


# To Compile FISHR2: 
- Navigate to folder where the makefile exists and type 'make' in your command line. 
- This will create a binary file FISHR2 in the directory.
- To execute the program refer to the examples below or the R example file above.


Known Issues:
See Instructions.txt
- If you don't have git installed, see the "To Download" section above for instructions.
- If the boost libraries are missing on your system, FISHR2 will not properly compile. 
FIX: sudo apt-get install libboost-all-dev.
This will install boost in your Unix/mac system.
- If there is a problem installing boost, there is a binary file called 'FISHR2' already included. You needn't compile it. You can just start using it. The examples below use this precompiled version of FISHR2.



# Examples:

./binaries/FISHR2 -mapfile ./src/Beagle.Phased.Group2.1k.map -pedfile ./src/Beagle.Phased.Group2.1k.ped  -bits 200 -err_hom 0 -err_het 0  -min_cm_initial 3 -homoz  -w_extend -h_extend -min_cm_final 3 -min_snp 64 -window 50 -gap 100 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -ma-threshold 0.2 -log-file logs |gzip > FISHR_OUT.gz


./binaries/FISHR2 -geneticmapfile ./src/big.gens -samplefile ./src/big.sample -hapsfile ./src/big.hap   -bits 200 -err_hom 0 -err_het 0  -min_cm_initial 1.5 -homoz  -w_extend -h_extend -min_cm_final 3 -min_snp 64 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -ma-threshold 0.2 -log-file logs |gzip > FISHR_OUT.gz


./binaries/FISHR2  -mapfile ./src/Beagle.Phased.Group2.1k.map -pedfile ./src/Beagle.Phased.Group2.1k.ped  -bits 200 -err_hom 0 -err_het 0  -min_cm_initial 1.5 -homoz  -w_extend -h_extend -min_cm_final 3 -min_snp 64 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -ma-threshold 0.2 -log-file logs -germline_output ./newFolderHere/ |gzip > FISHR_OUT.gz


./binaries/FISHR2 -mapfile ./src/Beagle.Phased.Group2.1k.map -pedfile ./src/Beagle.Phased.Group2.1k.ped  -bits 200 -err_hom 0 -err_het 0  -min_cm_initial 1.5 -homoz  -w_extend -h_extend -min_cm_final 3 -min_snp 64 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -ma-threshold 0.2 -ibd2 2.0   -log-file  logs -germline_output ./src/TestFolder |gzip > FISHR_OUT.gz


./binaries/FISHR2 -mapfile ./src/Beagle.Phased.Group2.1k.map -pedfile ./src/Beagle.Phased.Group2.1k.ped  -bits 200-err_hom 0 -err_het 0  -min_cm_initial 3 -homoz  -w_extend -h_extend -min_cm_final 3 -min_snp 64 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -ma-threshold 0.2 -ibd2 2.0   -log-file  logs -germline_output ./src/TestFolder |gzip > FISHR_OUT.gz

./binaries/FISHR2 -mapfile ./src/big.map -hapsfile ./src/big.hap -samplefile ./src/big.sample  -bits 200 -err_hom 0 -err_het 0  -min_cm_initial 3 -homoz  -w_extend -h_extend -min_cm_final 3 -min_snp 64 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -ma-threshold 0.2 -ibd2 2.0   -log-file  logs -germline_output ./src/TestFolder |gzip > FISHR_OUT.gz
