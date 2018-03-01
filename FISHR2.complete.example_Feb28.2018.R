
#Example of preparing data for FISHR2, finding optimal parameters for FISHR2, and using FISHR2 on a small (n=8k, #Snps=1185, cM=16) simulated SNP dataset on chromosome 15 that is in PLINK binary format to start with. 

#Written by Matthew Keller
#Dec 7, 2016
#Updated Feb 12, 2018


#*************************************************************#
# THIS SCRIPT IS INTENDED TO BE RUN ON UNIX OPERATING SYSTEM. #
# FISHR2 HAS NOT BEEN TESTED OR COMPILED FOR NON-UNIX OS's.   #
# WE FURTHER ASSUME YOU ARE RUNNING THIS SCRIPT USING R.      #
#*************************************************************#




##############################
# PROGRAMS USED IN THIS SCRIPT

#REQUIRED PROGRAMS:
#1) FISHR2 - https://github.com/matthew-c-keller/FISHR2
#This script assumes that you are using the precomiled FISHR2 binary located in <path.where.you.installed.it>/FISHR2/binaries. If the precompiled version doesn't work for you or if you move FISHR2 to a new location, make sure you change the script (place the correct path to FISHR2)accordingly. We also assume that your working directory is the one where all the example data is located (FISHR2/src).


#2) R (obviously)
#Why write this script in R rather than in a shell script? Most of the commands below are system() commands and can be run in unix terminal or a shell script. However, other commands (e.g., that look at the data using hist()) require R. So to keep things simple and all in a single script, we've chosen to write this up as an R script. We anticipate that most users using FISHR2 will run it in a shell script and look at their results in R.



#OPTIONAL PROGRAMS:
#1) SHAPEIT - for phasing the data. We have already phased the data for you here, but you can phase it yourself in the optional part of the script. See http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html

#2) GERMLINE2 - https://github.com/matthew-c-keller/GERMLINE2 - NOTE: this is a modified version of the original GERMLINE software. This does NOT need to be installed for FISHR2 to work, but it does need to be installed if you wish to use parameter_finder utility. Install it from github using, e.g., git clone https://github.com/matthew-c-keller/GERMLINE2.git. Here, we assume that GERMLINE2 is installed in your /FISHR2/utilities folder. If you install it somewhere else, you'll need to modify the script below to provide the correct file path to where you've installed it.

#3) parameter_finder - installed when you installed FISHR2 and exists in the utilities folder, although it needs to be compiled using "make". It is used for finding optimal FISHR2 parameters.

#4) gap - https://github.com/rtahmasbi/GAP - installed when you installed FISHR2 and exists in the utilities folder, although it needs to be compiled using "make". Formats from SHAPEIT2 to phased PED file expected by GERMLINE2 and FISHR_Low_Ram

#Optional parts of this script assume that you have installed and compiled SHAPEIT, GERMLINE2, parameter_finder and gap in the utilities directory. If you install these programs in other locations, change the script below accordingly to indicate where those programs are (or put them in a directory that is in your $PATH and just call them directly).
##############################






##############################
#DATA USED IN THIS SCRIPT (in the src folder)

#REQUIRED FILES:
#Test.SI.haps & Test.SI.sample - files created by SHAPEIT from unphased Test.bed/bim/fam files (see below); to save you time (and not have to wait for SHAPEIT to finish phasing the data), we've supplied these phased haps & sample files for you. This script assumes you are using the SHAPEIT formatted files because they are much more common. For ease of use, we recommend you use FISHR2 in this way.

#Test.SI.map - the map file (for running FISHR2) MUST have a cM distance in the 3rd column for FISHR2 to work properly. Most map files you inherit will not include cM distances (they will be 0 for all markers). Thus, YOU the user needs to figure out what cM distance corresponds to each bp distance in your particular dataset/build. This script assumes that you have properly generated cM distances in the map file. If you want to be able to use FISHR2 without cM distances, you can place Mb distances in the 3rd column of the map file instead, but be aware that this is a quick and dirty approximation and IBD calling accuracy will suffer.



#OPTIONAL FILES FOR PHASING DATA:
#Test.bed/bim/fam  - these are unphased files that we then used as input to SHAPEIT, which phased this data to create Test.SI.haps and Test.SI.sample above. NOTE: the 3rd column of the Test.bim file has cM positions, which the user will typically need to figure out for their particular SNP set themselves. This column is not required for SHAPEIT (which only requires the genetic map file below), but we've included it here because we created Test.SI.map (used in FISHR2, which does require cM distances) from this file.

#hapmap.genmap15.txt - genetic map file. This is required for phasing using SHAPEIT. For additional information on the genetic map files, see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gmap



#OPTIONAL FILES FOR RUNNING GERMLINE2, FISHR_LOW_RAM, OR PARAMETER_FINDER
#Test.SI.ped & Test.SI.map - these are PHASED PED files - see http://www.cs.columbia.edu/~gusev/germline/ for file specifications. A/C C/T T/A G/G G/A means that ACTGG is one estimated phase and CTAGA is the other. You can convert a SHAPEIT phased *haps and *sample file to a phased ped file using our gap program, which is supplied with FISHR2. Note that if you use phased PED files (#2 above), you should not use PLINK to manipulate your phased PED file in any way. If you do, your phase will be lost! So do all PLINK data cleaning first, then phase, then transform your data to phased ped file format, then use FISHR2. Never use PLINK on a phased ped file unless you are willing to lose phase information! In general, be very careful aboutmanipulating phased data, and always ensure that your phase hasn't been destroyed.

##############################







##############
#1 Set directory, phase using SHAPEIT (commented out for speed) on PLINK BED format files, and create map file

#Note: you can of course use whatever program you want to phase it, but it's up to you to convert the output of that program into SHAPEIT output or phased PED files. For a phasing pipeline from BEAGLE to phased PED file, see http://www.cs.columbia.edu/~gusev/germline/.

#Set working directory
setwd("/path/to/directory/FISHR2/src") # obviously, change /path/to/directory here to wherever you put FISHR2 folder
#setwd("/work/KellerLab/mmkeller/FISHR/FISHR2.Feb2018/FISHR2/src") #for example, the path on my unix box

#OPTIONAL****
#We begin with Test.bed/bim/fam files as well as a genetic map file (hapmap.genmap15.txt)
#For additional information on the genetic map files, see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gmap. You will need to create your own genetic map file (a file 

#NOTE: the following 2 commands are commented out. To save you time, the files created by these commands (Test.SI.haps & Test.SI.sample) are included in the directory

#system("shapeit.v2.r727.linux.x64 --input-bed Test.bed Test.bim Test.fam -M hapmap.genmap15.txt --thread 20 --output-max Test.SI.haps.gz Test.SI.sample") #phase using shapeit; took 50 minutes on my machine

#system("gunzip Test.SI.haps.gz") #gunzip the hap file

#create a corresponding map file by changing the bim file. NOTE: This map file has cM distances in 3rd column 
system("cut -f1,2,3,4 Test.bim > Test.SI.map")

##############






###############
#2) RUN FISHR2 ON PHASED SHAPEIT FORMATED DATA

#Inputs below are to speed things up for demonstration purposes; see next for recommended inputs
system("../binaries/FISHR2 -mapfile Test.SI.map -samplefile Test.SI.sample -hapsfile Test.SI.haps -bits 120 -err_hom 0 -err_het 0 -min_cm_initial 4 -homoz -w_extend -h_extend -min_cm_final 5 -min_snp 100 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -emp-ma-threshold 0.045 -log-file SI.Test.8k |gzip > Test.8k.SIFormat.FISHR2.gz")

#Recommended inputs; this took ~ 5 minutes on my machine. This is commented out.
#system("../binaries/FISHR2 -mapfile Test.SI.map -samplefile Test.SI.sample -hapsfile Test.SI.haps -bits 60 -err_hom 0 -err_het 0 -min_cm_initial 1.5 -homoz -w_extend -h_extend -min_cm_final 3 -min_snp 64 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -emp-ma-threshold 0.045 -log-file SI.Test.8k |gzip > Test.8k.SIFormat.FISHR2.gz")


#DESCRIPTION OF INPUT PARAMETER
#-mapfile - name of the mapfile with cM distances in 3rd column. Note that SHAPEIT does not produce this file - you the user must create this (it's easily created if you phased your data using SHAPEIT from the *bim file)
#-samplefile - name of sample file produced by SHAPEIT
#-hapsfile - name of haps file produced by SHAPEIT
#-bits 60 - the # SNPs in each fixed window that GERMLINE2 uses for initial matches
#-err_hom 0 - allow 0 mismatching homozygous markers
#-err_het 0 - allow 0 mismatching heterozygous markers
#-min_cm_initial - minimum length in cM for an initial candidate segment,  which FISHR2 will then process (potentially join with other segments, truncate, or drop)
# -homoz - output runs of homozygosity (within-person IBD)
#-w_extend - tells GERMLINE2 (run internally) to extend match beyond the "bits" window until the first mismatch occurs
#-h_extend - tells GERMLINE2 (run internally) to use phase information as well as just opposite homozygosity in determining the SH length. Here, we DO use phase information to increase the accuracy of our initial calls.
#-min_cm_final - the minimum length in cM of the final segments outputed by FISHR2
#-min_snp - the minimum length in SNPs of the final segments outputed by FISHR2
#-window 50 - the window used to calculate the moving average (MA) is 50 SNPs wide. The SNP in the middle of the window takes the average value of those 50 SNPs.
#-gap 5 - at the intial step, consolidate two SHs(each at least 1.5 cM long, as dictated by the min_cm_initial parameter) between the same two individuals if separated by fewer than 5 SNPs. If these are actually two separate IBD segments, they will often be broken up again at the MA step
#-output.type finalOutput - the only output most users will ever need. finalOutput [default] (soon to be "Full") - Outputs the final called SHs. The columns of this file are pers1, pers2, bp.distance.start, bp.distance.end, no.of.snps.in.match, and cm.distance
#-count-gap-errors TRUE - says to count IEs in the gap whenever two SHs are joined (via the gap argument above). If this is FALSE (not recommended), then once joined, the MA step will not break up SHs at the gap junctures.
#-emp-ma-threshold .045 - we break up SHs at the point where the MA value > .045. If the remaining SHs are < 3cM, the entire SH is dropped.
#-emp-pie-threshold .015 - at the final step after the MA trimming, remove any SH that the proportion of IEs > .015 across the whole called segment.
#-log.file - name of log file.
#-germline_output - *OPTIONAL* - name of the folder where GERMLINE2 (intermediate) output goes; if this parameter is not used, then no GERMLINE2 output will be created. Note that here we do not use the -germline_output parameter and therefore produce no GERMLINE2 intermediate output here

#If you have an error about GLIBC (e.g., "GLIBC_2.14 not found...", you should make sure you have gcc installed on your system and in your $PATH. If that's the case and you still get the error, you probably need to recompile FISHR2. Just type "make" in the main folder (/path/to/directory/FISHR2/)


#Look at log file
system("less SI.Test.8k.log")
#you should have identified 75,046 IBD segments > 5cM in length and that passed the PIE threshold of .015 and MA threshold of .045 suggested above. GERMLINE2 originally detected 201,792 segments > 4 cM (the min_cm_initial parameter), and these were passed on to FISHR2. 60,385 of the total segments were dropped straight off because they were less than 5cM (after any merging 5382 separated by fewer than 5 SNPs), then 60,979 more were dropped that were originally > 5cM but that trimming them based on the Moving Average of IEs made them < 5cM, and finally 0 (in this case) were removed because they had an overall proportion of IEs (PIE) > .015. We are left with 75,046 segments that passed all QC. Note that the "Total time" output by the algorithm in the log file is that used by the FISHR2 part of the program and does not consider the initial and longer GERMLINE part of the algorithm (this is a bug we'll try to fix).

#Look at segments
segs <- read.table(gzfile("Test.8k.SIFormat.FISHR2.gz"),header=FALSE,stringsAsFactors=FALSE)
names(segs) <- c('p1','p2','start','end','snps','cm')
head(segs)
dim(segs) #75,046 segments detected
summary(segs) #we can see that the -min_snp argument of 100 had no effect; the smallest segment is 242 SNPs in length
hist(segs$cm) #our genome is only 16 cM long, so some SHs stretch the entire genome.
plot(segs$snps,segs$cm)
sa <- segs[segs$p1==segs$p2,] #note there are 4 runs of homozygosity detected

#look at where the segments occur - plot the midpoints
segs$mdpt <- segs$start + (segs$end-segs$start)/2
hist(segs$mdpt,breaks=100,col='blue') #quite a bit of variability across the chromosome
###############







###############
#3) RUN FISHR2 ON PHASED SHAPEIT FORMATED DATA TO GET IBD2 INFORMATION

#Inputs below are to speed things up for demonstration purposes; see next for recommended inputs
system("../binaries/FISHR2 -mapfile Test.SI.map -samplefile Test.SI.sample -hapsfile Test.SI.haps -bits 120 -err_hom 0 -err_het 0 -min_cm_initial 4 -w_extend -h_extend -min_cm_final 5 -min_snp 100 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -emp-ma-threshold 0.045 -ibd2 4 -log-file SI.Test.8k.IBD2 |gzip > Test.8k.SIFormat.FISHR2.IBD2.gz")

#Recommended inputs; this took ~ 5 minutes on my machine. This is commented out.
#system("../binaries/FISHR2 -mapfile Test.SI.map -samplefile Test.SI.sample -hapsfile Test.SI.haps -bits 60 -err_hom 0 -err_het 0 -min_cm_initial 1.5 -homoz -w_extend -h_extend -min_cm_final 3 -min_snp 64 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -emp-ma-threshold 0.045 -ibd 2.5 -log-file SI.Test.8k |gzip > Test.8k.SIFormat.FISHR2.gz")

#Note the use of -ibd2 4 argument, which tells FISHR2 to output IBD2 (or ibd4) segments that are >= 4cM
#IMPORTANT - when using ibd2, do not use the -homoz argument; detecting runs of homozygosity will interact poorly with the -ibd2 flag. When you do not use -homoz, no IBD4 segments will be detected.

#Look at segments
segs2 <- read.table(gzfile("Test.8k.SIFormat.FISHR2.IBD2.gz"),header=FALSE,stringsAsFactors=FALSE)
names(segs2) <- c('p1','p2','start','end','snps','cm','ibd')
head(segs2)
dim(segs2) #75,179 segments detected
summary(as.factor(segs2$ibd)) #we can see that there are 137 IBD2 segments and 0 IBD4 segments.
segs2$pairs <- paste(segs2$p1,".",segs2$p2,sep='')
ibd2.pairs <- segs2$pairs[segs2$ibd=="IBD2" | segs2$ibd=="IBD4"]
segs2.ibd2 <- segs2[segs2$pairs %in% ibd2.pairs,]
segs2.ibd2 #note that the IBD2 segments are always as long or shorter than the IBD1 segments. This is a expected, based on how IBD2 detection works.

#what kind of numbers of IBD2 & IBD4 do we expect given IBD1 numbers?
npair <- (8000*7999)/2
prob.ibd <- 75046/npair
(prob.ibd^2)*npair #=176=the number of IBD2 we expect, which is close to the number observed (137)
###############








##############
#4) ***OPTIONAL*** Convert the SHAPEIT file format to phased PED files and run FISHR2 on phased PED file
#Using Phased PED files is useful only in two situations:
#1) if you're using gl.parameter finder & GERMLINE2 to find optimal parameters for running FISHR2
#2) if you're running FISHR_Low_Ram (which only inputs phased PED files) because FISHR2 RAM's usage is too high 
#3) (typically not useful) You want to run FISHR2 using Phased PED format. That's OK I guess, but one would usually just run them on SHAPEIT formatted data

#Use one of our utility programs to convert SHAPEIT output to what GERMLINE2 and FISHR_Low_Ram wants. The following two commands will over-write the existing Test.SI.ped and Test.SI.map files
#See also https://github.com/rtahmasbi/GAP for newest versions
#This assumes you have compiled gap program in the ../utilities/gap folder!
system("../utilities/gap/gap --haps2ped --file Test.SI --code01 --out Test.SI") #convert Test.SI.sample & Test.SI.haps to phased PED. --code01 tells gap that the SNPs are in 0/1 format.


#Run FISHR2 on the phased PED data
system("../binaries/FISHR2 -pedfile Test.SI.ped -mapfile Test.SI.map -bits 120 -err_hom 0 -err_het 0 -min_cm_initial 4 -w_extend -homoz -h_extend -min_cm_final 5 -min_snp 100 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -emp-ma-threshold 0.045 -log-file SI.Test.8k.PedFormat |gzip > Test.8k.PedFormat.FISHR2.gz")


#DESCRIPTION OF PARAMETER
#-pedfile - a PHASED pedfile. Once again, this is NOT a regular PED file as created by, e.g., PLINK!
#-mapfile - a map file with cM distances in 3rd column
#The rest of parameters are identical to those in #2 above


#Look at segments
segs.ped <- read.table(gzfile("Test.8k.PedFormat.FISHR2.gz"),header=FALSE,stringsAsFactors=FALSE)
names(segs.ped) <- c('p1','p2','start','end','snps','cm')
head(segs.ped)
dim(segs.ped) #75,046 segments detected, exactly as with using the SHAPEIT input
summary(segs.ped) 

##############







##############
#5) ***OPTIONAL*** Run FISHR_Low_Ram on Phased PED data

#Run GERMLINE2 as a precursor to FISHR_Low_Ram. When using the low RAM version of FISHR, you must first run GERMLINE2, and then run FISHR_Low_Ram on that data created by GERMLINE2
#NOTE: this assumes that you have downloaded GERMLINE2 from github (https://github.com/matthew-c-keller/GERMLINE2) and have installed it in the same directory FISHR2 is in (i.e., 2 directories up from the src/ working directory we're in right now).

system("../../GERMLINE2/GERMLINE2 -pedfile Test.SI.ped -mapfile Test.SI.map -outfile Test2.SI.GL -bin_out -err_hom 0 -err_het 0 -reduced -bits 120 -min_m 4 -w_extend -h_extend")  #ignore a warning that "stream is not good" if you get it

#Parameters:
#-pedfile - a PHASED pedfile. 
#-mapfile - a map file with cM distances in 3rd column
#-bin_out - a compressed output necessary for being read by gl.parameter.finder & FISHR
#-err_hom 0 - allow 0 mismatching homozygous markers
#-err_het 0 - allow 0 mismatching heterozygous markers
#-reduced - a flag telling GERMLINE to reduce the columns of the output; necessary for being read by gl.parameter.finder & FISHR
#-bits 120 - the # SNPs in each fixed window that GERMLINE uses for initial matches
#-min_m - minimum length in cM for match to be output; here 4
#-w_extend extend match beyond the "bits" window until the first opposite homozygote (OH) occurs.
#-h_extend - tells GERMLINE to use phase information as well as just opposite homozygosity in determining the SH length.

#If you want to look at the GERMLINE2 segments before FISHR QC's them, you can use parse_bmatch to do that; it may need to be compiled
#system("../../GERMLINE2/parse_bmatch Test2.SI.GL.bmatch Test2.SI.GL.bsid Test2.SI.GL.bmid | gzip > Test.SI.GL.match.gz")
#gl.segs <- read.table(gzfile("Test.SI.GL.match.gz"),header=FALSE,stringsAsFactors=FALSE)



#Now run FISHR_Low_Ram. Note that this must be downloaded from here (https://github.com/matthew-c-keller/FISHR) and compiled in the folder /ErrorFinder23.3_Low_Ram. This requires MUCH less RAM than usual FISHR
system("../../FISHR/ErrorFinder23.3_Low_Ram/FISHR_Low_Ram -ped-file Test.SI.ped -bmatch Test2.SI.GL.bmatch -bsid Test2.SI.GL.bsid -bmid Test2.SI.GL.bmid -reduced 100 5 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold .015 -empirical-ma-threshold .045 -log-file Test.SI.FLR | gzip > Test.SI.FLR.gz")

#Look at segments
segs.flr <- read.table(gzfile("Test.SI.FLR.gz"),header=FALSE,stringsAsFactors=FALSE)
names(segs.flr) <- c('p1','p2','start','end','snps','cm')
head(segs.flr)
dim(segs.flr) #76382 segments (compared to 75,046 from FISHR2; the algorithm has changed slightly since FISHR_Low_Ram)
summary(segs.flr)

###############








###############
#6) ***OPTIONAL*** RUN GERMLINE2 AND PARAMTER_FINDER TO FIND PIE & MA THRESHOLDS TO USE
#This section is for users who want to user parameter_finder to choose optimal -emp-ma-threshold and -emp-pie-threshold parameters. For people who already know their parameters or want to use the ones we recommend, these steps would be skipped.

#Note again: the PED file must be phased for this to work
#run GERMLINE to identify regions that are almost certainly IBD, at least in the middle
#If you get a library missing error, Navigate to the Germline2 folder and compile the program by issuing the make command.
system("../../GERMLINE2/GERMLINE2 -pedfile Test.SI.ped -mapfile Test.SI.map -outfile Test.SI -bin_out -err_hom 0 -err_het 0 -reduced -bits 120 -min_m 8 -w_extend")  #ignore a warning that "stream is not good" if you get it

#-h_extend - tells GERMLINE to use phase information as well as just opposite homozygosity in determining the SH length. It is important NOT to use this option when running GERMLINE to find optimal parameters. This is because using phase information will bias the detected SHs to have fewer phase errors than randomly chosen truly IBD segments, and will make our sensitivity & specificity values based on the chosen thresholds below appear better than they really would be in real circumstances. Although it's true that OHs also cause IEs, they are a large minority of IEs, most of which are caused by phase & SNP errors.

system("less Test.SI.log")
#you should have identified 48367 IBD segments; almost all of these are truly IBD, even if often over-extended. For getting PIE and MA distributions of truly IBD segments, we'll use the middlemost 50% of these IBD segments to ensure that they're truly IBD.




#Use parameter_finder. This assumes you have compiled the parameter_finder utility. You can also use a precompiled version in ./utilities/parameter_finder_binaries
system("../utilities/parameter_finder/parameter_finder -bmatch Test.SI.bmatch -bsid Test.SI.bsid -bmid Test.SI.bmid -ped-file Test.SI.ped -window 50 -cut-value 0.5 -reduced 500 8 -output-type Error1 -log-file Test.SI.paramfind | gzip > Test.SI.PF.gz")

#Inputs to gl_parameter_finder:
#-bmatch - the binary output (from -bin_out) of GERMLINE; has one SH per row
#-bsid - the subject IDs output from GERMLINE
#-bmid - the marker IDs output from GERMLINE
#-ped-file - the path to the phased data; should be identical to the input -pedfile used in GERMLINE
#-window 50 - says to use a 50 SNP window for figuring moving average of IEs (MA)
#-cut-value .5 - says to trim the SH to the middlemost 50%: 25% from the left and 25% from the right are trimmed. This ensure that the remaining segment is IBD
#-reduced 500 8 - recommended paramers in real data; only use initial SHs that are at least 500 SNPs long and > 8cM. The outputed middlemost segments will typically be < 8cM (e.g., ~ 4+ cM)
#-output.type Error1 - this is the typical output format to be requested.
# | gzip > this pipes the standard out to gzip so the final output is a gzipped file.



#Look at parameter_finder output to determine PIE & MA thresholds for use with FISHR2
#We will read in different data than generated above - this is data phased with 8k individuals and is much higher quality
x <- read.table(gzfile("Test.SI.PF.gz"),header=TRUE) #Note: you can use fread() function from data.table library for MUCH faster reading
nrow(x)
#Note, the file has slightly fewer (48210) segments than the 48367 segments outputed by GERMLINE2. This is because a small proportion of the SHs from GERMLINE2 are not > 500 SNPs in length, even though all are at least 8cM long.

#The columns of the outputed gzip file are as follows:
#id1 - ID of the first person 
#id2 - ID of the second person who shares a SH with the first person
#marker_id1 - bp position for the first SNP of the middlemost section as defined by -cut-value (e.g., the middle 50% of the original SH if -cut-value is .50)
#marker_id2 - bp position for the last SNP of the middlemost section as defined by -cut-value
#snp_length - number of SNPs in the middlemost section
#cm_length - cM length of the middlemost section
#start.end - two values. first value is how many SNPs in the middlemost section started. So a value of 100 would mean that the original SH started at 1 (always) and the middlemost section started 100 SNPs into this. The second value is the final SNP of the middlemost section.
#pie - the proportion of IEs in this middlemost section
#ma_max - the maximum moving average (MA) of IEs in the middlemost section. The MA is averaged across the number of SNPs defined in -window.
#random_pie - the PIE for two random people at the same location (and length) as the middlemost section.
#random_ma_max - the MA for two random people at the same location (and length) as the middlemost section
#errors - where the IEs (separated by /) occurred in the middlemost section for the two original people in id1 and id2
#moving_averages - a vector (separated by /) of the MA values across all the SNPs in the middlemost section for the two original people in id1 and id2

x2 <- x[,1:11] #remove the moving average column

#histograms to visualize thresholds of .015 -emp-pie-threshold and .045 -emperical-ma-threshold
op <- par(mfrow=c(2,2))
hist(x2$pie,xlim=c(0,.12),breaks=50,col='lightblue',main='histogram of PIE in truly IBD segments')
abline(v=.015,col='red')
hist(x2$ma_max,xlim=c(0,.3),breaks=50,col='lightgreen',main='histogram of max MA in truly IBD segments')
abline(v=.045,col='red')
hist(x2$random_pie,xlim=c(0,.12),breaks=50,col='darkblue',main='histogram of PIE in non-IBD segments')
abline(v=.015,col='red')
hist(x2$random_ma_max,xlim=c(0,.3),breaks=50,col='darkgreen',main='histogram of max MA in non-IBD segments')
abline(v=.045,col='red')
par(op)


#how many truly IBD segments would be dropped given thresholds of .015 -emp-pie-threshold and .045 -emperical-ma-threshold?
x2$would.drop1 <- x2$ma_max > .045 | x2$pie > .015
summary(x2$would.drop1) #2022 of 48210
sum(x2$would.drop1)/nrow(x2) #2022/48210 = .042 is estimate of 1-sensitivity for this threshold for long SHs, so sensitivity ~ .956

#how many non-IBD segments would be retained given cutoffs above?
x2$would.keep1 <- x2$random_ma_max < .045 & x2$random_pie < .015
summary(x2$would.keep1) #72/48210, estimate of 1-specificity
sum(x2$would.keep1)/nrow(x2) #72/48210 = .00149 is estimate of 1-specificity for this threshold for long SHs. So we would make very few wrong calls proportionately. Nevertheless, because the base rate of non-IBD is so much (~200 times) higher than IBD, the PPV can still be <<1, even with such high sensitivity and specificity values. So with these values, we might predict a PPV of ~.80 with a 200 folder higher rate of non-IBD vs. IBD.

#These seem like decent thresholds. Note that this estimate doesn't account for the increasing uncertainty that occurs at SH endpoints, where IEs begin to accumulate. In real data and for short IBD segments, the IEs ocurring at endpoints make up a greater and greater share of the total IBD segment length and increase the false negative and false positive rate over what's estimated here.

###############


































