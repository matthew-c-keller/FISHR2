
#Example of preparing data for FISHR2, finding optimal parameters for FISHR2, and using FISHR2 on a small (n=8k, #Snps=1185, cM=16) simulated SNP dataset on chromosome 15 that is in PLINK binary format to start with

#Written by Matthew Keller
#Dec 7, 2016


#Why write this script in R? Most of the commands below are system() commands and can be run in unix terminal or a shell script. However, other commands (e.g., that look at the data using hist()) require R. So to keep things simple and all in a single script, we've chosen to write this up as an R script. We anticipate that most users using FISHR2 will run it in a shell script and look at their results in R.





##############################
# PROGRAMS REQUIRED BY THIS SCRIPT
#We assume that you download and compile (or use the precompiled binaries) the following software for this script to run. It goes without saying that you have R installed on your machine and that you're using this script in R!

#1) FISHR2 - https://github.com/matthew-c-keller/FISHR2

#2) GERMLINE2 - https://github.com/matthew-c-keller/GERMLINE - NOTE: this is a modified version of the original GERMLINE software

#3) parameter_finder - installed when you installed FISHR2

#4) gap - https://github.com/rtahmasbi/GAP -  - formats from SHAPEIT2 to phased PED file expected by GERMLINE2

#This script assumes that you have installed FISHR2 and your working directory is the one created when you clone FISHR2. We assume that you install GERMLINE2 and gap in that directory. If you install these programs in other locations, change the script below accordingly to indicate where those programs are (or put them in a directory that is in your $PATH and just call them directly).

##############################






##############################
#DATA INCLUDED IN THE FOLDER FOR THIS EXAMPLE SCRIPT
#Test.bed/bim/fam files. NOTE: the 3rd column of the bim file has cM positions, which the user will typically need to figure out for their particular SNP set themselves.
#hapmap.genmap15.txt - genetic map file
#For additional information on the genetic map files, see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gmap


#Test.SI.haps & Test.SI.sample - files created by SHAPEIT; to save you time (and not have to wait for SHAPEIT to finish phasing the data), we've supplied these files for you.
##############################






##############################
# SOME CRUCIAL ISSUES RE FORMATING DATA & USING PHASING DATA

# FISHR2 requires that SNP data is in one of the following two PHASED format types:
#1) PHASED PED file format. I.e., A/C C/T T/A G/G G/A means that ACTGG is one estimated phase and CTAGA is the other
#2) SHAPEIT file format (*.hap/*.sample)

# ***IMPORTANT IF YOU ARE USING PHASED PED FILES (#1 above) ***: If you use PLINK to manipulate your phased PED file in any way, your phase will be lost! So do all PLINK data cleaning first, then phase, then transform your data to phased ped file format, then use FISHR2. Never use PLINK on a phased ped file unless you are willing to lose phase information!

#Also note that the *map files MUST have a cM distance in the 3rd column for FISHR2/GERMLINE2 to work properly



# USEFUL UTILITIES & SOFTWARE

# SHAPEIT phase software: https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#home

# BEAGLE phase software: http://faculty.washington.edu/browning/beagle/beagle.html 

# Utilities for dealing with BEAGLE formated data: http://faculty.washington.edu/browning/beagle_utilities/utilities.html 
#beagle2vcf.jar
#vcf2beagle.jar

# Utility for dealing with GERMLINE & GERMLINE2 data & output: http://www.cs.columbia.edu/~gusev/germline/
#ped_to_bgl

##############################







##############
#1 Phase using SHAPEIT on PLINK BED format files

#Note: you can of course use whatever program you want to phase it, but it's up to you to convert the output of that program into SHAPEIT output or phased PED files. For a phasing pipeline from BEAGLE to phased PED file, see http://www.cs.columbia.edu/~gusev/germline/.

#Set working directory
setwd("/path/to/directory/FISHR2") # obviously, change /path/to/directory here to wherever you put FISHR2 & the utilitites

#We begin with Test.bed/bim/fam files as well as a genetic map file (hapmap.genmap15.txt)
#For additional information on the genetic map files, see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gmap


#NOTE: the following 2 commands are commented out. To save you time, the files created by these commands (Test.SI.haps & Test.SI.sample) are included in the directory

#system("shapeit.v2.r727.linux.x64 --input-bed Test.bed Test.bim Test.fam -M hapmap.genmap15.txt --thread 20 --output-max Test.SI.haps.gz Test.SI.sample") #phase using shapeit; took 50 minutes on my machine

#system("gunzip Test.SI.haps.gz") #gunzip the hap file
##############








###########NOTE#############
#Sections 2 through 4 are OPTIONAL, for users who want to user parameter_finder to choose optimal parameters for running FISHR2. For people who already know their parameters, these steps would be skipped.
############################



##############
#2 ***OPTIONAL*** Convert the SHAPEIT file format to phased PED files
#This is useful only if you're using gl.parameter finder & GERMLINE2 to find optimal parameters for running FISHR2
#Otherwise, it's fine to leave the files in the SHAPEIT format for input into FISHR2

#Use one of our utility programs to convert SHAPEIT output to what GERMLINE wants
#See also https://github.com/rtahmasbi/GAP for newest versions
system("./utilities/gap/gap --haps2ped --file Test.SI --code01 --out Test.SI")

#create a corresponding map file by changing the bim file. NOTE: This map file has cM distances in 3rd column 
system("cut -f1,2,3,4 Test.bim > Test.SI.map")
###############






###############
#3 ***OPTIONAL*** RUN GERMLINE2 FOR FINDING OPTIMAL THRESHOLDS USING PARAMETER_FINDER

#NOTE: this assumes that you have downloaded & compiled GERMLINE2, available at https://github.com/matthew-c-keller/GERMLINE2

#Note again: the PED file must be phased for this to work
#run GERMLINE to identify regions that are almost certainly IBD, at least in the middle
#If you get a library missing error, Navigate to the Germline2 folder and compile the program by issuing the make command.
system("./GERMLINE2/GERMLINE2 -pedfile Test.SI.ped -mapfile Test.SI.map -outfile Test.SI -bin_out -err_hom 0 -err_het 0 -reduced -bits 120 -min_m 8 -w_extend")  #ignore a warning that "stream is not good" if you get it

#-pedfile - a PHASED pedfile. 
#-mapfile - a map file with cM distances in 3rd column
#-bin_out - a compressed output necessary for being read by gl.parameter.finder & FISHR
#-err_hom 0 - allow 0 mismatching homozygous markers
#-err_het 0 - allow 0 mismatching heterozygous markers
#-reduced - a flag telling GERMLINE to reduce the columns of the output; necessary for being read by gl.parameter.finder & FISHR
#-bits 120 - the # SNPs in each fixed window that GERMLINE uses for initial matches
#-min_m - minimum length in cM for match to be output; here 8
#-w_extend extend match beyond the "bits" window until the first opposite homozygote (OH) occurs.
#-h_extend - tells GERMLINE to use phase information as well as just opposite homozygosity in determining the SH length. It is important NOT to use this option when running GERMLINE to find optimal parameters. This is because using phase information will bias the detected SHs to have fewer phase errors than randomly chosen truly IBD segments, and will make our sensitivity & specificity values based on the chosen thresholds below appear better than they really would be in real circumstances. Although it's true that OHs also cause IEs, they are a large minority of IEs, most of which are caused by phase & SNP errors.

system("less Test.SI.log")
#you should have identified 48367 IBD segments; almost all of these are truly IBD, even if often over-extended. For getting PIE and MA distributions of truly IBD segments, we'll use the middlemost 50% of these IBD segments to ensure that they're truly IBD.
###############




###############
#3) ***OPTIONAL*** GL.PARAMTER.FINDER TO FIND PIE & MA THRESHOLDS TO USE
#This assumes you have compiled the parameter_finder utility. You can also use a precompiled version in ./utilities/parameter_finder_binaries

system("./utilities/parameter_finder/parameter_finder -bmatch Test.SI.bmatch -bsid Test.SI.bsid -bmid Test.SI.bmid -ped-file Test.SI.ped -window 50 -cut-value 0.5 -reduced 500 8 -output-type Error1 -log-file Test.SI.paramfind | gzip > Test.SI.PF.gz")

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

###############





###############
#4) ***OPTIONAL*** LOOK AT GL.PARAMTER.FINDER OUTPUT AND DETERMINE PIE & MA THRESHOLDS FOR FISHR

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

#These seem like decent thresholds. Note that this method doesn't account for the increasing uncertainty that occurs at SH endpoints, where IEs begin to accumulate. In real data and for short IBD segments, the IEs ocurring at endpoints make up a greater and greater share of the total IBD segment length and increase the false negative and false positive rate over what's estimated here.

###############












###############
#6) ***ALTERNATIVE 1*** RUN FISHR2 ON PHASED PED FILE

system("./binaries/FISHR2 -pedfile Test.SI.ped -mapfile Test.SI.map -bits 60 -err_hom 0 -err_het 0 -min_cm_initial 1.5 -homoz -w_extend -h_extend -min_cm_final 3 -min_snp 64 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -emp-ma-threshold 0.045 -log-file Test.8k -germline_output ./ped.germ |gzip > Test.8k.PedFormat.FISHR2.gz")


#DESCRIPTION OF PARAMETER
#-pedfile - a PHASED pedfile. Once again, this is NOT a regular PED file as created by, e.g., PLINK!
#-mapfile - a map file with cM distances in 3rd column
#-bits 60 - the # SNPs in each fixed window that GERMLINE2 uses for initial matches
#-err_hom 0 - allow 0 mismatching homozygous markers
#-err_het 0 - allow 0 mismatching heterozygous markers
#-min_cm_initial - minimum length in cM for an initial candidate segment,  which FISHR2 will then process (potentially join with other segments, truncate, or drop)
# -homoz - output runs of homozygosity (within-person IBD)
#-w_extend - tells GERMLINE2 to extend match beyond the "bits" window until the first mismatch occurs
#-h_extend - tells GERMLINE2 to use phase information as well as just opposite homozygosity in determining the SH length. Here, we DO use phase information to increase the accuracy of our initial calls.
#-min_cm_final - the minimum length in cM of the final segments outputed by FISHR2
#-min_snp - the minimum length in SNPs of the final segments outputed by FISHR2
#-window 50 - the window used to calculate the moving average (MA) is 50 SNPs wide. The SNP in the middle of the window takes the average value of those 50 SNPs.
#-gap 5 - at the intial step, consolidate two SHs(each at least 1.5 cM long, as dictated by the min_cm_initial parameter) between the same two individuals if separated by fewer than 5 SNPs. If these are actually two separate IBD segments, they will often be broken up again at the MA step
#-output.type finalOutput - the only output most users will ever need. finalOutput [default] (soon to be "Full") - Outputs the final called SHs. The columns of this file are pers1, pers2, bp.distance.start, bp.distance.end, no.of.snps.in.match, and cm.distance
#-count-gap-errors TRUE - says to count IEs in the gap whenever two SHs are joined (via the gap argument above). If this is FALSE (not recommended), then once joined, the MA step will not break up SHs at the gap junctures.
#-emp-ma-threshold .045 - we break up SHs at the point where the MA value > .045. If the remaining SHs are < 3cM, the entire SH is dropped.
#-emp-pie-threshold .015 - at the final step after the MA trimming, remove any SH that the proportion of IEs > .015 across the whole called segment.
#-log.file - name of log file.
#-germline_output - *OPTIONAL* - name of the folder where GERMLINE2 (intermediate) output goes; if this parameter is not used, then no GERMLINE2 output will be created


#CHANGE
system("less Test.8k.log")
#you should have identified 276,556 IBD segments > 3cM in length and that passed the PIE threshold of .015 and MA threshold of .045 suggested above. GERMLINE2 originally detected 2,236,826 segments > 1.5 cM, and these were passed on to FISHR. 935,114 of these were dropped straight off (after any merging 80,842 separated by fewer than 5 SNPs), then 944,265 were dropped that were originally > 3cM but that trimming them based on the Moving Average of IEs made them < 3cM, and finally 49 (only) were removed because they had an overall proportion of IEs (PIE) > .015.

###############






###############
#7) ***ALTERNATIVE 2*** RUN FISHR2 ON PHASED SHAPEIT FORMATED DATA


system("./binaries/FISHR2 -mapfile Test.SI.map -samplefile Test.SI.sample -hapsfile Test.SI.haps -bits 60 -err_hom 0 -err_het 0 -min_cm_initial 1.5 -homoz -w_extend -h_extend -min_cm_final 3 -min_snp 64 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -emp-ma-threshold 0.045 -log-file SI.Test.8k |gzip > Test.8k.SIFormat.FISHR2.gz")

#PARAMETERS AS ABOVE EXCEPT:
#-mapfile - name of the mapfile; identical to the one used above when running FISHR2 on ped/map file. Note that SHAPEIT does not produce this file - you the user must create this (it's easily created if you phased your data using SHAPEIT from the *bim file)
#-samplefile - name of sample file produced by SHAPEIT
#-hapsfile - name of haps file produced by SHAPEIT
#NOTE: here we do not use the -germline_output parameter and therefore produce no GERMLINE2 intermediate output here

system("less SI.Test.8k.log")
#the output is identical to the above

###############





###############
#8) ***ALTERNATIVE 3*** RUN FISHR2 ON PHASED SHAPEIT FORMATED DATA AS ABOVE, BUT GET IBD2 AND IBD4 INFORMATION


system("./binaries/FISHR2 -mapfile Test.SI.map -samplefile Test.SI.sample -hapsfile Test.SI.haps -bits 60 -err_hom 0 -err_het 0 -min_cm_initial 1.5 -homoz -w_extend -h_extend -min_cm_final 3 -min_snp 64 -window 50 -gap 5 -output-type finalOutput -count.gap.errors TRUE -emp-pie-threshold 0.015 -emp-ma-threshold 0.045 -ibd2 2 -log-file SI.Test.8k.ibd2 |gzip > Test.8k.SIFormat.FISHR2.ibd2.gz")

#PARAMETERS AS ABOVE EXCEPT:
#-ibd2 2 - this is the minimum cM threshold for IBD2 & IBD4 segments; typically, we want this about the same length or a bit shorter than the -min_cm_final argument. The moving average for detecting IBD2 and IBD4 segments is set at the value of emp-ma-threshold

system("less SI.Test.8k.ibd2.log")
#the output is identical to the above

###############







###############
#9) Look at FISHR2 output

ff <- read.table(gzfile("Test.8k.SIFormat.FISHR2.gz"),header=FALSE,colClasses=c("character","character",rep("numeric",4)))
nrow(ff) #276556
names(ff) <- c('id1','id2','start','end','snps','cm')

#look at distribution of lengths
hist(ff$cm) #our genome is only 16 cM long, so some SHs stretch the entire genome.

#look at where they occur - plot the midpoints
ff$mdpt <- ff$start + (ff$end-ff$start)/2
hist(ff$mdpt,breaks=100,col='blue') #quite a bit of variability across the chromosome



# Look at IBD2/4 output
ffibd2 <- read.table(gzfile("Test.8k.SIFormat.FISHR2.ibd2.gz"),header=FALSE,colClasses=c("character","character",rep("numeric",4)))
nrow(ffibd2) #277804, so 277804-276556 = 1248 IBD2 or IBD4 segments
names(ffibd2) <- c('id1','id2','start','end','snps','cm','ibd')
summary(ffibd2)
summary(as.factor(ffibd2$ibd)) #276556 IBD1 segments, 1228 IBD2, and 20 IBD4


#what kind of numbers of IBD2 & IBD4 do we expect given IBD1 numbers?
npair <- (8000*7999)/2
prob.ibd <- 276556/npair
(prob.ibd^2)*npair #=2390=the number of IBD2 we expect, which is somewhat smaller than the number observed (1248), perhaps due to a higher false negative rate for detecting IBD2.
(prob.ibd^3)*npair #=20=the number of IBD4 we expect, which is exactly what is observed for IBD4
#This is to the 3rd power because 1) probability that 2 haplotypes between two people are IBD; 2) probability that the other 2 haplotypes are IBD given the first 2 are and; 3) probability that those two pair of haplotypes are the same (i.e., are themselves IBD). It would be a mistake to think this should be to the 4th power.
###############











































