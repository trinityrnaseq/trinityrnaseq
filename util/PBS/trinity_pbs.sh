#!/bin/bash 
set -e  # turn on exit on error
##################################################################################################################################
##########################                                                                ########################################
##########################     Trinity PBS job submission with multi part dependencies    ########################################
##########################                                                                ########################################
##################################################################################################################################
### Author: Josh Bowden, Alexie Papanicolaou, CSIRO
### Version 1.0
###
###		
###		Script to split the Trinity workflow into multiple stages so as to efficiently request 
###		and use appropriate resources (walltime and number of cores) on a computer cluster / supercomputer.
###		Currently creates scripts for PBS Torque or PBSpro
###	    
###		trinity_pbs script install instructions:
###			1. Copy all trinity_pbs.* files into a directory (we will call it "TRINITY_PBS_DIR").
###			2. Add TRINITY_PBS_DIR to the PATH i.e. export or set PATH=TRINITY_PBS_DIR:$PATH  (perhaps export PATH in .bashrc file)
###			3. Change the "TRINITYPBSPATH" variable found below to point to the directory also. i.e. TRINITYPBSPATH=TRINITY_PBS_DIR
###			4. Set MEMDIRIN to name of a node-local filesystem so a network drive is not needed unecesarily for Scripts 4b and 5b
###			5. Set MODTRINITY to any modules that need to be loaded so Trinity.pl can be run.
###			6. Set TRINITYPATH to the path to Trinity.pl executable
###			7. Set PBSTYPE to --pbspro or --pbs, dependent on the system present.
###		That should be all that is needed from an admin perspective (besides making scripts accessible and exectable for users)			
###
###		Users need make a copy of TRINITY.CONFIG.template and then modify variables in it. See TRINITY.CONFIG.template for further details.
###		
###		The current script does the following.
###		Part 1. Reads data from TRINITY.CONFIG and creates the input directory, data file names, output data directory and Trinity.pl command line 
###			User inputs from TRINITY.CONFIG file :
###			JOBPREFIX				A string of less than 11 characters long. PBS will use this as a jobname prefix.
###			DATADIRECTORY  			Where input data exists
###			OUTPUTDIR  				Where user wants output data to go - requires a lot of space even for small datatsets 
###			STANDARD_JOB_DETAILS 	the Trinity.pl command line
###			ACCOUNT  				Account details of user (if required by PBS system being used)
###	   		
###		Part 2. Writes scripts to run Trinity.pl in 6 stages: 
###			3 intial (Inchworm, and 2 x Chrysalis stages: Chrysalis::GraphFromFasta and Chrysalis::ReadsFromTranscripts) 
###			2 parallel stages (Chrysalis::QuantifyGraph and Butterfly) which are executed in parallel.
###			1 collection of results as Trinity.Fasta.
###	      			
###			Information input from command line filename for stage 'x' : 
###			WALLTIME_Px  			Amount of time stage requires
###			MEM_Px					The amount of memory the stage requires
###			NCPU_Px     			The number of CPUs the stage may use
###			PBSNODETYPE _Px				The PBS (for --pbspro only) queue name
###			NUMPERARRAYITEM_Px		The number of massively parallel jobs in each parallel satge.
###	          		                                     
###		Part 3. Runs scripts dependant upon what stage has been detected as completed, using PBS job dependencies  
###                                                                                                      
###			Command line usage:
###				To start (or re-start) an analysis:
###					>trinity_pbs.sh TRINITY.CONFIG.template
###				To stop previously started PBS jobs on the queue:
###					>trinity_kill.pl OUTPUTDIR 
###				Where:
###					TRINITY.CONFIG.template = user specific job details
###					OUTPUTDIR   = is path to output data directory
### 
### Output job script submission files. These are saved in the output directory (OUTPUTDIR) and can be modified/re-run if any job fails.
### *_run.sh       Runs all the following scripts - with job dependencies and only the jobs that still need to be run.
### *_p1.sh        Runs Inchworm stage. Does not scale well past a single socket. Only request at most the number of cores on a single CPU.
### *_p2.sh        Runs Chrysalis::GrapghFromFasta clustering of Inchworm output. Should scale to number of cores on node
### *_p3.sh        Runs Chrysalis::ReadsToTranscripts. I/O limited. Try to use local filesystem (not implemented)
### *_p4a.sh       Creates jobs to run Chrysalis::QuantifyGraph and Butterfly parallel tasks
### *_p4b.sh       QuantifyGraph job. "NUMPERARRAYITEM" tasks from the file /chrysalis/quantifyGraph_commands are run for each job 
### *_p5b.sh       Butterfly job. "NUMPERARRAYITEM" tasks from the file /chrysalis/butterfly_commands are run for each job 
###                to start off at last completed stage. At present leaves all data on temporary area of shared network drive and 
###                copies Trinity.fatsa to home directory (with specific job prefix in filename). 
### * = $JOBPREFIX. $JOBPREFIX should not be > 10 characters long


##################################################################################################################################
####################### SET TRINITY INSTALLATION PATH (where Trinity.pl resides #################################################
##################################################################################################################################
# We need the path even if loaded using module
TRINITYPATH="/home/pap056/software/trinity_2013_08_14/"
# If you are loading using module, you can make the next variable blank
NEWPATH=
NEWPATH="export PATH=$PATH:$TRINITYPATH"

##################################################################################################################################
########    Set TRINITYPBSPATH to the directory where trinity_pbs.sh scripts are installed
##################################################################################################################################
TRINITYPBSPATH=`dirname "$0"`; # set to location of this script

##################################################################################################################################
########  Set cluster specific name for compute node local filesystem
##################################################################################################################################
MEMDIRIN="\$TMPDIR"  # available on Barrine  


##################################################################################################################################
########     Set system specific PATHS and load system specific modules (if available) 
##################################################################################################################################
# Example for for Barrine:
PBSTYPE="--pbspro"
# Here we ensure that Java 1.6 is used and Java 1.7 is removed (Butterfly dependency)
MODTRINITY="
 module load mpt/2.00 perl/5.15.8 bowtie/12.7 jellyfish/1.1.5 samtools/1.18 java/1.6.0_22-sun;
 module rm java/1.7.0_02
"

## That should be all the admin modifications needed.

# Append Trinity path data to env variables loaded by every script
MODTRINITY="
 $NEWPATH;
 export TRINITYPATH="/home/pap056/software/trinity_2013_08_14";
 $MODTRINITY
"
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

##################################################################################################################################
#########    Load files that contains functions 
##################################################################################################################################
# Modify function F_GETNODESTRING in file trinity_pbs.header so that a correct PBS header is returned to suit your PBS cluster
if [ -e "$TRINITYPBSPATH"/trinity_pbs.header ] ; then
  	source "$TRINITYPBSPATH"/trinity_pbs.header 
else
	echo "$1 requires file \"trinity_pbs.header\" to be present in: "
	echo "$TRINITYPBSPATH"
	exit 1
fi

##################################################################################################################################
#########    Load input config file
##################################################################################################################################
if [ -e "$1" ] ; then
  	source "$1"
else
	echo "Error: Input file does not exist: "$1"  "
	exit 1
fi


##################################################################################################################################
## Common variables to PBS and PBSpro
## and other needed variables that a user should not need to modify
##################################################################################################################################
if [ $UEMAIL ]; then PBSUSER="#PBS -M "$UEMAIL"" ; fi
HASHBANG="#!/bin/bash"

##################################################################################################################################
## PBS torque and PBSpro have some differences. 
## Organise these here and also check further on (line 184) and change NODETYPE to match cluster system
## MODTRINITY will also be different on different clusters - it sets up the paths to the required executables
##################################################################################################################################
if [[ "$PBSTYPE" = "--pbspro" ]] ; then	
	JOBARRAY="-J"
	JOBARRAY_ID="\$PBS_ARRAY_INDEX"
	AFTEROKARRAY="afterok"	
elif [[ "$PBSTYPE" = "--pbs" ]] ; then
	## PBS torque:
	JOBARRAY="-t"
	JOBARRAY_ID="\$PBS_ARRAYID"
	AFTEROKARRAY="afterokarray"	
else # no paramaters present
	F_USAGE	
	exit 0
fi

##############################################################################################################################################
##########     Part 1: Set up file names for input directory and for output data dir and Trinity.pl command line          ####################
##############################################################################################################################################
	echo ""
	## Ensure JOBPREFIX is not greatr than 11 characters as PBS-pro can not handle > 15 characters for total job name length
	echo "submitting trinity jobs with prefix: "
	echo "		$JOBPREFIX"	

	###### Set input data directory - $DATADIR is CSIRO specific	
	echo "Input directory: "
	echo "		$DATADIRECTORY"

	###### Set output data directory (OUTPUTDIR)
	echo "Output directory: Scripts and output data will be written to:"
	echo "		$OUTPUTDIR" 
	mkdir -p "$OUTPUTDIR"
	cd "$OUTPUTDIR"

	
	### Modify STANDARD_JOB_DETAILS for analysis specific input to Trinity.pl 
	echo "The following trinity command line will be run:"
	echo "$STANDARD_JOB_DETAILS" 
	echo "" 
	if [[ "$1" = "--pbs" ]] ; then
		echo "		Use: \"pbs_check.pl -t PBS_JOBID\" "
		echo "		To view stdout and stderr from each separate job while they are running" 
	fi
###########################################################################################################################################################		
###		Do some checking that files exist etc. (User should not modify)  
###		This sets the $DS variable

SET_DS "$FILENAMEINPUT"

###########################################################################################################################################################	
#################################                                                                                       ###################################
#################################       Part 2: Create the shell scripts to be run via the PBS batch system             ###################################
#################################       Users should modify WALLTIME and MEM dependent upon dataset size                ###################################
#################################       and NCPU to appropriate value for compute node cpu resources                    ###################################
#################################       Check: "Trinity RNA-seq Assembler Performance Optimisation" (Henschel 2012)     ###################################
#################################       for current best practice.                                                      ###################################
#################################       N.B. On busy clusters it may be best not to try to request                      ###################################
#################################       a full nodes resources. i.e. If 8 cores per node are present,                   ###################################
#################################       only request half of these.                                                     ###################################
###########################################################################################################################################################

###########################################################################################################################################################
##############################  Script 1: Write script to run Inchworm                                         ############################################
#############################   MEM should equal JFMEM, which is the amount of memory requested for Jellyfish ############################################

JOBNAME1="$JOBPREFIX"_p1
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" "$MEM_P1" "$NCPU_P1" "$PBSNODETYPE_P1" "$WALLTIME_P1" "$JOBNAME1" "$ACCOUNT" "$PBSUSER" "$MODTRINITY" "$JOBPREFIX")

F_WRITESCRIPT "$0" ""$TRINITYPBSPATH"/trinity_pbs.p1"

#############################################################################################################################################################
##############################  Script 2: Chrysalis::GraphFromFasta                                             #############################################
##############################  This script has a dependency on part 1 completion without error.                #############################################
# It would be good to force an exit(0) before ReadsToTranscripts after checkpoint file /chrysalis/GraphFromIwormFasta.finished is 
# written (i.e. add --no_run_readstotrans to Trinity and pass through to Chrysalis) as the script 3  can be started directly after.
# This may be less of an issue with the new (fast) version of Trinity::GraphFromFasta (since version 2012-06-08).


JOBNAME2="$JOBPREFIX"_p2 
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" "$MEM_P2" "$NCPU_P2" "$PBSNODETYPE_P2" "$WALLTIME_P2" "$JOBNAME2" "$ACCOUNT" "$PBSUSER" "$MODTRINITY" "$JOBPREFIX")

F_WRITESCRIPT "$0" ""$TRINITYPBSPATH"/trinity_pbs.p2"

###########################################################################################################################################################
##############################  Script 3: Script to run Chrysalis::ReadsToTranscripts                                                   ###################
##############################  ReadsToTranscripts can be slow due to reads from disk,                                                  ###################
##############################  This script has a dependency on part 2 completion with error.                                            ###################
##############################  Section script is skipped if enough time was given in Part 2.                                           ###################


JOBNAME3="$JOBPREFIX"_p3 
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" "$MEM_P3" "$NCPU_P3" "$PBSNODETYPE_P3" "$WALLTIME_P3" "$JOBNAME3" "$ACCOUNT" "$PBSUSER" "$MODTRINITY" "$JOBPREFIX")

F_WRITESCRIPT "$0" ""$TRINITYPBSPATH"/trinity_pbs.p3"


##########################################################################################################################################################
############################## Script 4a: Write script to call the Chrysalis QuantifyGraph array job                                    ##################
############################## N.B. SLOTLIMIT="%x" indicates 'slot limit' i.e. the number of concurrent jobs to execute in an array     ##################
############################## (Not available in PBSpro )                                                                               ##################
#### NB Disabling emails for arrays

#SLOTLIMIT="%64"   # available for PBS Torque
JOBNAME4="$JOBPREFIX"_p4a
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" 1gb 1 "$PBSNODETYPE_P4" 00:30:00 "$JOBNAME4" "$ACCOUNT" "$PBSUSER" "$MODTRINITY" "$JOBPREFIX" )

F_WRITESCRIPT "$0" ""$TRINITYPBSPATH"/trinity_pbs.p4a"

###########################################################################################################################################################
##############################   Script Array part 4b: Write script to be run as an array Job . Runs Chrysalis::QuantifyGraph            ##################
##############################   This scipt has a dependency on part 4a being run. If an array component fails it will email user.       ##################
##############################   Files named quantifyGraph_commands_X are written with subset of total commands (X is the array ID from the PBS system)  ##


JOBNAME4B="$JOBPREFIX"_p4b
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" "$MEM_P4" "$NCPU_P4" "$PBSNODETYPE_P4" "$WALLTIME_P4" "$JOBNAME4B" "$ACCOUNT" " " "$MODTRINITY" "$JOBPREFIX")

F_WRITESCRIPT "$0" ""$TRINITYPBSPATH"/trinity_pbs.p4b"


###########################################################################################################################################################
##############################   Script Array 5b: Write script to run Butterfly Array Job                                                ##################
##############################   This script has a dependency on part 4a 4b being run.                                                   ##################
##############################   Files named butterfly_commands_X are written with subset of total commands (X is the array ID from the PBS system) ###

JOBNAME5B="$JOBPREFIX"_p5b
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" "$MEM_P5" "$NCPU_P5" "$PBSNODETYPE_P5" "$WALLTIME_P5" "$JOBNAME5B" "$ACCOUNT" " " "$MODTRINITY" "$JOBPREFIX")

F_WRITESCRIPT "$0" ""$TRINITYPBSPATH"/trinity_pbs.p5b"

############################################################################################################################################################
###############                                                                                                                          ###################
###############  Part 3:  Write main control script that executes scripts that were created above.                                       ###################
###############                                                                                                                          ###################
############################################################################################################################################################
F_WRITESCRIPT "$0" ""$TRINITYPBSPATH"/trinity_pbs.cont"

###############  Run the script written in Part 3
###############  Checks to see what is current stage of calculation and executes scripts created in above code                           ###################
bash ""$JOBPREFIX"_run.sh"

exit 0

