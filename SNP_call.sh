#!/bin/bash

echo ''
echo $(date)
echo 'SNP_call.sh Version 1.1.0'

#-------------------------------------------------------------------------------
# Read script arguments and check if files/dirs exist

while [[ $# > 0 ]]
do
key="$1"

case $key in
    --sj)
	SJ="$2" # File path to all folders containing splice junction files (STAR folder)
	shift   # if more than 1 folder seperate by comma.
	;;
    --GenFa)      # OPTIONAL: path to genome
	GENFA="$2"
	shift # past argument
	# check if .fasta file exists
	if [ ! -f "$GENFA" ]; then
	    echo 'ERROR: File' $GENFA 'Does not exist!'
	    exit 1
	else
	    echo 'FOUND: File' $GENFA
	fi
	;;
    --GenDir)      # OPTIONAL: path to the newgenome
	GENDIR="$2"
	shift # past argument
	;;
    --FaDir)      # .fastq path
	FADIR="$2"
	shift # past argument
	;;
    --bamDir)     # OPTIONAL: Path to Bam files
	BAMDIR="$2"
	shift
	;;
    --read)      # Read length
	READLEN="$2"
	shift # past argument
	;;
    --gtf)            # OPTIONAL: path to .gtf
	GTF="$2"
	shift # past argument
	# check if .gtf file exists
	if [ ! -f "$GTF" ]; then
	    echo 'ERROR: File' $GTF 'Does not exist!'
	    exit 1
	else
	    echo 'FOUND: File' $GTF
	fi
	;;
    -m|--mastersheet) # file name of the mastersheet
	MASTER="$2"
	shift # past argument
	# check if MASTER file exists
	;;
    --execute)        # Only used for testing; use --execute no
	EXECUTE="$2"
	shift # past argument
	;;
    --bed)        # Bed file for subsetting .bam
	BED="$2"
	shift # past argument
	;;
    --RawDir)
	DIRIN="$2"
	shift # past argument
	;;
    --RunType)
        RTYPE="$2"  # Accepts makeIdx/mapp
	shift
	;;
esac
shift # past argument or value
done

#===============================================================================

if [ -z "$READLEN" ]; then
    READLEN="125"
    echo "SETTING: Read length to" $READLEN	    
else
    echo "FOUND: Read length" $READLEN
fi


if [ "$RTYPE" = "makeIdx" ]
then
    if [ -z "$READLEN" ]; then
	echo 'ERROR: You have to specify read length at --read'
	exit 1
    fi
    mkdir -p {slurm,bash,$GENDIR}
    # echo all input variables
    echo ''
    echo 'Executing' $EXECUTE
    echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    echo 'Input arguments:'
    echo '-----------------------'
    echo 'Genome .fasta =' $GENFA
    echo '.gtf file =' $GTF
    echo 'Read length =' $READLEN
    echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
	
else
    ## 1st check for MASTER
    if [ ! -f "$MASTER" ]
    then
	echo 'ERROR: File' $MASTER 'Does not exist!'
	exit 1
    else
	echo 'FOUND: File' $MASTER
    fi
    
    ## 2ndCheck line terminators in MASTER
    if cat -v $MASTER | grep -q '\^M' 
    then
	echo 'Converting line terminators'
	sed 's/\r/\n/g' $MASTER > sheet.tmp
	mv sheet.tmp $MASTER
    else
	echo 'Line terminators are correct'
    fi

    ## Determine which STAR command to use (e.g STAR/STARlong) 
    if [ "$READLEN" == "300" ] 
    then
	STAR="STARlong"
    else
	STAR="STAR"
    fi
    
    ## Determine the number of samples
    END=$(sed '1d' $MASTER | wc -l) # skip hearder line

    ##
    echo ''
    echo 'Input arguments:'
    echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    echo 'Number of samples= ' $END
    echo 'FIRST sample:' $(awk ' NR=='2' {OFS="\t"; print; }' $MASTER)
    echo 'LAST sample:' $(awk ' NR=='$END+1' {OFS="\t"; print; }' $MASTER)

    case "$RTYPE" in
	"mapp")
	    ## Create the folder tree if it does not exist
	    mkdir -p {slurm,bash,star_2pass,mapp_summary,vcf}
	    
	    ## echo all input variables
	    echo '-----------------------'
	    echo 'mapping'
	    echo 'Location of .fastq files =' $FADIR
	    echo 'Genome .fasta =' $GENFA
	    echo '.gtf file =' $GTF
	    echo 'STAR command =' $STAR
	    echo '-----------------------'    	    
	    ;;
	"ReRun")
	    ## Create the folder tree if it does not exist
	    mkdir -p {slurm,bash,star_2pass,mapp_summary,vcf,fastq_trim}

	    	    ## echo all input variables
	    echo '-----------------------'
	    echo 'ReRunning raw .fastq'
	    echo 'Location of .fastq files =' $FADIR
	    echo 'Genome .fasta =' $GENFA
	    echo '.gtf file =' $GTF
	    echo 'STAR command =' $STAR
	    echo '-----------------------'
	    
	    ;;
	"BamSub")
	    if [ ! -d "$BAMDIR"  ] || [ -z "$BAMDIR" ]
	    then
		echo "ERROR: .bam dir Misssing or could not be found" $BAMDIR
		exit 1
	    fi

	    if [ ! -f "$BED"  ] || [ -z "$BED" ]
	    then
		echo "ERROR: .bed file Misssing or could not be found" $BED
		exit 1
	    fi
	    	    
	    ## Create the folder tree if it does not exist
	    mkdir -p {slurm,bash,bam_sub,vcf}

	    ## echo all input variables
	    echo '-----------------------'
	    echo 'Using .bam subset'
	    echo 'Location of .bam files =' $BAMDIR
	    echo '.bed file =' $BED
	    echo 'Genome .fasta =' $GENFA
	    echo '.gtf file =' $GTF
	    echo 'STAR command =' $STAR
	    echo '-----------------------'
	    ;;
	*)
	    echo "ERROR: Unknown --RunType:" $RTYP
	    exit 1
	    ;;
    esac
fi
		    
#-------------------------------------------------------------------------------
# Path to the R scripts
SCRIPTPATH=$(readlink -f "$0") # finds the path where the script resides
HELPPATH=$(dirname "$SCRIPTPATH")/helper_scripts
SHPATH=$(dirname "$SCRIPTPATH")/bash

echo $SCRIPTPATH
echo $SHPATH
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

##==============================================================================
## PRINTING THE BASH SCRIPTS 2 FILE
##==============================================================================
## sourcing the scriupts ensures that variables are passed down

case "$RTYPE" in
    "makeIdx")
	## make genome index
	. $SHPATH/snp-StarIndex.sh
	;;
    "ReRun")
	. $SHPATH/snp-Trim.sh
	## ReRun STAR
	. $SHPATH/snp-Star.sh
	## Run GATK (BAM) pipeline
	. $SHPATH/snp-gatk-bam.sh
	## Run GATK (BAM) pipeline
	. $SHPATH/snp-gatk-vcf.sh
	;;
    "mapp")
	## ReRun STAR
l	. $SHPATH/snp-Star.sh
	## Run GATK (BAM) pipeline
	. $SHPATH/snp-gatk-bam.sh
	## Run GATK (BAM) pipeline
	. $SHPATH/snp-gatk-vcf.sh
	;;
    "BamSub")
	## Subset .bam by .bed
	. $SHPATH/snp-bam-sub.sh
	BAMDIR=bam_sub
	## Run GATK (BAM) pipeline
	. $SHPATH/snp-gatk-vcf.sh
	;;
esac
   
##==============================================================================
## SCRIPT SUBMISSION (INCLUDING DEPENDCY)
##==============================================================================

if [ "$EXECUTE" == "no" ]
then
    echo ''
    echo '=================='
    echo 'Nothing submitted!'
    echo '=================='
else
    case "$RTYPE" in
	"makeIdx")
	    ## 1) splice junctions
	    command='sbatch bash/snp_call-SJ.sh'
	    job1=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' Splice junctions:' $job1
	    ##-------------------------------------------------------------------------------
	    ## 2) STAR index
	    command="sbatch --dependency=afterok:$job1 bash/snp_call-StarIdx.sh"
	    job2=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' Making STARindex:' $job2
	    ;;
	"ReRun")
	    command="sbatch bash/snp_call-trim.sh"
	    job1=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' Trimming:' $job1
	    # job 2
	    command="sbatch --dependency=afterok:$job1 bash/sbatch-star.sh"
	    job2=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' 2nd round mapping' $job2
	    # job 2.1
	    command="sbatch --dependency=afterok:$job2 bash/snp_call-StarCheck.sh"
	    job21=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' star check processing' $job21
	    ## job 3 Picard
	    command="sbatch --dependency=afterok:$job2 bash/snp_call-mDupl.sh"
	    job3=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' picard processing' $job3
	    ## job 4
	    command="sbatch --dependency=afterok:$job4 bash/snp_call-splitNtrim.sh"
	    job5=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' splitNtrim processing' $job5
	    ## job 6
	    command="sbatch --dependency=afterok:$job5 bash/snp_call-HaploCall.sh"
	    job6=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' HaploTypeCalling' $job6
	    ## job 7
  	    command="sbatch --dependency=afterok:$job6 bash/snp_call-SNPcall.sh"
	    job7=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' Variant Joining/filtering' $job7
	    ;;
	"mapp")
	    ;;
	"BamSub")
	    
	    command="sbatch bash/snp_call-BamSub.sh"
	    job5=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' BAM subsetting' $job5

	    ## job 6
	    command="sbatch --dependency=afterok:$job5 bash/snp_call-HaploCall.sh"
	    job6=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' HaploTypeCalling' $job6

	    ## job 7
  	    command="sbatch --dependency=afterok:$job6 bash/snp_call-myList.sh"
	    job7=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' Listing g.vcf' $job7	    

	    ## job 8
  	    command="sbatch --dependency=afterok:$job7 bash/snp_call-SNPcall.sh"
	    job8=$($command | awk ' { print $4 }')
	    echo '---------------'
	    echo ' Variant Joining/filtering' $job8
	    ;;
    esac

    echo ' '
    echo '==>> ALL SUBMITTED <<=='
fi

