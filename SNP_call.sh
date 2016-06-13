
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
	RLEN="$2"
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
	RTarray=(makeIdx mapp ReRun BamSub)
	if ! [[ "${RTarray[@]}" =~ "$RTYPE" ]]
	then
	    echo "ERROR: Unknown --RunType:" $RTYP
	    exit 1
	fi
	;;
esac
shift # past argument or value
done


if [ "$RTYP" == "makeIdx" ]
then
    	if [ -z "$RLEN" ]; then
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
	echo 'Read length =' $RLEN
	echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
	;;
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
	    mkdir -p {slurm,bash,star_2pass,mapp_summary,gatk,gVCF}
	    
	    ## Determine which STAR command to use (e.g STAR/STARlong)
	    if [ "$RLEN" == "300" ]
	    then
		STAR="STARlong"
	    else
		STAR="STAR"
	    fi

	    ## echo all input variables
	    echo '-----------------------'
	    echo 'Location of .fastq files =' $FADIR
	    echo 'Genome .fasta =' $GENFA
	    echo '.gtf file =' $GTF
	    echo 'STAR command =' $STAR
	    echo '-----------------------'    
	    
	    ;;
	"ReRun")
	    ## Create the folder tree if it does not exist
	    mkdir -p {slurm,bash,star_2pass,mapp_summary,gatk,gVCF}

	    	    ## echo all input variables
	    echo '-----------------------'
	    echo 'Location of .fastq files =' $FADIR
	    echo 'Genome .fasta =' $GENFA
	    echo '.gtf file =' $GTF
	    echo 'STAR command =' $STAR
	    echo '-----------------------'
	    
	    ;;
	"BamSub")
	    if [ ! -d "$BAMDIR"  ] || [ -z "$BAMDIR" ]
	    then
		echo "ERROR: Misssing or could not be found" $BAMDIR
		exit 1
	    fi

	    if [ ! -f "$BED"  ] || [ -z "$BED" ]
	    then
		echo "ERROR: Misssing or could not be found" $BED
		exit 1
	    fi
	    	    
	    ## Create the folder tree if it does not exist
	    mkdir -p {slurm,bash,bam_sub,gatk,vcf}

	    ## echo all input variables
	    echo '-----------------------'
	    echo 'Location of .bam files =' $BAMDIR
	    echo 'Genome .fasta =' $GENFA
	    echo '.gtf file =' $GTF
	    echo 'STAR command =' $STAR
	    echo '-----------------------'
	    ;;
	
    esac


fi
		    
#-------------------------------------------------------------------------------
# Path to the R scripts
SCRIPTPATH=$(readlink -f "$0") # finds the path where the script resides
HELPPATH=$(dirname "$SCRIPTPATH")/helper_scripts
SHPATH=$(dirname "$SCRIPTPATH")/bash

echo $SCRIPTPATH
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

##==============================================================================
## PRINTING THE BASH SCRIPTS (executing BASH in BASH: Variables are passed)
##==============================================================================

if [ "$RTYP" == "makeIdx"  ]  ## Only use the commands to make the STAR index
then
    sh $SHPATH/snp-StarIndex.sh
fi

if [ "$RTYP" == "ReRun"  ]  ## Only use the commands to make the STAR index
then
    sh $SHPATH/snp-Trim.sh
fi

if [ "$RTYP" == "mapp"  ] || [ "$RTYP" == "ReRun" ]
then
    ## ReRun STAR
    sh $SHPATH/snp-Star.sh
    ## Run GATK (BAM) pipeline
    sh $SHPATH/snp-gatk-bam.sh
    ## Run GATK (BAM) pipeline
    sh $SHPATH/snp-gatk-vcf.sh
fi

if [ "$RTYP" == "BamSub" ]
then
    
    

fi
   
# ==============================================================================

# SCRIPT SUBMISSION

if [ "$EXECUTE" != "no" ]
then
    if [ "$RTYP" == "makeIdx" ]
    then
	#-------------------------------------------------------------------------------
	# 1) splice junctions
	command='sbatch bash/snp_call-SJ.sh'
	SJjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' Splice junctions'
	echo ' slurm ID:' $SJjob

	#-------------------------------------------------------------------------------
	# 2) STAR index
	command="sbatch --dependency=afterok:$SJjob bash/snp_call-StarIdx.sh"
	IDXjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' Making STARindex'
	echo ' slurm ID' $IDXjob
    fi
    if [ "$RTYP" == "ReRun" ]
    then
	mkdir -p fastq_trim
	command="sbatch bash/snp_call-trim.sh"
	TRjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' Trimming'
	echo ' slurm ID:' $SJjob
    fi
    if [ "$RTYP" == "mapp"  ] || [ "$RTYP" == "ReRun" ]
    then
	#-------------------------------------------------------------------------------
	# 1) 2nd Round STAR
	if [ -z "$TRjob" ] # if no Trim job exists
	then
	    command="sbatch bash/snp_call-Star.sh"
	else
	    command="sbatch --dependency=afterok:$TRjob bash/sbatch-star.sh"
	fi

	STARjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' 2nd round mapping'
	echo ' slurm ID' $STARjob

	#-------------------------------------------------------------------------------
	command="sbatch --dependency=afterok:$STARjob bash/snp_call-StarCheck.sh"
	CHECKjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' star check processing'
	echo ' slurm ID' $CHECKjob

	#-------------------------------------------------------------------------------
	# 1) Picard
	command="sbatch --dependency=afterok:$STARjob bash/snp_call-mDupl.sh"
	PICjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' picard processing'
	echo ' slurm ID' $PICjob


	#-------------------------------------------------------------------------------
	command="sbatch --dependency=afterok:$PICjob bash/snp_call-splitNtrim.sh"
	SPLITjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' splitNtrim processing'
	echo ' slurm ID' $SPLITjob

        #-------------------------------------------------------------------------------
	command="sbatch --dependency=afterok:$SPLITjob bash/snp_call-HaploCall.sh"
	VARjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' slurm ID' $VARjob

        command="sbatch --dependency=afterok:$VARjob bash/snp_call-myList.sh"
	
        #-------------------------------------------------------------------------------
	command="sbatch --dependency=afterok:$VARjob bash/snp_call-SNPcall.sh"
	VARjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' Variant filtering'
	echo ' slurm ID' $VARjob
	
    fi
else
    echo '=================='
    echo 'Nothing submitted!'
    echo '=================='
fi
