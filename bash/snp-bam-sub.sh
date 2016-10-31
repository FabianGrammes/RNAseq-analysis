#!/bin/bash

cat > bash/snp_call-BamSub.sh << EOF
#!/bin/bash
#SBATCH --job-name=BAMsub
#SBATCH --nodes=1
#SBATCH --array=1-$END%30
#SBATCH --output=slurm/snp_call-BamSub-%A_%a.out

module load samtools
module list
date

set -o errexit   # Exit if error occurs

## Set variables
FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$2 ; }' $MASTER)
BAM_in=\$(ls $BAMDIR/*.bam | grep \$FILEBASE)
BAM_sub=bam_sub/\$FILEBASE'-DeGenesSubset.dedupped.splitCig.bam'

samtools view -b -L $BED -o \$BAM_sub \$BAM_in

samtools index -b \$BAM_sub

## success report
awk -v var="\$FILEBASE" 'BEGIN {print var, "pass"}' >> report_bamsub.txt 

echo "FINISHED"

EOF
