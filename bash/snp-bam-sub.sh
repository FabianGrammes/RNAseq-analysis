cat > bash/snp_call-BamSub.sh << EOF
#!/bin/bash
#SBATCH --job-name=BAMsub
#SBATCH --nodes=1
--array=1-$END%20
#SBATCH --output=slurm/snp_call-BamSub-%A_%a.out

module load samtools
module list
date

# set -o nounset   # Prevent unset variables from been used.
set -o errexit   # Exit if error occurs

## Set variables
FILEBASE=$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print $2 ; }' S1-HK_rerunn.txt)
BAM_cig=gatk/\$FILEBASE'-2pass-Aligned.out.sortedByCoord.dedupped.splitCig.bam'

BAM_sub=gatk/\$FILEBASE'-DeGenesSubset.dedupped.splitCig.bam'

BED=HRvsLR_DE-loci.bed

samtools view -b -L \$BED -o \$BAM_sub \$BAM_cig

samtools index -b \$BAM_sub

EOF
