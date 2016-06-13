## ==============================================================================

cat > bash/snp_call-HaploCall.sh << EOF
#!/bin/bash
#SBATCH --job-name=gVCF
#SBATCH -n 4
#SBATCH --nodes=1
#SBATCH --array=1-$END%20
#SBATCH --mem=15G
#SBATCH --output=slurm/snp_call-gVCF-%A_%a.out

module load gatk/3.5
module list
date

# set -o nounset   # Prevent unset variables from been used.
set -o errexit   # Exit if error occurs

## Set variables
FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$2 ; }' $MASTER)
BAM_cig=gatk/\$FILEBASE'-2pass-Aligned.out.sortedByCoord.dedupped.splitCig.bam'
VCF_out=gatk/\$(basename \$BAM_cig .bam).g.vcf.gz

echo "==> HaplotypeCaller: gVCF"

gatk -T HaplotypeCaller \
-R $GENFA \
-I \$BAM_cig \
-o \$VCF_out \
-ERC GVCF \
-nct 4 \
-dontUseSoftClippedBases \
-stand_call_conf 20.0 \
-stand_emit_conf 20.0 \
-variant_index_type LINEAR \
-variant_index_parameter 128000

echo '==> Done HaplotypeCaller'

EOF

## ==============================================================================

cat > bash/snp_call-myList.sh << EOF
#!/bin/bash
#SBATCH --job-name=SNPcall
#SBATCH -n 1
#SBATCH --output=slurm/snp_call-SNPcall-%A_%a.out

ls gatk/*.g.vcf > gVCFs.list

EOF

## ==============================================================================

cat > bash/snp_call-SNPcall.sh << EOF
#!/bin/bash
#SBATCH --job-name=SNPcall
#SBATCH -n 4
#SBATCH --nodes=1
#SBATCH --array=1-29
#SBATCH --mem=15G
#SBATCH --output=slurm/snp_call-SNPcall-%A_%a.out

module load gatk/3.5
module list
date

## chromosome id
CHROM=\$(printf "ssa%02d" \$SLURM_ARRAY_TASK_ID)
COMBI='gVCF/'\$CHROM'_combined.g.vcf.gz'
FILT='gVCF/'\$CHROM'_combined.filterd.g.vcf.gz'



echo '==> GenotypeGVCFs CHROM:' \$CHROM
-R $GENFA \
-V gVCFs.list \
-o \$COMBI \
-L \$CHROM \
-nt 4
echo '==> DONE: GenotypeGVCFs'


echo "==> VariantFiltration"

gatk -T VariantFiltration \
-R $GENFA \
-V \$COMBI \
-window 35 -cluster 3 \
-filterName FS -filter "FS>30.0" \
-filterName QD -filter "QD<2.0" \
-o \$FILT

echo '==> DONE: VariantFiltration'

EOF

## ==============================================================================
