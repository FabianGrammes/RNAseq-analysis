#!/bin/bash

## ==============================================================================

cat > bash/snp_call-HaploCall.sh << EOF
#!/bin/bash
#SBATCH --job-name=gVCF
#SBATCH -n 1
#SBATCH --array=1-$END%48
#SBATCH --mem=5G
#SBATCH --output=slurm/snp_call-gVCF-%A_%a.out

module load gatk/3.5
module list
date

## Set variables
FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$2 ; }' $MASTER)
BAM=\$(ls $BAMDIR/*.bam | grep \$FILEBASE)
VCF=vcf/\$(basename \$BAM .bam).g.vcf

echo "==> HaplotypeCaller: gVCF"

gatk -T HaplotypeCaller \
-R $GENFA \
-I \$BAM \
-o \$VCF \
-ERC GVCF \
-dontUseSoftClippedBases \
-stand_call_conf 20.0 \
-stand_emit_conf 20.0 \
-variant_index_type LINEAR \
-variant_index_parameter 128000

echo '==> Done HaplotypeCaller'

EOF

cat > bash/snp_call-myList.sh << EOF
#!/bin/bash
#SBATCH --job-name=List
#SBATCH -n 1
#SBATCH --output=slurm/snp_call-gVCF-listl-%A_%a.out

ls vcf/*.g.vcf > gVCFs.list

EOF


cat > bash/snp_call-SNPcall.sh << EOF
#!/bin/bash
#SBATCH --job-name=SNPcall
#SBATCH -n 4
#SBATCH --nodes=1
#SBATCH --mem=15G
#SBATCH --output=slurm/snp_call-SNPcall-%A_%a.out

module load gatk/3.5
module list
date

echo '==> GenotypeGVCFs'
gatk -T GenotypeGVCFs -R $GENFA \
-V gVCFs.list \
-o vcf/VCFgenotypes.raw.vcf \
-nt 4
echo '==> DONE: GenotypeGVCFs'

echo "==> VariantFiltration"
gatk -T VariantFiltration -R $GENFA \
-V vcf/VCFgenotypes.raw.vcf \
-window 35 -cluster 3 \
-filterName FS -filter "FS>30.0" \
-filterName QD -filter "QD<2.0" \
-o vcf/VCFgenotypes.filtered.vcf
echo '==> DONE: VariantFiltration'

EOF
