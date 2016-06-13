
# ==============================================================================
# 4) Add read groups, sort, mark duplicates, and create index

cat > bash/snp_call-mDupl.sh << EOF
#!/bin/bash
#SBATCH --job-name=picard
#SBATCH -n 5
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --array=1-$END%20
#SBATCH --output=slurm/snp_call-picard-%A_%a.out

module load samtools picard
module list
date

#-------------------------------------------------------------------------------
FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$2 ; }' $MASTER)
BAM=star_2pass/\$FILEBASE'-2pass-Aligned.out.bam'
BAMC=star_2pass/\$FILEBASE'-2pass-Aligned.out.sortedByCoord.bam'
PIC=star_2pass/\$FILEBASE'-2pass-Aligned.out.sortedByCoord.dedupped.bam'
PIM=star_2pass/\$FILEBASE'.dedupped.metrics'
#-------------------------------------------------------------------------------

samtools sort -@5 -m 2G -o \$BAMC -T \$BAM'.temp' -O bam \$BAM
rm \$BAM
echo '==> Done sorting'


picard MarkDuplicates I=\$BAMC O=\$PIC CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=\$PIM
rm \$BAMC
echo '==> Done marking duplicates'
EOF

## ==============================================================================
## 5) Split'N'Trim and reassign mapping qualities

cat > bash/snp_call-splitNtrim.sh << EOF
#!/bin/bash
#SBATCH --job-name=splitNtrim
#SBATCH -n 1
#SBATCH --array=1-$END
#SBATCH --output=slurm/snp_call-splitNtrim-%A_%a.out

module load gatk/3.5
module list
date

FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$2 ; }' $MASTER)
BAM_in=star_2pass/\$FILEBASE'-2pass-Aligned.out.sortedByCoord.dedupped.bam'
BAM_cig=gatk/\$FILEBASE'-2pass-Aligned.out.sortedByCoord.dedupped.splitCig.bam'

gatk SplitNCigarReads -R $GENFA \
-I \$BAM_in \
-o \$BAM_cig \
-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
-U ALLOW_N_CIGAR_READS

echo '==> Done splitting cigars'
rm \$BAM_in
EOF

