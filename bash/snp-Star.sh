#!/bin/bash
## ==============================================================================
## 3) ReRun STAR

cat > bash/snp_call-Star.sh << EOF
#!/bin/bash
#SBATCH --job-name=STAR2nd
#SBATCH -n 15
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --array=1-$END%8
#SBATCH --output=slurm/snp_call-Star-%A_%a.out

module load star
module list
date

SAMPLE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$1 ; }' $MASTER)
FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$2 ; }' $MASTER)

R1A=( \$( ls $FADIR\$FILEBASE*_R1_*.fastq.gz ))
R2A=( \$( ls $FADIR\$FILEBASE*_R2_*.fastq.gz ))

R1=\$(printf ",%s" "\${R1A[@]}")
R2=\$(printf ",%s" "\${R2A[@]}")
R1=\${R1:1}
R2=\${R2:1}

#-------------------------------------------------------------------------------
# BAM read groups:
RGa=() # Array to hold the read group string
for ((i=0; i<\${#R1A[@]}; i++))
do
    printf "%s\t%s\n" \$( printf "ID:L%03d" \$((\$i+1)) ) \$( basename \${R1A[\$i]} .trim.fastq.gz ) >> ReadGroup_summary_2pass.txt
    RGa+=\$(printf "ID:L%03d PL:illumina LB:\$SAMPLE SM:\$SAMPLE , " \$((\$i+1)))
done

RG=\$( printf "%s" "\${RGa[@]}" )
RG=\$(echo \$RG | sed 's/ ,\$//g' )
#-------------------------------------------------------------------------------

OUT=star_2pass/\$FILEBASE'-2pass-'

echo "Running  --> " \$R1 \$R2

# Run STAR
$STAR --limitGenomeGenerateRAM 62000000000 \
--genomeDir $GENDIR \
--readFilesCommand zcat \
--readFilesIn \$R1 \$R2 \
--outFileNamePrefix \$OUT \
--outSAMmode Full \
--outSAMtype BAM Unsorted \
--runThreadN 20 \
--readMatesLengthsIn NotEqual \
--outSAMattrRGline \$RG

echo "FILE --> " \$OUT " PROCESSED"

EOF

## ==============================================================================

cat > bash/snp_call-StarCheck.sh << EOF
#!/bin/bash
#SBATCH --job-name=STARstat
#SBATCH -n 1
#SBATCH --output=slurm/snp_call-StarStats%j.out

module load R
module list
date

cd star_2pass
R CMD BATCH /mnt/users/fabig/cluster_pipelines/RnaMapping/helper_scripts/STAR_Log.R
cd ..

EOF


# ==============================================================================
