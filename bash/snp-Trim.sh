#!/bin/bash
cat > bash/snp_call-trim.sh << EOF
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --array=1-$END
#SBATCH --job-name=TRIMMER
#SBATCH --output=slurm/trim-%A_%a.out

module load anaconda
module list
date

FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$2 ; }' $MASTER)

R1A=( \$( ls $DIRIN'/'\$FILEBASE*_R1_*.fastq.gz ))
R2A=( \$( ls $DIRIN'/'\$FILEBASE*_R2_*.fastq.gz ))

# Loop in parallel
count=\${#R1A[@]}
for i in \`seq 1 \$count\`
do
    R1=\${R1A[\$i-1]}
    O1='fastq_trim/'\$(basename \$R1 | sed 's/.fastq.gz//')'.trim.fastq.gz'
    R2=\${R2A[\$i-1]}
    O2='fastq_trim/'\$(basename \$R2 | sed 's/.fastq.gz//')'.trim.fastq.gz'
    echo "==>>" \$R1 \$R2
    cutadapt -q 20 -O 8 --minimum-length 40 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o \$O1 -p \$O2 \$R1 \$R2
done

EOF
