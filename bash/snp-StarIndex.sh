#!/bin/bash
## ==============================================================================
## 1) Junctions
cat > bash/snp_call-SJ.sh << EOF
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --job-name=JUNCTIONS
#SBATCH --output=slurm/snp_call-SJ-%A.out

module list
date

# convert 2 array
SJa=($SJ)
SJa=\$(echo \$SJa | tr "," "\n")

# Join all SJ files
for i in \$SJa
do
    echo '++>>' $i
    for ii in \$(ls \$i/*SJ.out.tab)
    do
	echo \$ii
	awk '\$7>1' \$ii >> $GENDIR/SJ_all.tab
    done
done

# Filter the joined file
awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if(\$5>0){print \$1,\$2,\$3,strChar[\$4]}}' $GENDIR/SJ_all.tab | sort | uniq > $GENDIR/SJ_in.tab

cat $GENDIR/SJ_in.tab | grep "^ssa" > $GENDIR/SJ_in_genome.tab # Only utilize genome 

echo '==>>FINISHED'
EOF

## ==============================================================================
## 2) Remake the STAR genome index considering the new splice junctions

cat > bash/snp_call-StarIdx.sh << EOF
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --job-name=StarIdx
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --mem=60G
#SBATCH --output=slurm/snp_call-StarIdx-%A.out

module load star
module list
date


STAR --runMode genomeGenerate \
--runThreadN 10 \
--limitGenomeGenerateRAM 62000000000 \
--genomeChrBinNbits 12 \
--genomeDir $GENDIR \
--sjdbGTFtagExonParentTranscript Parent \
--genomeFastaFiles $GENFA \
--sjdbOverhang \$(($READLEN-1)) \
--sjdbFileChrStartEnd $GENDIR/SJ_in_genome.tab \
--sjdbGTFfile $GTF

echo '==>>FINISHED'

EOF
## ==============================================================================
