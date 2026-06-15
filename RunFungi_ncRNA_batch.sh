#!/bin/bash
#SBATCH --job-name=Fun_ncRNA
#SBATCH --account=chyte02p
#SBATCH --partition=aoraki
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=40
#SBATCH --time=4-00:00:00

DATA_DIR="/projects/health_sciences/bms/biochemistry/brown_lab/DB/Trichocosom_cleaned"
OUTDIR="/projects/health_sciences/bms/biochemistry/brown_lab/DB/Trichocosom_ncRNA"
CLASS="Sordariomycetes"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate /projects/health_sciences/bms/biochemistry/brown_lab/davidc/Fun_ncRNA_dependencies

for FUNGI in $(find $DATA_DIR -type f -name '*.fasta' | perl -pe 's/\S+\/([^\/]+).fasta$/\1/g')
do

SUBDIR=$OUTDIR/$FUNGI

perl ./Fungi_ncRNA.pl \
-fna $DATA_DIR/$FUNGI.fasta \
-gff $DATA_DIR/$FUNGI.fixed.gff \
-fungal_class $CLASS \
-out $SUBDIR \
-run_rnammer_trnascan \
-merge_separately \
-gffcmp \
-gcf $FUNGI

done

conda deactivate
