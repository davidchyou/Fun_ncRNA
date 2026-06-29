#!/bin/bash
#SBATCH --job-name=Fun_ncRNA
#SBATCH --account=chyte02p
#SBATCH --partition=aoraki
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=40
#SBATCH --time=1-00:00:00

#module load perl-bioperl/1.7.6-c6yd4fm

source ~/miniforge3/etc/profile.d/conda.sh
conda activate /projects/health_sciences/bms/biochemistry/brown_lab/davidc/Fun_ncRNA_dependencies

fna_path="./data/GCF_024613125.1_Pisori2_genomic.fna"
gff_path="./data/GCF_024613125.1_Pisori2_genomic.gff"
class="Agaricomycetes"
outdir="example_output"

rm -rf $outdir

perl ./Fungi_ncRNA.pl \
-fna $fna_path \
-gff $gff_path \
-fungal_class $class \
-out $outdir \
-run_rnammer_trnascan \
-merge_separately \
-gffcmp

conda deactivate