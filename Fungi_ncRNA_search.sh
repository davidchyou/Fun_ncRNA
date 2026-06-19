SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
CURDIR=`pwd -P`
INFILE=$1 #"GCF_000146045.fna"
OUTDIR=$2 #"out"
CLASS=$3 #"Fungi"
RNAMMER_TRNASCAN=$4
CORE_CM_CUTOFF=$5
SNO_CM_CUTOFF=$6
SRP_CM_CUTOFF=$7
EXCLUDE_SN_VARIANTS=$8
ALL_RFAM_CMPATH=$9

CORE_CMDIR=$SCRIPTPATH/core_cm
CORE_CMPATH=$CORE_CMDIR/$CLASS.core.cm
SNO_CMPATH=$SCRIPTPATH/sno_u3_cm/sno.cm
U3_CMPATH=$SCRIPTPATH/sno_u3_cm/u3.cm
SRP_CMPATH=$SCRIPTPATH/SRP_Basidio/SRP_Basidio.cm
RNAMMER_OUT=$OUTDIR/RNAMMER
TRNA_OUT=$OUTDIR/TRNA
CORE_CM_OUT=$OUTDIR/CM_core
SNO_CM_OUT=$OUTDIR/CM_sno
U3_CM_OUT=$OUTDIR/CM_U3
SRP_CM_OUT=$OUTDIR/CM_SRP
ALL_RFAM_OUT=$OUTDIR/CM_all_rfam

if [ ! -e $ALL_RFAM_CMPATH ]
then
	ALL_RFAM_CMPATH=$SCRIPTPATH/Rfam/Rfam_f.cm #"/Volumes/scratch/brownlab/chrisbr/DB/Rfam/Rfam.cm"
fi

if [ ! -e $CORE_CMPATH ]
then
	CORE_CMPATH=$CORE_CMDIR/Fungi.core.cm
fi

rm -rf $OUTDIR
mkdir $OUTDIR

if [[ $RNAMMER_TRNASCAN -gt 0 ]]
then

rm -rf $RNAMMER_OUT
mkdir $RNAMMER_OUT
cp $INFILE $RNAMMER_OUT/input_genome.fna
cd $RNAMMER_OUT
#rnammer -S euk -m tsu,ssu,lsu -gff input_genome.rnammer.gff input_genome.fna -multi > input_genome.rnammer.log.txt 2> input_genome.rnammer.err.txt
perl $SCRIPTPATH/bin/rnammer/rnammer -S euk -m tsu,ssu,lsu -gff input_genome.rnammer.gff input_genome.fna -multi > input_genome.rnammer.log.txt 2> input_genome.rnammer.err.txt
rm -f input_genome.fna
cd $CURDIR
cat $RNAMMER_OUT/input_genome.rnammer.gff | grep -v '^#' | awk -F'\t' '{print $1, $2, $9, $4, $5, $6, $7, $8, $3}' OFS='\t' > $OUTDIR/input_genome.rnammer.gff
#cat $RNAMMER_OUT/input_genome.rnammer.gff | grep -v '^#' | cut -f 1,2,9,4,5,6,7,8,3 > $OUTDIR/input_genome.rnammer.gff

fi

if [[ $RNAMMER_TRNASCAN -gt 0 ]]
then

rm -rf $TRNA_OUT
mkdir $TRNA_OUT
cp $INFILE $TRNA_OUT/input_genome.fna
cd $TRNA_OUT
#tRNAscan-SE -j input_genome.trna.gff input_genome.fna > input_genome.trna.log.txt 2> input_genome.trna.err.txt
#/projects/health_sciences/bms/biochemistry/brown_lab/davidc/bioinfo_tools.v5.sif tRNAscan-SE -b input_genome.trna.bed input_genome.fna > input_genome.trna.log.txt 2> input_genome.trna.err.txt
tRNAscan-SE -b input_genome.trna.bed input_genome.fna > input_genome.trna.log.txt 2> input_genome.trna.err.txt
cat input_genome.trna.bed | awk -F'\t' '{print $1, "tRNAScan-SE", "tRNA", ($2>0)?$2:1, $3, $5, $6, ".", "ID="$4}' OFS='\t' | sort -k1,1 -k4,4n > input_genome.trna.gff 2> input_genome.trna.err.txt
rm -f input_genome.fna
cd $CURDIR
cat $TRNA_OUT/input_genome.trna.gff | grep -v '^#' > $OUTDIR/input_genome.trna.gff

fi

rm -rf $CORE_CM_OUT
mkdir $CORE_CM_OUT
cmsearch --cpu 128 --tblout $CORE_CM_OUT/input_genome.cmcore.tbl $CORE_CMPATH $INFILE > $CORE_CM_OUT/input_genome.cmcore.log.txt 2> $CORE_CM_OUT/input_genome.cmcore.err.txt
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl -T $CORE_CM_CUTOFF --all $CORE_CM_OUT/input_genome.cmcore.tbl | grep -v '^#' > $OUTDIR/input_genome.cmcore.gff

rm -rf $SNO_CM_OUT
mkdir $SNO_CM_OUT
cmsearch --cpu 128 --tblout $SNO_CM_OUT/input_genome.sno.tbl $SNO_CMPATH $INFILE > $SNO_CM_OUT/input_genome.sno.log.txt 2> $SNO_CM_OUT/input_genome.sno.err.txt
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl -T $SNO_CM_CUTOFF  --all $SNO_CM_OUT/input_genome.sno.tbl | grep -v '^#' > $OUTDIR/input_genome.sno.gff

rm -rf $U3_CM_OUT
mkdir $U3_CM_OUT
cmsearch --cpu 128 --tblout $U3_CM_OUT/input_genome.u3.tbl $U3_CMPATH $INFILE > $U3_CM_OUT/input_genome.u3.log.txt 2> $U3_CM_OUT/input_genome.u3.err.txt
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl -T $SNO_CM_CUTOFF --all $U3_CM_OUT/input_genome.u3.tbl | grep -v '^#' > $OUTDIR/input_genome.u3.gff

rm -rf $SRP_CM_OUT
mkdir $SRP_CM_OUT
cmsearch --cpu 128 --tblout $SRP_CM_OUT/input_genome.srp.tbl $SRP_CMPATH $INFILE > $SRP_CM_OUT/input_genome.srp.log.txt 2> $SRP_CM_OUT/input_genome.srp.err.txt
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl -T $SRP_CM_CUTOFF --all $SRP_CM_OUT/input_genome.srp.tbl | grep -v '^#' > $OUTDIR/input_genome.srp.gff

rm -rf $ALL_RFAM_OUT
mkdir $ALL_RFAM_OUT
cmsearch --cpu 128 --tblout $ALL_RFAM_OUT/input_genome.rfam.tbl --cut_tc $ALL_RFAM_CMPATH $INFILE > $ALL_RFAM_OUT/input_genome.rfam.log.txt 2> $ALL_RFAM_OUT/input_genome.rfam.err.txt
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl -E 0.001 --all $ALL_RFAM_OUT/input_genome.rfam.tbl | grep -v '^#' > $OUTDIR/input_genome.rfam.gff

if [[ $EXCLUDE_SN_VARIANTS -gt 0 ]]
then

rm -f $OUTDIR/input_genome.cmcore.gff
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl -T $CORE_CM_CUTOFF --all $CORE_CM_OUT/input_genome.cmcore.tbl | grep -v '^#' | awk -F'\t' '(!($3~/U11/||$3~/U12/||$3~/U4atac/||$3~/U6atac/))' > $OUTDIR/input_genome.cmcore.gff

fi

cat $OUTDIR/input_genome.*.gff | sort -k1,1 -k4,4n > $OUTDIR/input_genome.all.gff
bedtools cluster -i $OUTDIR/input_genome.all.gff -d 0 -s | sort -k10,10 -k6,6nr | sort -k10,10 -u | sort -k1,1 -k4,4n | cut -f 1-9 > $OUTDIR/input_genome.all.nr.gff

