SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
CURDIR=`pwd -P`
BASEDIR=$1
OUTDIR=$2
RNAMMER_TRNASCAN=$3
CORE_CM_CUTOFF=$4
SNO_CM_CUTOFF=$5
SRP_CM_CUTOFF=$6
EXCLUDE_SN_VARIANTS=$7

RNAMMER_OUT=$OUTDIR/RNAMMER
TRNA_OUT=$OUTDIR/TRNA
CORE_CM_OUT=$OUTDIR/CM_core
SNO_CM_OUT=$OUTDIR/CM_sno
U3_CM_OUT=$OUTDIR/CM_U3
SRP_CM_OUT=$OUTDIR/CM_SRP
ALL_RFAM_OUT=$OUTDIR/CM_all_rfam

rm -rf $OUTDIR
mkdir $OUTDIR

if [ -e $BASEDIR/RNAMMER/input_genome.rnammer.gff ]
then

cp -r $BASEDIR/RNAMMER $RNAMMER_OUT
cat $RNAMMER_OUT/input_genome.rnammer.gff | grep -v '^#' | awk -F'\t' '{print $1, $2, $9, $4, $5, $6, $7, $8, $3}' 'OFS=\t' > $OUTDIR/input_genome.rnammer.gff
#cat $RNAMMER_OUT/input_genome.rnammer.gff | grep -v '^#' | cut -f 1,2,9,4,5,6,7,8,3 > $OUTDIR/input_genome.rnammer.gff

else

if [[ $RNAMMER -gt 0 ]]
then

rm -rf $RNAMMER_OUT
mkdir $RNAMMER_OUT
cp $INFILE $RNAMMER_OUT/input_genome.fna
cd $RNAMMER_OUT
perl $SCRIPTPATH/bin/rnammer/rnammer -S euk -m tsu,ssu,lsu -gff input_genome.rnammer.gff input_genome.fna -multi > input_genome.rnammer.log.txt 2> input_genome.rnammer.err.txt
#rnammer -S euk -m tsu,ssu,lsu -gff input_genome.rnammer.gff input_genome.fna -multi > input_genome.rnammer.log.txt 2> input_genome.rnammer.err.txt
rm -f input_genome.fna
cd $CURDIR
cat $RNAMMER_OUT/input_genome.rnammer.gff | grep -v '^#' | awk -F'\t' '{print $1, $2, $9, $4, $5, $6, $7, $8, $3}' 'OFS=\t' > $OUTDIR/input_genome.rnammer.gff
#cat $RNAMMER_OUT/input_genome.rnammer.gff | grep -v '^#' | cut -f 1,2,9,4,5,6,7,8,3 > $OUTDIR/input_genome.rnammer.gff

fi

fi

if [ -e $BASEDIR/TRNA/input_genome.trna.gff ]
then

cp -r $BASEDIR/TRNA $TRNA_OUT
cat $TRNA_OUT/input_genome.trna.gff | grep -v '^#' > $OUTDIR/input_genome.trna.gff

else

if [[ $TRNASCAN -gt 0 ]]
then

rm -rf $TRNA_OUT
mkdir $TRNA_OUT
cp $INFILE $TRNA_OUT/input_genome.fna
cd $TRNA_OUT
#tRNAscan-SE -j input_genome.trna.gff input_genome.fna > input_genome.trna.log.txt 2> input_genome.trna.err.txt
tRNAscan-SE -b input_genome.trna.bed input_genome.fna > input_genome.trna.log.txt 2> input_genome.trna.err.txt
cat input_genome.trna.bed | awk -F'\t' '{print $1, "tRNAScan-SE", "tRNA", ($2>0)?$2:1, $3, $5, $6, ".", "ID="$4}' OFS='\t' | sort -k1,1 -k4,4n > input_genome.trna.gff 2> input_genome.trna.err.txt
rm -f input_genome.fna
cd $CURDIR
cat $TRNA_OUT/input_genome.trna.gff | grep -v '^#' > $OUTDIR/input_genome.trna.gff

fi

fi

if [ -e $BASEDIR/CM_core/input_genome.cmcore.tbl ]
then

cp -r $BASEDIR/CM_core $CORE_CM_OUT
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl  -T $CORE_CM_CUTOFF --all $CORE_CM_OUT/input_genome.cmcore.tbl | grep -v '^#' > $OUTDIR/input_genome.cmcore.gff

else

touch $OUTDIR/input_genome.cmcore.gff

fi

if [ -e $BASEDIR/CM_sno/input_genome.sno.tbl ]
then

cp -r $BASEDIR/CM_sno $SNO_CM_OUT
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl  -T $SNO_CM_CUTOFF  --all $SNO_CM_OUT/input_genome.sno.tbl | grep -v '^#' > $OUTDIR/input_genome.sno.gff

else

touch $OUTDIR/input_genome.sno.gff

fi

if [ -e $BASEDIR/CM_U3/input_genome.u3.tbl ]
then

cp -r $BASEDIR/CM_U3 $U3_CM_OUT
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl  -T $SNO_CM_CUTOFF --all $U3_CM_OUT/input_genome.u3.tbl | grep -v '^#' > $OUTDIR/input_genome.u3.gff

else

touch $OUTDIR/input_genome.u3.gff

fi

if [ -e $BASEDIR/CM_SRP/input_genome.srp.tbl ]
then

cp -r $BASEDIR/CM_SRP $SRP_CM_OUT
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl -T $SRP_CM_CUTOFF --all $SRP_CM_OUT/input_genome.srp.tbl | grep -v '^#' > $OUTDIR/input_genome.srp.gff

else

touch $OUTDIR/input_genome.srp.gff

fi

if [ -e $BASEDIR/CM_all_rfam/input_genome.rfam.tbl ]
then

cp -r $BASEDIR/CM_all_rfam $ALL_RFAM_OUT
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl -E 0.001 --all $ALL_RFAM_OUT/input_genome.rfam.tbl | grep -v '^#' > $OUTDIR/input_genome.rfam.gff

else

touch $OUTDIR/input_genome.rfam.gff

fi

if [[ $EXCLUDE_SN_VARIANTS -gt 0 ]]
then

rm -f $OUTDIR/input_genome.cmcore.gff
perl $SCRIPTPATH/bin/infernal-tblout2gff.pl  -T $CORE_CM_CUTOFF --all $CORE_CM_OUT/input_genome.cmcore.tbl | grep -v '^#' | awk -F'\t' '(!($3~/U11/||$3~/U12/||$3~/U4atac/||$3~/U6atac/))' > $OUTDIR/input_genome.cmcore.gff

fi

cat $OUTDIR/input_genome.*.gff | sort -k1,1 -k4,4n > $OUTDIR/input_genome.all.gff
bedtools cluster -i $OUTDIR/input_genome.all.gff -d 0 -s | sort -k10,10 -k6,6nr -u | sort -k1,1 -k4,4n | cut -f 1-9 > $OUTDIR/input_genome.all.nr.gff

