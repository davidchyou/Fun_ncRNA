SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
CURDIR=`pwd -P`
BASEDIR=$1 #"/Volumes/scratch/brownlab/davidc/RefPref_fungi_JGI_default/*"
OUTDIR=$2 #"test_out"

rm -rf $OUTDIR
mkdir $OUTDIR

TRNA_GFF=$OUTDIR/input_genome.rfamtr.nr.gff
TRNA_GFF_2=$OUTDIR/input_genome.alltrna.nr.gff

cat $BASEDIR/input_genome.rfam.gff | awk -F'\t' '($3~/tRNA/)' | sort -k1,1 -k4,4n > $OUTDIR/input_genome.rfamtr.sorted.gff
$SCRIPTPATH/bin/bedtools cluster -i $OUTDIR/input_genome.rfamtr.sorted.gff -d 0 -s | sort -k10,10 -k6,6nr | sort -k10,10 -u | sort -k1,1 -k4,4n | cut -f 1-9 > $TRNA_GFF

if [ -e $BASEDIR/input_genome.trna.gff ]
then

cat $BASEDIR/input_genome.trna.gff | sort -k1,1 -k4,4n > $OUTDIR/input_genome.trna.sorted.gff
$SCRIPTPATH/bin/bedtools cluster -i $OUTDIR/input_genome.trna.sorted.gff -d 0 -s | sort -k10,10 -k6,6nr | sort -k10,10 -u | sort -k1,1 -k4,4n | cut -f 1-9 > $OUTDIR/input_genome.trna.nr.gff

cat $OUTDIR/input_genome.trna.nr.gff $OUTDIR/input_genome.rfamtr.nr.gff | sort -k1,1 -k4,4n > $OUTDIR/input_genome.alltrna.sorted.gff
$SCRIPTPATH/bin/bedtools cluster -i $OUTDIR/input_genome.alltrna.sorted.gff -d 0 -s | awk -F'\t' '{print $0, ($2~/tRNAscan-SE/)?1:0}' 'OFS=\t' | sort -k10,10 -k11,11nr | sort -k10,10 -u | sort -k1,1 -k4,4n | cut -f 1-9 > $TRNA_GFF_2

fi

RRNA_GFF=$OUTDIR/input_genome.rfamrr.nr.gff
RRNA_GFF_2=$OUTDIR/input_genome.allrrna.nr.gff

cat $BASEDIR/input_genome.rfam.gff | awk -F'\t' '($3~/rRNA/)' | sort -k1,1 -k4,4n > $OUTDIR/input_genome.rfamrr.sorted.gff
$SCRIPTPATH/bin/bedtools cluster -i $OUTDIR/input_genome.rfamrr.sorted.gff -d 0 -s | sort -k10,10 -k6,6nr | sort -k10,10 -u | sort -k1,1 -k4,4n | cut -f 1-9 > $RRNA_GFF

if [ -e $BASEDIR/input_genome.rnammer.gff ]
then

cat $BASEDIR/input_genome.rnammer.gff | sort -k1,1 -k4,4n > $OUTDIR/input_genome.rnammer.sorted.gff
$SCRIPTPATH/bin/bedtools cluster -i $OUTDIR/input_genome.rnammer.sorted.gff -d 0 -s | sort -k10,10 -k6,6nr | sort -k10,10 -u | sort -k1,1 -k4,4n | cut -f 1-9 > $OUTDIR/input_genome.rnammer.nr.gff

cat $OUTDIR/input_genome.rnammer.nr.gff $OUTDIR/input_genome.rfamrr.nr.gff | sort -k1,1 -k4,4n > $OUTDIR/input_genome.allrrna.sorted.gff
$SCRIPTPATH/bin/bedtools cluster -i $OUTDIR/input_genome.allrrna.sorted.gff -d 0 -s | awk -F'\t' '{print $0, ($2~/RNAmmer/)?1:0}' 'OFS=\t' | sort -k10,10 -k11,11nr | sort -k10,10 -u | sort -k1,1 -k4,4n | cut -f 1-9 > $RRNA_GFF_2

fi

OTHER_GFF=$OUTDIR/input_genome.allother.nr.gff

cat $BASEDIR/input_genome.rfam.gff | awk -F'\t' '(!($3~/tRNA/||$3~/rRNA/))' | sort -k1,1 -k4,4n > $OUTDIR/input_genome.rfam.sorted.gff
cat $BASEDIR/input_genome.cmcore.gff $BASEDIR/input_genome.sno.gff $BASEDIR/input_genome.u3.gff $BASEDIR/input_genome.srp.gff $OUTDIR/input_genome.rfam.sorted.gff | sort -k1,1 -k4,4n > $OUTDIR/input_genome.allother.sorted.gff
$SCRIPTPATH/bin/bedtools cluster -i $OUTDIR/input_genome.allother.sorted.gff -d 0 -s | sort -k10,10 -k6,6nr | sort -k10,10 -u | sort -k1,1 -k4,4n | cut -f 1-9 > $OTHER_GFF


if [ -e $TRNA_GFF_2 ] && [ -e $RRNA_GFF_2 ]
then

cat $TRNA_GFF_2 $RRNA_GFF_2 $OTHER_GFF | sort -k1,1 -k4,4n > $OUTDIR/input_genome.final.nr.gff

elif [ -e $TRNA_GFF_2 ] && [ ! -e $RRNA_GFF_2 ]
then

cat $RRNA_GFF $OTHER_GFF | sort -k1,1 -k4,4n > $OUTDIR/input_genome.final.1.gff
$SCRIPTPATH/bin/bedtools cluster -i $OUTDIR/input_genome.final.1.gff -d 0 -s | sort -k10,10 -k6,6nr | sort -k10,10 -u | sort -k1,1 -k4,4n | cut -f 1-9 > $OUTDIR/input_genome.final.2.gff
cat $TRNA_GFF_2 $OUTDIR/input_genome.final.2.gff | sort -k1,1 -k4,4n > $OUTDIR/input_genome.final.nr.gff
rm -f $OUTDIR/input_genome.final.1.gff
rm -f $OUTDIR/input_genome.final.2.gff

elif [ ! -e $TRNA_GFF_2 ] && [ -e $RRNA_GFF_2 ]
then

cat $TRNA_GFF $OTHER_GFF | sort -k1,1 -k4,4n > $OUTDIR/input_genome.final.1.gff
$SCRIPTPATH/bin/bedtools cluster -i $OUTDIR/input_genome.final.1.gff -d 0 -s | sort -k10,10 -k6,6nr | sort -k10,10 -u | sort -k1,1 -k4,4n | cut -f 1-9 > $OUTDIR/input_genome.final.2.gff
cat $RRNA_GFF_2 $OUTDIR/input_genome.final.2.gff | sort -k1,1 -k4,4n > $OUTDIR/input_genome.final.nr.gff
rm -f $OUTDIR/input_genome.final.1.gff
rm -f $OUTDIR/input_genome.final.2.gff

else

cat $TRNA_GFF $RRNA_GFF $OTHER_GFF | sort -k1,1 -k4,4n > $OUTDIR/input_genome.final.1.gff
$SCRIPTPATH/bin/bedtools cluster -i $OUTDIR/input_genome.final.1.gff -d 0 -s | sort -k10,10 -k6,6nr | sort -k10,10 -u | sort -k1,1 -k4,4n | cut -f 1-9 > $OUTDIR/input_genome.final.2.gff
cat $OUTDIR/input_genome.final.2.gff | sort -k1,1 -k4,4n > $OUTDIR/input_genome.final.nr.gff
rm -f $OUTDIR/input_genome.final.1.gff
rm -f $OUTDIR/input_genome.final.2.gff

fi
