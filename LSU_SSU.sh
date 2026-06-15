INFILE=$1
TMPDIR=$2

mkdir $TMPDIR

cat $INFILE | awk -F'\t' '($3~/SSU_rRNA/||$3~/18s_rRNA/){print $1, $4, $5, $10, $6, $7}' OFS='\t' > $TMPDIR/ssu.bed
cat $INFILE | awk -F'\t' '($3~/LSU_rRNA/||$3~/28s_rRNA/){print $1, $4, $5, $10, $6, $7}' OFS='\t' > $TMPDIR/lsu.bed

bedtools window \
-a $TMPDIR/ssu.bed \
-b $TMPDIR/lsu.bed \
-l 1500 \
-r 1500 \
-sm | perl -F'\t' -lane '{$strand=$F[5];
$start_1=($strand=="+")?$F[1]:$F[2];
$end_1=($strand=="+")?$F[2]:$F[1];
$start_2=($strand=="+")?$F[7]:$F[8];
$end_2=($strand=="+")?$F[8]:$F[7];
$length=($start_2-$end_1>0)?($start_2-$end_1):($start_1-$end_2);
$start=($start_1<$end_2)?$start_1:$end_2;
$end=($start_1<$end_2)?$end_2:$start_1;
print "$F[0]\t$start\t$end\tSSU_LSU:$start-$end\t$length\t$strand"}' | \
sort -k1,1 -k2,2n | awk -F'\t' '($5>=200){print $0}' | \
bedtools cluster -i stdin -d 0 -s | \
sort -k7,7 -u | wc -l

rm -rf $TMPDIR
