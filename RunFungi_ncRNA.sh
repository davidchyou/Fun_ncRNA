fna_path="./data/GCA_000171015.2_TRIAT_v2.0_genomic.fna"
gff_path="./data/GCA_000171015.2_TRIAT_v2.0_genomic.gff"
class="Sordariomycetes"
outdir="TAIMI_ncRNA"

rm -rf $outdir

perl ./Fungi_ncRNA.pl \
-fna $fna_path \
-gff $gff_path \
-fungal_class $class \
-out $outdir \
-gffcmp
