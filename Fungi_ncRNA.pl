use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $fna_in = "NA";
my $gff_in = "NA";
my $outdir = "Fungi_ncRNA";
my $fungal_class = "Fungi";
my $run_rnammer_trnascan = 0;
my $gffcmp = 0;
my $assm = "NA";
my $core_cm_cutoff = 60;
my $sno_cm_cutoff = 30;
my $srp_cm_cutoff = 30;
my $exclude_sn_variants = 1;
my $precomp_dir = "NA";
my $recomp = 0;
my $merge_separately = 0;
my $rfam_path = "$cd_path/Rfam/Rfam_f.cm";

my $ind = 0;
foreach(@ARGV) {
	
	#Required!! Input genome (fasta format)
	if (@ARGV[$ind] eq '-fna') {
		$fna_in = @ARGV[$ind + 1];
		if (! (-e $fna_in)) {
			die "cannot open genomic file: $fna_in\n";
		}
	}
	
	#Optional, but recommended. Annotation file in GFF format
	if (@ARGV[$ind] eq '-gff') {
		$gff_in = @ARGV[$ind + 1];
	}
	
	#Optional, taxomonical class of fungus. Class-specific covariation models will be used. Otherwise, it will be one fine-tuned for all fungi (class unspecific).
	if (@ARGV[$ind] eq '-fungal_class') {
		$fungal_class = @ARGV[$ind + 1];
	}
	
	#Run optional dependencies RNAMMER and tRNAScan-SE. Not run by default. 
	if (@ARGV[$ind] eq '-run_rnammer_trnascan') {
		$run_rnammer_trnascan = 1;
	}
	
	#Run GFFCompare to test whether predict ncRNAs overlap with features on the same strand. Not run by default. Skipped if a GFF annotation is not provided.
	if (@ARGV[$ind] eq '-gffcmp') {
		$gffcmp = 1;
	}
	
	#Optional, the output directory path. Fungi_ncRNA if unspecified.
	if (@ARGV[$ind] eq '-out') {
		$outdir = @ARGV[$ind + 1];
	}
	
	#Optional, assembly ID for multi-contig genomic files. NA if unspecified.
	if (@ARGV[$ind] eq '-gcf') {
		$assm = @ARGV[$ind + 1];
	}
	
	#CMSearch cutoff used for fungi fine-tuned or class-specific core ncRNA covariation models. Default 60.
	if (@ARGV[$ind] eq '-core_cm_cutoff') {
		$core_cm_cutoff = @ARGV[$ind + 1];
	}
	
	#CMSearch cutoff used for fungi fine-tuned snoRNA covariation models. Default 30.
	if (@ARGV[$ind] eq '-sno_cm_cutoff') {
		$sno_cm_cutoff = @ARGV[$ind + 1];
	}
	
	#CMSearch cutoff used for fungi SRP covariation models. Default 30.
	if (@ARGV[$ind] eq '-srp_cm_cutoff') {
		$srp_cm_cutoff = @ARGV[$ind + 1];
	}
	
	#Also predict variants of snRNAs. Default no.
	if (@ARGV[$ind] eq '-include_sn_variants') {
		$exclude_sn_variants = 0;
	}
	
	#Reanalyse an output directory, running RNAMMER and tRNAScan-SE for additional tRNA and tRNA if not ran previously.
	if (@ARGV[$ind] eq '-reanalyse') {
		$precomp_dir = @ARGV[$ind + 1];
		$recomp = 1;
	}
	
	#When overlapping occurs, take the best LSU, SSU or 5S rRNA predicted by RNAMMER if any, otherwise the best LSU, SSU and 5S predicted by CMSearch using the Rfam covariation model. Similarly, take the best tRNA predicted by tRNAScan-SE if any, otherwise the best tRNA predicted by CMSearch using the Rfam covariation model. For the remaining, take the one with the best score. Default behaviour is to always take the one with the best score when overlapping occurs. 
	if (@ARGV[$ind] eq '-merge_separately') {
		$merge_separately = 1;
	}
	
	#Specify a different path to the Rfam covariation models, useful when Rfam is updated. Default is Rfam/Rfam_f.cm (a subset containing models with at least one hit over the RefSeq v219 fungal genomes).
	if (@ARGV[$ind] eq '-rfam_path') {
		$rfam_path = @ARGV[$ind + 1];
	}
	
	#print help.
	if (@ARGV[$ind] eq '-h') {
		system("cat $cd_path/help.txt");
	}
	
	$ind++;
}

if (! -e $rfam_path) {
	$rfam_path = "$cd_path/Rfam/Rfam_f.cm";
}

my %core = ();
if (-e "$cd_path/data/core_list.txt") {
	open(CORE_CM, "$cd_path/data/core_list.txt");
	while(my $line = <CORE_CM>) {
		chomp $line;
	
		my $mdl = $line;
		if (not exists $core{$mdl}) {
			$core{$mdl} = 1;
		}
	}
	close(CORE_CM);
}

my %sno = ();
if (-e "$cd_path/data/sno_list.txt") {
open(SNO_CM, "$cd_path/data/sno_list.txt");
	while(my $line = <SNO_CM>) {
		chomp $line;
	
		my $mdl = $line;
		if (not exists $sno{$mdl}) {
			$sno{$mdl} = 1;
		}
	}
	close(SNO_CM);
}

if (-d $outdir) {
	system("rm -rf $outdir");
}
mkdir($outdir);

if ((-d "$precomp_dir/NCRNA") and ($recomp > 0)) {
	system("sh $cd_path/Fungi_ncRNA_recomp.sh $precomp_dir/NCRNA $outdir/NCRNA $run_rnammer_trnascan $core_cm_cutoff $sno_cm_cutoff $srp_cm_cutoff $exclude_sn_variants >> $outdir/log.txt 2>> $outdir/log.err.txt");
} else{
	system("sh $cd_path/Fungi_ncRNA_search.sh $fna_in $outdir/NCRNA $fungal_class $run_rnammer_trnascan $core_cm_cutoff $sno_cm_cutoff $srp_cm_cutoff $exclude_sn_variants $rfam_path >> $outdir/log.txt 2>> $outdir/log.err.txt");
}

if (($merge_separately > 0) and ((-e "$outdir/NCRNA/input_genome.rnammer.gff") or (-e "$outdir/NCRNA/input_genome.trna.gff"))) {
	system("sh $cd_path/Fungi_ncRNA_sep_merge.sh $outdir/NCRNA $outdir/NCRNA/merge");
	system("rm -f $outdir/NCRNA/input_genome.all.nr.gff");
	system("cp $outdir/NCRNA/merge/input_genome.final.nr.gff $outdir/NCRNA/input_genome.all.nr.gff");
	system("rm -rf $outdir/NCRNA/merge");
}

my $new_fun_gff = "$outdir/input_genome.ncrna.gff";
my $new_fun_gff_assm = "$outdir/NCRNA/input_genome.all.nr.assm.gff";

my $count = 0;
open(NEWGFF, ">$new_fun_gff");
open(ASSMGFF, ">$new_fun_gff_assm");
open(FUNGFF, "$outdir/NCRNA/input_genome.all.nr.gff");
while(my $line = <FUNGFF>) {
	chomp $line;
	
	my @toks = split(/[\t]/, $line);
	my $contig = $toks[0];
	my $app = $toks[1];
	my $type = $toks[2]; 
	$count++;
	
	my $id = "$contig.$type.$count";
	my $comment_1 = "ID=$id;Name=$id;Dbxref=$type";
	my $comment_2 = "ID=$id.1;Name=$id.1;Parent=$id";
	my $comment_3 = "ID=$id.1.exon.1;Parent=$id.1";
	
	$toks[2] = "gene"; $toks[8] = $comment_1; print NEWGFF join("\t",@toks) . "\n";
	$toks[2] = "mRNA"; $toks[8] = $comment_2; print NEWGFF join("\t",@toks) . "\n";
	$toks[2] = "exon"; $toks[8] = $comment_3; print NEWGFF join("\t",@toks) . "\n";
	
	my $query = $type;
	
	if (index($id, "basidiomycota.U3.precursor") >= 0) {
		$query = "basidiomycota_U3_precursor";
	}
	
	$query =~ s/_Eurotiomycetes$//g;
    $query =~ s/_Sordariomycetes$//g;
	$query =~ s/_Saccharomycetes$//g;
	$query =~ s/_Dothideomycetes$//g;
	$query =~ s/_Agaricomycetes$//g;
	$query =~ s/_Leotiomycetes$//g;
	$query =~ s/_Tremellomycetes$//g;
	$query =~ s/_Microsporidia$//g;
	$query =~ s/_Fungi$//g;
	
	my $class = "Other";
	
	if(index($query,"rRNA") >= 0){
		$class = "rRNA";
	}
	
	if(index($query,"tRNA") >= 0){
		$class = "tRNA";
	}
	
	if (exists $sno{$query} or (index($query, "basidiomycota_U3_precursor") >= 0)) {
		$class = "snoRNA";
	}
	
	if ((index($app, "tRNAscan-SE") >= 0) and (index($query, "pseudogene") >= 0)) {
		$class = "tRNA";
	}
	
	if (index($query, "SRP_Basido") >= 0) {
		$class = "Fungi_SRP";
	}
	
	if (exists $core{$query}) {
		$class  = $query;
	}
	
	print ASSMGFF "$line\t$id\t$assm\t$class\n";
	
}
close(FUNGFF);
close(ASSMGFF);
close(NEWGFF);

if (($gffcmp > 0) and (-e $gff_in)) {
	system("perl $cd_path/GFFCompare_wrapper.pl $new_fun_gff $gff_in $fna_in $outdir/GFFCMP $assm 0");
}

if (-e $gff_in) {
	system("cp $gff_in $outdir/input_genome.gff");
	system(qq~cat $gff_in $new_fun_gff | grep -v '^#' > $outdir/input_genome.combined.gff.tmp~);
	system(qq~sort -k1,1 -k4,4n $outdir/input_genome.combined.gff.tmp > $outdir/input_genome.combined.gff~);
	unlink("$outdir/input_genome.combined.gff.tmp");
}

my %gffcmp_code = ();
my %gffcmp_interp = ();

my %code_interp = ();
$code_interp{"="} = "Reference match";
$code_interp{"c"} = "Reference overlap";
$code_interp{"e"} = "Reference overlap";
$code_interp{"i"} = "Intronic";
$code_interp{"k"} = "Reference overlap";
$code_interp{"m"} = "Reference overlap";
$code_interp{"n"} = "Reference overlap";
$code_interp{"o"} = "Reference overlap";
$code_interp{"p"} = "Novel loci";
$code_interp{"r"} = "Repeat overlap";
$code_interp{"u"} = "Novel loci";
$code_interp{"x"} = "Antisense to reference";

if (-e "$outdir/GFFCMP/input_genome.fun.gffcmp.annotated.summ.txt") {
	open(GFFCMP_SUM, "$outdir/GFFCMP/input_genome.fun.gffcmp.annotated.summ.txt");
	while(my $line = <GFFCMP_SUM>) {
		chomp $line;
	
		my @toks = split(/[\t]/, $line);
		my $id = $toks[8];
		my $loc = $toks[5];
	
		my $interp = "Novel loci";
		if (exists $code_interp{$loc}) {
			$interp = $code_interp{$loc};
		}
	
		if (not exists $gffcmp_code{$id}) {
			$gffcmp_code{$id} = $loc;
		}
	
		if (not exists $gffcmp_interp{$id}) {
			$gffcmp_interp{$id} = $interp;
		}
	}
	close(GFFCMP_SUM);
}

my $new_fun_gff_assm_gffcmp = "$outdir/NCRNA/input_genome.all.nr.assm.gffcmp.gff";

open(R_ASSMGFF, $new_fun_gff_assm);
open(ASSMGFFCMPGFF, ">$new_fun_gff_assm_gffcmp");
while(my $line = <R_ASSMGFF>) {
	chomp $line;
	
	my @toks = split(/[\t]/, $line);
	my $id = $toks[9];
	
	my $loc = "u";
	if (exists $gffcmp_code{$id}) {
		$loc = $gffcmp_code{$id};
	}
	
	my $interp = "Novel loci";
	if (exists $gffcmp_interp{$id}) {
		$interp = $gffcmp_interp{$id};
	}
	
	print ASSMGFFCMPGFF "$line\t$loc\t$interp\n";
}
close(ASSMGFFCMPGFF);
close(R_ASSMGFF);

system("perl $cd_path/PredictionSummary.pl $new_fun_gff_assm_gffcmp $outdir/tmp > $outdir/summary.txt 2> $outdir/log.err.txt");

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}
