use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $fun_gff_in = @ARGV[0]; #"out/NCRNA/input_genome.all.nr.gff";
my $ref_gff_in = @ARGV[1]; #"/Volumes/archive/brownlab/DB/RefSeq219/fungi/C/Coprinopsis_cinerea_okayama7#130-240176#GCF_000182895.1/GCF_000182895.1_CC3_genomic.gff";
my $ref_fna_in = @ARGV[2]; #"/Volumes/archive/brownlab/DB/RefSeq219/fungi/C/Coprinopsis_cinerea_okayama7#130-240176#GCF_000182895.1/GCF_000182895.1_CC3_genomic.fna";
my $outdir = @ARGV[3]; #"gffcmp";
my $assm = @ARGV[4]; #"GCF_000182895.1";
my $reformat_i2gff = @ARGV[5];

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

system("cp $fun_gff_in $outdir/input_genome.fun.gff");
system("cp $ref_fna_in $outdir/input_genome.fna");
system("cp $ref_gff_in $outdir/input_genome.gff");

my $new_fun_gff = "$outdir/input_genome.all.nr.gffcmp.gff";

if ($reformat_i2gff < 1) {
	system("cp $fun_gff_in $new_fun_gff");
} else {
	my $count = 0;
	open(NEWGFF, ">$new_fun_gff");
	open(FUNGFF, "$outdir/input_genome.fun.gff");
	while(my $line = <FUNGFF>) {
		chomp $line;
	
		my @toks = split(/[\t]/, $line);
		my $contig = $toks[0];
		my $type = $toks[2]; 
		$count++;
	
		my $id = "$contig.$type.$count";
		my $comment_1 = "ID=$id;Name=$id;Dbxref=$type";
		my $comment_2 = "ID=$id.1;Name=$id.1;Parent=$id";
		my $comment_3 = "ID=$id.1.exon.1;Parent=$id.1";
	
		$toks[2] = "gene"; $toks[8] = $comment_1; print NEWGFF join("\t",@toks) . "\n";
		$toks[2] = "mRNA"; $toks[8] = $comment_2; print NEWGFF join("\t",@toks) . "\n";
		$toks[2] = "exon"; $toks[8] = $comment_3; print NEWGFF join("\t",@toks) . "\n";
	}
	close(FUNGFF);
	close(NEWGFF);
}

my $new_fun_gtf = "$outdir/input_genome.fun.gffcmp.gtf";
system("gffread $new_fun_gff -T -o $new_fun_gtf");
system("gffcompare --debug -d 10 -e 10 -o $outdir/input_genome.fun.gffcmp.out -r $outdir/input_genome.gff -s $outdir/input_genome.fna $new_fun_gtf >> $outdir/log.txt 2>> $outdir/log.err.txt");
#system("$cd_path/bin/gffread $new_fun_gff -T -o $new_fun_gtf");
#system("$cd_path/bin/gffcompare --debug -d 10 -e 10 -o $outdir/input_genome.fun.gffcmp.out -r $outdir/input_genome.gff -s $outdir/input_genome.fna $new_fun_gtf >> $outdir/log.txt 2>> $outdir/log.err.txt");

my $final_gffcmp_gtf = "$outdir/input_genome.fun.gffcmp.annotated.gtf";
my $summ = "$outdir/input_genome.fun.gffcmp.annotated.summ";

open(SUMM, ">$summ.txt");
open(FGTF, $final_gffcmp_gtf);
while(my $line = <FGTF>) {
	chomp $line;
	
	my @toks = split(/[\t]/, $line);
	
	my $app = $toks[1];
	my $feat = $toks[2];
	my $comment = $toks[8];
	
	if ($feat ne "transcript") {
		next;
	}
	#print "$line\n";
	
	my ($gene_id) = ($comment =~ /gene_id \"([^"]+)\";/);
	my ($class_code) = ($comment =~ /class_code \"([^"]+)\";/);
	my ($contig,$query,$id_short) = ($gene_id =~ /(\S+)\.(\S+)\.(\S+)/);
	
	if (index($gene_id, "basidiomycota.U3.precursor") >= 0) {
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
	
	my $super_class = "Other";
	
	if(index($query,"rRNA") >= 0){
		$super_class = "rRNA";
	}
	
	if(index($query,"tRNA") >= 0){
		$super_class = "tRNA";
	}
	
	if (exists $sno{$query} or (index($query, "basidiomycota_U3_precursor") >= 0)) {
		$super_class = "snoRNA";
	}
	
	if ((index($app, "tRNAScan-SE") >= 0) and (index($query, "pseudogene") >= 0)) {
		$super_class = "tRNA";
	}
	
	if (index($query, "SRP_Basido") >= 0) {
		$super_class = "Fungi_SRP";
	}
	
	if (exists $core{$query}) {
		$super_class  = $query;
	}
	
	my $comment_new = $gene_id;
	$toks[2] = $query;
	$toks[5] = $class_code;
	$toks[7] = $assm;
	$toks[8] = $comment_new;
	
	my $strout = join("\t", @toks) . "\t$super_class";
	print SUMM "$strout\n";
	
}
close(FGTF);
close(SUMM);

unlink("$outdir/input_genome.fun.gff");
unlink("$outdir/input_genome.gff");
unlink("$outdir/input_genome.fna");
unlink("$outdir/input_genome.fna.fai");

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}

#cat PV_rfam.txt | perl -F'\t' -lane 'BEGIN{$count=0;}{$type=$F[2];$count++;$contig=$F[0];$id="$contig.Rfam.$count.$type";$comment_1="ID=$id;Name=$id;Dbxref=$type";$comment_2="ID=$id.1;Name=$id.1;Parent=$id";$comment_3="ID=$id.1.exon.1;Parent=$id.1";$F[2]="gene";$F[8]=$comment_1;print join("\t",@F);$F[2]="mRNA";$F[8]=$comment_2;print join("\t",@F);$F[2]="exon";$F[8]=$comment_3;print join("\t",@F);}' > PV_rfam.gff3
#$path_gff=$path;$path_gff=~s/\.fna$/.gff/g;

# cat $OUTDIR_GFFCMP/input_genome.allother.annotated.gtf | awk -F'\t' '($3~/transcript/&&$9~/class_code/)' | \
# perl -F'\t' -lane '{$comment=$F[8]; ($transcript_id,$class_code)=($comment=~/transcript_id\s\"(\S+)\"\;\s.*\sclass_code\s\"(\S+)\"\;\s*.*/);
#                     $name=$transcript_id; $name=~s/\_\d+$//g; $query=$name; 
#                     if(index($query,"U3")>=0){$query="U3";}
# 					$query=~s/_Eurotiomycetes$//g;
# 					$query=~s/_Sordariomycetes$//g;
# 					$query=~s/_Saccharomycetes$//g;
# 					$query=~s/_Dothideomycetes$//g;
# 					$query=~s/_Agaricomycetes$//g;
# 					$query=~s/_Other_Basidiomycota$//g;
# 					$query=~s/_Leotiomycetes$//g;
# 					$query=~s/_Tremellomycetes$//g;
# 					$query=~s/_Other_Ascomycota$//g;
# 					$query=~s/_Other_Microsporidia$//g;
# 					$query=~s/_Other_Fungi$//g;
# 					$query=~s/_Fungi$//g;
# 					$comment_new=$transcript_id;
# 					$F[2]=$query;$F[5]=$class_code;$F[8]=$comment_new;$F[7]='$ASSM';
# 					$strout=join("\t",@F);print $strout;
#                     }' > $OUTDIR/input_genome.allother.annotated.summ.gff
# 
# cat $OUTDIR_GFFCMP/input_genome.allrtrna.annotated.gtf | awk -F'\t' '($3~/transcript/&&$9~/class_code/)' | \
# perl -F'\t' -lane '{$comment=$F[8]; ($transcript_id,$class_code)=($comment=~/transcript_id\s\"(\S+)\"\;\s.*\sclass_code\s\"(\S+)\"\;\s*.*/);
#                     $name=$transcript_id; $name=~s/\_\d+$//g; $query="tRNA"; 
#                     if(index($name,"rRNA")>=0){$query="rRNA";}
# 					$comment_new=$transcript_id;
# 					$F[2]=$query;$F[5]=$class_code;$F[8]=$comment_new;$F[7]='$ASSM';
# 					$strout=join("\t",@F);print $strout;
#                     }' > $OUTDIR/input_genome.allrtrna.annotated.summ.gff
