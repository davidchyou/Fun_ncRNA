my $gff = @ARGV[0]; #"AFA1163_2/NCRNA/input_genome.all.nr.assm.gffcmp.gff";

my $assembly = "";
my @ncrna_types = ("U1", "U2", "U4", "U5", "U6", "U1_yeast", "RNaseP_nuc", "RNase_MRP", "Fungi_SRP", "snoRNA", "tRNA", "LSU_rRNA", "SSU_rRNA", "5S_rRNA", "5.8S_rRNA");

my %counts = ();
foreach my $t (@ncrna_types) {
	$counts{$t} = 0;
}

open(GFF, $gff);
while(my $line = <GFF>) {
	chomp $line;
	
	my @toks = split(/[\t]/, $line);
	my $ncrna_name = $toks[2];
	my $ncrna_type = $toks[11];
	my $assm = $toks[10];
	
	if ($assembly eq "") {
		$assembly = $assm;
	}
	
	if ($ncrna_type eq "rRNA") {
		if (($ncrna_name eq "28s_rRNA") or ($ncrna_name eq "LSU_rRNA_eukarya") or ($ncrna_name eq "LSU_rRNA_bacteria") or ($ncrna_name eq "LSU_rRNA_archaea")) {
			$counts{"LSU_rRNA"} += 1;
		} elsif (($ncrna_name eq "18s_rRNA") or ($ncrna_name eq "SSU_rRNA_eukarya") or ($ncrna_name eq "SSU_rRNA_bacteria") or ($ncrna_name eq "SSU_rRNA_archaea") or ($ncrna_name eq "SSU_rRNA_microsporidia")) {
			$counts{"SSU_rRNA"} += 1;
		} elsif (($ncrna_name eq "5S_rRNA") or ($ncrna_name eq "5s_rRNA") or ($ncrna_name eq "8s_rRNA")) {
			$counts{"5S_rRNA"} += 1;
		} elsif ($ncrna_name eq "5_8S_rRNA") {
			$counts{"5.8S_rRNA"} += 1;
		} else {
			next;
		}
	} else {
		$counts{$ncrna_type} += 1;
	}
}

my $n_snrna = $counts{"U1"} + $counts{"U2"} + $counts{"U4"} + $counts{"U5"} + $counts{"U6"} + $counts{"U1_yeast"} + 0.0;
my $n_snorna = $counts{"snoRNA"} + 0.0;
my $n_rnasep = $counts{"RNaseP_nuc"} + 0.0;
my $n_rnasemrp = $counts{"RNase_MRP"} + 0.0;
my $n_rnasesrp = $counts{"Fungi_SRP"} + 0.0;
my $n_trna = $counts{"tRNA"} + 0.0;
my $n_lsu = $counts{"LSU_rRNA"} + 0.0;
my $n_ssu = $counts{"SSU_rRNA"} + 0.0;
my $n_5s = $counts{"5S_rRNA"} + 0.0;
my $n_58s = $counts{"5.8S_rRNA"} + 0.0;
my $n_rrna = $n_lsu + $n_ssu + $n_5s + $n_58s;
my $n_repeat = ($n_lsu + $n_ssu + $n_58s) / 3;

print "# Core ncRNAs were predicted with Fun_ncRNA (Chyou et al 2026). There were $n_snrna snRNA, $n_snorna snoRNA, $n_rnasep RNase P, $n_rnasemrp RNase MRP, $n_rnasesrp SRP, $n_trna tRNA and $n_rrna rRNA predicted. Of the $n_rrna rRNA, there were $n_5s 5S rRNA, $n_lsu LSU, $n_58s 5.8S rRNA, and $n_ssu SSU predicted. The number of complete or partial rRNA repeat units of SSU, 5.8S and LSU in this assembly is estimated to be $n_repeat.\n";

foreach my $t (@ncrna_types) {
	my $count = $counts{$t};
	print "$t\t$count\t$assembly\n"; 
}
