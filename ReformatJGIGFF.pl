my $infile = @ARGV[0]; #"Tatroviridev2_FrozenGeneCatalog_20100319.gff";
my $outfile = @ARGV[1];

my %all_tids = ();
my %contigs = ();
my %f_start = ();
my %f_end = ();
my %strands = ();
my %names = ();

open(GFF, $infile);
while(my $line = <GFF>) {
	chomp $line;
	
	my @toks = split(/[\t]/, $line);
	my $contig = $toks[0];
	my $type = $toks[2];
	my $start = $toks[3];
	my $end = $toks[4];
	my $strand = $toks[6];
	my $comment = $toks[8];
	
	my ($tid) = ($comment =~ /transcriptId\s([^;]+)/);
	my ($pid) = ($comment =~ /proteinId\s([^;]+)/);
	my ($name) = ($comment =~ /name\s([^;]+)/);
	
	if ($type eq "CDS") {
		$tid = $pid;
		if (not exists $all_tids{$tid}) {
			$all_tids{$tid} = 1;
			$contigs{$tid} = $contig;
			$f_start{$tid} = $start;
			$f_end{$tid} = $end;
			$strands{$tid} = $strand;
			$names{$tid} = $name;
		} else {
			if ($f_start{$tid} > $start) {
				$f_start{$tid} = $start;
			}
		
			if ($f_end{$tid} < $end) {
				$f_end{$tid} = $end;
			}
		}
	}
}
close(GFF);

open(OUT, ">$outfile");
my %processed = ();
open(GFF, $infile);
while(my $line = <GFF>) {
	chomp $line;
	
	my @toks = split(/[\t]/, $line);
	my $contig = $toks[0];
	my $type = $toks[2];
	my $start = $toks[3];
	my $end = $toks[4];
	my $strand = $toks[6];
	my $comment = $toks[8];
	
	my ($tid) = ($comment =~ /transcriptId\s([^;]+)/);
	my ($pid) = ($comment =~ /proteinId\s([^;]+)/);
	my ($name) = ($comment =~ /name\s([^;]+)/);
	
	if ($type eq "CDS") {
		$tid = $pid;
		
		if (not exists $processed{$tid}) {
			print OUT "$contigs{$tid}\tJGI\tgene\t$f_start{$tid}\t$f_end{$tid}\t.\t$strands{$tid}\t.\tID=gene-$tid;Name=$names{$tid};locus_tag=$tid;Dbxref=$tid\n";
			print OUT "$contigs{$tid}\tJGI\tmRNA\t$f_start{$tid}\t$f_end{$tid}\t.\t$strands{$tid}\t.\tID=rna-$tid;Name=$names{$tid};Parent=gene-$tid;locus_tag=$tid;Dbxref=$tid\n";
			$processed{$tid} = 1;
		}
		
		my ($exon_number) = ($comment =~ /exonNumber\s([^;]+)/);
		my $comment_new = "ID=exon-$tid-$exon_number;Name=$names{$tid};Parent=rna-$tid;locus_tag=$tid;Dbxref=$tid";
		
		$toks[2] = "exon";
		$toks[8] = $comment_new;
		print OUT join("\t", @toks) . "\n";
		
		$comment_new = "ID=cds-$tid;Name=$names{$tid};Parent=rna-$tid;locus_tag=$tid;Dbxref=$tid";
		$toks[2] = "CDS";
		$toks[8] = $comment_new;
		print OUT join("\t", @toks) . "\n";
	} else {
		next;
	}
}
close(GFF);
close(OUT);
