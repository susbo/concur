#!/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;

my $version = "0.1";

GetOptions("version" => sub { VersionMessage() }
          , "help" => \my $help
          , "man" => \my $man
          , 'input=s' => \my $input
          , 'genome=s' => \my $genome
          , 'out=s' => \my $out
          , 'name=s' => \my $name
          , 'threads=s' => \my $threads
          ) or pod2usage(2);

pod2usage(2) if $help or $man;
unless ($input && $genome && $out) { pod2usage("Missing mandatory arguments."); pod2usage(2); }
unless ($threads) { $threads=1; }

if (!$name) { # Set name to input file name by default.
	my @tmp = split "/",$input;
	$name = $tmp[$#tmp]; $name =~ s/.bam//;
}

#### Info #####
# Print welcome message and information.

print "Running CONCUR v$version\n\n";
print "### Parameters ###\n";
print "Input file: $input\n";
print "Genome: $genome\n";
print "Output folder: $out\n";
print "Output file name: $name\n";
print "Number of threads: $threads\n";
print "##################\n";

#### Running CONCUR ####
my $tmp = "$out/tmp";
mkdir $out;
mkdir $tmp;

# Run sub-routins
#GenomeToTranscript();
#PeriodicityAtTSS();
PredictFrame();
CodonFrequency();
CodonFrequency2();
PlotCodonFrequency("codon_frequency/$name","correlation");
CorrelateToBest();
PlotCodonFrequency("codon_frequency/$name.selected","correlation.selected");
#FinalCodonFrequency();
#Validate();

sub GenomeToTranscript {
	print "Mapping genomic reads to transcripts...\n";

	# Convert alignment bam to bed format.
	mkdir "$tmp/genome_to_transcript";
	open IN,"samtools view $input |";
	open OUT,"| gzip -c > $tmp/genome_to_transcript/$name.bed.gz";
	while (my $row = <IN>) {
		chomp $row;
		my @line = split "\t",$row;
		next if $line[5] =~ /[NID]/; # TODO Handle reads with insertions and deletions
		my ($pos,$chr) = ($line[3]-1,$line[2]);
		next if $chr eq "chrM" or $chr eq "MT";
		my $length = length($line[9]);
		print OUT "$chr\t$pos\t";
		print OUT "".($pos+$length)."\t";
		print OUT "$line[0]\t$line[9]\n";
	}
	close IN;
	close OUT;
	# Extract reads with at least 50% overlap to exons.
	SystemBash("bedtools intersect -a <(zcat $tmp/genome_to_transcript/$name.bed.gz) -b data/$genome.bed -wo | gzip -c > $tmp/genome_to_transcript/$name.intersect.gz");
}

# Calculate periodicity at TSS per read length
sub PeriodicityAtTSS {
	print "Calculating periodicity...\n";
	mkdir "$tmp/periodicity";
	my $file = "$tmp/genome_to_transcript/$name.intersect.gz";
	foreach my $len (20..50) {
		SystemBash("zcat $file | awk -v len=\"$len\" '\$11==\"+\" && length(\$5)==len' | awk '{print \$10+\$2-\$7+12}' | awk '\$1<100' | sort | uniq -c | sort -k2,2 -n > $tmp/periodicity/$name.$len.txt");
	}
}

# Predict reading frame for each read length
sub PredictFrame {
	print "Predicting frame per read length...\n";
	open OUT,">$out/$name.predicted_frame.txt";
	my $max_ratio = 0;
	my $max_i;
	my $max_N = 0;

	print OUT "Len\tFrame\tShift\tScoreS\tReads\tFrameW\tScoreW\tCountsW\tFrameT\tScoreT\tCountsT\tScore\tSelected\n";

	foreach my $i (20..50) {
		my @counts;	# Counts per position relative to the TIS
		open IN,"$tmp/periodicity/$name.$i.txt";
		while (my $row = <IN>) {
			$row =~ /\s(\d+) ([0-9-]+)$/;
			next if $2 < -20;
			$counts[$2+20] += $1;
		}
		close IN;
		my @total = (0,0,0);	# Counts per reading frame in TIS region
		my @window = (0,0,0);	# Most abundant frame per three nucleotides
		foreach my $j (-20..99) {
			$total[$j%3] += $counts[$j+20]//0;
			my $b = $j;
			$b = $j+1 if (($counts[$b+20]//0) < ($counts[$j+1+20]//0));
			$b = $j+2 if (($counts[$b+20]//0) < ($counts[$j+2+20]//0));
			$window[$b%3]++ unless (($counts[$b+20]//0) < 4);
		}
		print OUT "$i\t"; # Fragment length

		my $frame_window = IndexOfMax(@window); # Best frame according to sliding window max
		my $frame_total = IndexOfMax(@total); # Best frame according to total counts

		my $N = 0; $N += $_ for @total;	# Total number of reads

		my $predicted_frame = "-"; my $predicted_shift = "-";

		my @shift_score = (0)x3;
		my @shift_diff = ("-")x3;
		my @reads = (0)x40;
		my @predicted_shift = ("-","-","-");
		# Make a frame prediction if sliding window and total agree and total count is greater than threshold.
#		if ($frame_window == $frame_total && $N>10000) {
		if ($N>1000) {
#				$frame = $frame_total;
#				foreach my $j (-17..20) {	# Find predicted shift as first local maxima position at predicted frame
#			foreach my $j (-17..9) {	# Find predicted shift as first local maxima position at predicted frame
			foreach my $j (-6..6) {	# Find predicted shift as first local maxima position at predicted frame
				my $frame = $j%3;
#					print "$i\t$j\tif ((".($counts[$j+20]//0).") > (".($counts[$j+20-3]//0).") && (".($counts[$j+20]//0).") > (".($counts[$j+20+3]//0)."))\t";
#$shift_score = Round(($counts[$j+20]//0) / (($counts[$j+20-3]//0) + ($counts[$j+20]//0) + ($counts[$j+20+3]//0) + 1),3);
#my $shift_score2 = ($counts[$j+20]//0) - ($counts[$j+20-3]//0)/2 - ($counts[$j+20+3]//0)/2;
#				$reads[$j+20] = $shift_score2;
#print "[$shift_score | $shift_score2]\t";
				if (($counts[$j+20]//0) > ($counts[$j+20-3]//0) && ($counts[$j+20]//0) > ($counts[$j+20+3]//0)) { # Only consider local maxima
#print "*";
					my $score = Round(($counts[$j+20]//0) / (($counts[$j+20-3]//0) + ($counts[$j+20]//0) + ($counts[$j+20+3]//0)+1),3);
#					my $diff = ($counts[$j+20]//0) - ($counts[$j+20-3]//0)/2 - ($counts[$j+20+3]//0)/2;
					my $diff = Round(($counts[$j+20]//0)*2 - ($counts[$j+20-9]//0)/9 - ($counts[$j+20-6]//0)/3 - ($counts[$j+20-3]//0) - ($counts[$j+20+3]//0)/3 - ($counts[$j+20+6]//0)/9 - ($counts[$j+20+9]//0)/9 ,0);
#					my $diff = ($counts[$j+20]//0) - ($counts[$j+20-3]//0) - ($counts[$j+30+3]//0);
#					if ($score > $shift_score[$frame]) { # Only replace it better!
					if ($diff ge $shift_diff[$frame] || $shift_diff[$frame] eq "-") { # Only replace it better or not set yet.
						$shift_score[$frame] = $score;
						$shift_diff[$frame] = $diff;
						$predicted_shift[$frame] = $j;
#						last;
					}
				}
#print "\n";
			}
		}
#		print "All Reads: @reads\n";
#		$predicted_shift = IndexOfMax(@reads)-20;
#		print "Best shift: $predicted_shift\n";
#		$shift_score = Round(($counts[$predicted_shift+20]//0) / (($counts[$predicted_shift+20-3]//0) + ($counts[$predicted_shift+20]//0) + ($counts[$predicted_shift+20+3]//0)+1),3);
#		$predicted_shift = "-" if $predicted_shift == "-20";
#print "\n";
		if ($N == 0) { $frame_window = "-"; $frame_total = "-"; }
		print OUT "$predicted_frame\t",join(",",@predicted_shift),"\t",join(",",@shift_score),"\t",join(",",@shift_diff)."\t";
		print OUT "$N\t";

		my $score_window = MaxPercent(@window);
		print OUT "$frame_window\t$score_window\t";
		print OUT join(",",@window)."\t";

#		my $score_total = MaxPercent(@total);
#		print OUT "$frame_total\t$score_total\t";
		my @score_total = (Round($total[0]/($N+1),3),Round($total[1]/($N+1),3),Round($total[2]/($N+1),3));
		print OUT "$frame_total\t",join(",",@score_total)."\t";
		print OUT join(",",@total),"\t";

		my $total_score = Round(1*$score_window*$score_total[0],3);
		print OUT "$total_score\t";

		my $total_score2 = Round(1*$score_window*$score_total[0]*$shift_score[IndexOfMax(@shift_score)],3);
		print OUT "$total_score2\t";

#		my $score = 0;
#		$score += 1 if $score_window > 0.6 && $score_total > 0.4 && $N>=5000 && $shift_score[IndexOfMax(@shift_score)] > 0.5;
		my @score = (0,0,0);
		foreach my $f (0..2) {
#			$score[$f] = 1 if $shift_score[$f] > 0.5 && $total[$f]>1000 && $score_total[$f]>0.15;
			$score[$f] = 1 if $shift_score[$f] > 0.5 && $total[$f]>1000 && $score_total[$f]>0;
		}
		
		print OUT "",join(",",@score)."\n";
	}
	close OUT;
}

# Calculate codon frequencies
sub CodonFrequency {
	print "Calculating codon frequency per read length...\n";
	# Retrieve codon frequencies
	mkdir "$tmp/codon_frequency";

	my @all_shifts;
	open IN,"$out/$name.predicted_frame.txt";
	<IN>; # Skip header line
	while (my $row = <IN>) {
		chomp $row;
		my @line = split "\t",$row;
		my @scores = split ",",$line[$#line];
		my @shift = split ",",$line[2];
print "$line[0]\t";
		foreach my $i (0..2) {
			$shift[$i] = "-" if $scores[$i] == 0;
			$all_shifts[$line[0]][$i] = $shift[$i];
print "$shift[$i]\t";
		}
print "\n";
#		$line[2] = "-" if $line[$#line] == 0;
#		$all_shifts[$line[0]] = $line[2];
	}
	close IN;

	open IN,"zcat $tmp/genome_to_transcript/$name.intersect.gz |" or die "Cannot open infile\n";
	open OUT,">$tmp/codon_frequency/$name.tmp" or die "Cannot open outfile\n";

	# Keep track of multiples with the same ID
	our @prev_frame;
	our @prev_shift;
	our @prev_codon;
	our @prev_pos;
	our $prev_seq;
	our $prev_gene;
	our $prev_strand;
	our $prev_id;
	our $prev_length;

	while (my $row = <IN>) {
		chomp $row;
		my @line = split "\t",$row;

		# Retrieve data
		# Change here if additional columns are added
		my ($a_chr,$a_start,$a_end,$a_id,$a_seq) = @line[0..4];
		my ($b_chr,$b_start,$b_end,$b_transcript,$b_shift,$b_strand,$b_gene,$b_frame,$b_overlap) = @line[5..$#line];

		# Calculate position in transcript
		my $pos;
		if ($b_strand eq "-") {
			$pos = $b_shift+$b_end-$a_end;
		} else {
			$pos = $b_shift+$a_start-$b_start;
		}

		# Calculate frame according to gff file
		my $frame;
		if ($b_strand eq "-") {
			$frame = $b_end-$a_end-$b_frame;
		} else {
			$frame = $a_start-$b_start-$b_frame;
		}
		my $frame_pos = $frame%3;

		# Determine read length and check if frame is correct for read length
		my $length = length($a_seq);	
		my $pred_shift = $all_shifts[$length][$frame_pos]//"-"; # If length is not in predicted_frame file
		if (!defined($pred_shift)) {
			print "length: $length, frame_pos: $frame_pos\n";
		}
		next if $pred_shift eq "-";
#		next unless $frame_pos eq ($pred_shift%3);	# Skip if not correct frame for length

		if ($b_strand eq "-") {
			$a_seq = rev($a_seq);
		}

		# Crop if flanking reads to the left of TIS
		my $number_of_shifts = $pred_shift;
		while ($number_of_shifts > 0) {
			$a_seq = "N".$a_seq;
			$number_of_shifts--;
		}

		$a_seq .= "NNNNNNNNNNNNNNNNNNNN";	# Avoid errors due to short sequence
		$a_seq = substr $a_seq,0-$number_of_shifts,29;
	
	
		my @codon;
		foreach my $i (-5 .. 3) {
			$codon[$i+5] = substr $a_seq,(15+3*$i),3;
		}

		if (defined($prev_id) && $prev_id ne $a_id) {
			output();	
		}
		# Save current data for later output
		push @prev_frame,$frame_pos;
		push @prev_shift,$pred_shift;
		push @prev_codon,\@codon;
		push @prev_pos,$pos;
	   ($prev_seq,$prev_gene,$prev_strand,$prev_id,$prev_length) = ($a_seq,$b_gene,$b_strand,$a_id,$length);
	}
	close IN;
	output();
	close OUT;	

	sub output {
		my $selected;
		$selected = rand @prev_frame; # Select random mapping if several alternatives

#		print OUT "$prev_seq\t$prev_gene\t$prev_strand\t$prev_pos[$selected]\t$prev_frame[$selected]\t$prev_id";
		foreach my $c (@{ $prev_codon[$selected] }) {
			print OUT "$c\t";
		}
		print OUT "$prev_length\t$prev_frame[$selected]\n";
		undef @prev_codon;
		undef @prev_frame;
		undef @prev_pos;
		undef @prev_shift;
	}

	sub rev {
		my $seq = shift;
		$seq = reverse $seq;
		my $ret = "";
		my @chars = split "",$seq;
		foreach my $c (@chars) {
			$ret .= "A" if $c eq "T";
			$ret .= "C" if $c eq "G";
			$ret .= "G" if $c eq "C";
			$ret .= "T" if $c eq "A";
		}
		return $ret;
	}
}


sub CodonFrequency2 {
	print "Calculating codon frequency (step 2) per read length...\n";

	mkdir "$tmp/codon_frequency/$name";
	SystemBash("rm -f $tmp/codon_frequency/$name/$name.codon.*.txt");

	foreach my $i (1..9) {
   	my @lengths=`cat $tmp/codon_frequency/$name.tmp | cut -f10 | sort | uniq`;
   	foreach my $length (@lengths) {
			chomp $length;
			next if $length > 50 or $length < 20;
			foreach my $frame (0..2) {
	      	SystemBash("cat $tmp/codon_frequency/$name.tmp | awk -v len=$length -v frame=$frame '\$10==len && \$11==frame' | cut -f$i | grep -v N | awk 'length>2' | sort | uniq -c | sort -k2,2 > $tmp/codon_frequency/$name/$name.codon.$i.$length-$frame.txt");
			}
		}
	}
}

sub PlotCodonFrequency {
	my ($folder,$out) = @_;
	print "Plotting codon per read length correlation...\n";
	mkdir "$tmp/$out";
	SystemBash("Rscript --vanilla correlation.R $tmp/$folder $name $genome $tmp/$out");
}

sub CorrelateToBest {
	print "Calculating correlation between read lengths...\n";
	SystemBash("mkdir -p $tmp/correlate_to_best");
	open OUT,">$tmp/correlate_to_best/CorrelateToBest.$name.txt";
	open IN,"$tmp/correlation/$name.correlation.csv" or die "Cannot open $name.correlation.csv\n";
	SystemBash("rm -f $tmp/codon_frequency/$name.selected/*");

	# Extract all names from first row
	chomp(my $row = <IN>);
	my @names = split ",",$row;
	shift @names; # First name is empty, remove it!

	# Extract shift for each read length
	my @shifts;
	foreach my $name (@names) {
		next if $name eq "BG";
		my @list = split "\\.",$name;
		$shifts[$list[1]]++;
	}
	$shifts[0] = "";

	# Save the whole correlation matrix
	my %cor;
	while (my $row = <IN>) {
		chomp $row;
		my @list = split ",",$row;
		my $name = shift @list;
		foreach my $i (0..$#list) {
			$cor{$name}{$names[$i]} = $list[$i];
#print "$name vs $names[$i] = $list[$i]\n";
		}
	}

	my @max_score;
	my @best_read;
	my %scores;
	foreach my $n (@names) {
		next if $n eq "BG";
		my @sorted = sort {$cor{$n}{$b} <=> $cor{$n}{$a}} keys %{ $cor{$n} }; # Sort each other read length by correlation to THIS
		my @tmp = split "\\.",$n;
		my $pos = $tmp[1]; # Extract position (e.g. 1-9; EPA)
		my $score = 0;
		foreach my $i (0..$shifts[$pos]-1) { # Check as many as there are read lenths are this position
			my $key = $sorted[$i];
#			print "$key ".Round($cor{$n}{$key},2)." ";
			$score++ if $key =~ /\.$pos$/;
		}
#		print " ".($score-1)."/".($shifts[$pos]-1)."\n"; # How many of the best correlating read lengths SHOULD correlate?
		if ($score >= ($max_score[$pos]//0)) {
			if ($cor{$best_read[$pos]//"BG"}{"BG"} > $cor{$n}{"BG"} || $score > $max_score[$pos]//0) { # Only replace ties if correlation is lower
				$max_score[$pos] = $score;
				$best_read[$pos] = $n;
			}
		}
		$scores{$pos}{$tmp[0]} = $score;
	}

	my %exclude;
	my %keep;
	
	my %cor_score;
	my @sum_per_pos;
	my @N_per_pos;
	my %flag;

	# Calculate cor_score as correlation at 2,3,7,8 - correlation at 5,6
	print OUT "Calculating score for each read...\n";
	my @all_reads = `cat $tmp/correlation/$name.correlation.csv | sed 1d | cut -d"." -f1 | sort | uniq | grep -v BG`;
	foreach my $read (@all_reads) {
		chomp $read;
		print OUT "$read";
		foreach my $pos (1..9) {
			my $correlation = $cor{"$read.$pos"}{"BG"};
			print OUT "\t".Round($correlation//1,3);

			$cor_score{$read} += ($correlation//0)/5 if ($pos == 2 || $pos == 3 || $pos == 4 || $pos == 7 || $pos == 8);
			$cor_score{$read} -= ($correlation//1)/2 if ($pos == 5 || $pos == 6);

			if (defined($correlation)) {
				$sum_per_pos[$pos] += $correlation;
				$N_per_pos[$pos]++;
			}
		}
		print OUT "\t".Round($cor_score{$read},3);
		print OUT "\n";
		$flag{$read} = "";
		if ($cor_score{$read} < 0) { $exclude{$read}++; $flag{$read} .= "BG,"; }
	}
	print OUT "MEAN";
	foreach my $pos (1..9) {
		print OUT "\t".Round($sum_per_pos[$pos]/$N_per_pos[$pos],3);
	}
	print OUT "\n";

	my @order = sort { $cor_score{$b} <=> $cor_score{$a} } (@all_reads);
	print OUT "Best score for: $order[0]\n";
	my $best = shift @order;
	my $best_bak = $best;
	while ($scores{"5"}{$best} < $max_score[5]*0.5 || $scores{"6"}{$best} < $max_score[6]*0.5) {
		if (@order == 0) { $best = $best_bak; last; } # Default to first in order if none is OK.
		$best = shift @order;
		print OUT "Best score for: $best\n";
	}

	foreach my $read (@all_reads) {
		next if $exclude{$read};
		print OUT "$read\t".Round($cor_score{$read},3);
		my ($scoreSum,$maxSum) = (0,0);
		foreach my $pos (5,6) {
			if (!defined($scores{$pos}{$read})) {
				$flag{$read} .= "MISSING$pos,";
				$exclude{$read}++;
				next;
			}
			my @tmp = split "\\.",$best_read[$pos];
#			$best = $tmp[0];
			my @cor_to_best;
			foreach my $i (1..9) {
				$cor_to_best[$i] = $cor{"$read.$i"}{"$best.$pos"}//0;
			}
			my @order = sort { $cor_to_best[$b] <=> $cor_to_best[$a] } (1..9);
			print OUT "\t[$best]$pos==$order[0],$order[1],$order[2]";
			print OUT "\t$scores{$pos}{$read}/$max_score[$pos]";
			$scoreSum += $scores{$pos}{$read};
			$maxSum += $max_score[$pos];
#			print "\t$scores{$pos}{$read}/$max_score[$pos]\tcor($pos,$names[$#names])=".Round($cor{"$read.$pos"}{"BG"},3)."\tcor($pos,$best_read[$pos])=".Round($cor{"$read.$pos"}{$best_read[$pos]},3)."\tcor($best_read[$pos],BG)=".Round($cor{$best_read[$pos]}{"BG"},3);
			if ($scores{$pos}{$read} < $max_score[$pos]*0.5) {
				$flag{$read} .= "SCORE$pos,";
#				$exclude{$read}++;
			}
			if ($cor{"$read.$pos"}{"BG"} > $cor{"$read.$pos"}{$best_read[$pos]}) {
				$flag{$read} .= "COR$pos,";
#				$exclude{$read}++;
			}
			if ($order[0] != $pos) {
				$flag{$read} .= "ORDER$pos,";
				$exclude{$read}++;
			}
		}
		print OUT "\t$scoreSum/$maxSum";
		if ($scoreSum/$maxSum < 0.5) {
			$flag{$read} .= "SCORE,";
#			$exclude{$read}++;
		}
		
		chop $flag{$read};
		print OUT "\t$flag{$read}\n";
	}

	my @exclude = sort keys %exclude;
	my @keep = @all_reads;
	foreach my $e (@exclude) { @keep = grep {$_ ne $e} @keep; }
#	print "@max_score\n";
	print OUT " .. Read lengths excluded: @exclude\n";
	print OUT " .. Read lengths kept: @keep\n";
	close IN;	
	close OUT;

	# Create links to the files that should be kept
	mkdir "$tmp/codon_frequency/$name.selected";
	foreach my $k (@keep) {
		my @files = `ls $tmp/codon_frequency/$name/$name.codon.*.$k.txt`;
		foreach my $file (@files) {
			chomp $file;
			$file =~ /codon\.(\d+)\.$k\.txt/;
			my $pos = $1;
			SystemBash("ln -s ../$name/$name.codon.$pos.$k.txt $tmp/codon_frequency/$name.selected/$name.codon.$pos.$k.txt");
		}
	}
}

sub FinalCodonFrequency {
	print "Calculate final codon frequency...\n";

	my %counts;
	my @files = `ls $tmp/codon_frequency/$name.selected/*`;
	foreach my $file (@files) {
		chomp $file;
		$file =~ /codon\.(\d+)\.(\d+)-(\d)\.txt/;
		my ($pos,$length,$frame) = ($1,$2,$3);
		open IN,$file;
		while (my $row = <IN>) {
			$row =~ /\s*(\d+) ([ACGT]+)$/;
			my ($count,$codon) = ($1,$2);
			$counts{$codon}{$pos} += $count;
		}
		close IN;
	}

	open OUT,">$out/$name.codons.txt" or die ("Cannot open outfile\n");
	print OUT "Codon";
	foreach my $pos (sort keys %{ $counts{"AAA"} }) { print OUT "\t$pos";	}
	print OUT "\n";
	foreach my $codon (sort keys %counts) {
		print OUT "$codon";
		foreach my $pos (sort keys %{ $counts{$codon} }) {
			print OUT "\t".($counts{$codon}{$pos}//0);
		}
		print OUT "\n";
	}
	close OUT;
}

sub Validate {
	print "Make final figures for validation...\n";
	SystemBash("Rscript --vanilla validate.R $out $name");
}

# Make sure that Bash is used on Ubuntu systems.
sub SystemBash {
  my @args = ( "bash", "-c", shift );
  system(@args);
}

# Print version message and then exit.
sub VersionMessage {
	print "Concur v$version\n";
	exit();
}

# Round x to y decimals.
sub Round {
	return int($_[0]*10**$_[1]+0.5)/10**$_[1];
}

# How common is the most common position?	
sub MaxPercent {
	my @sorted = sort {$a <=> $b} @_;
	my $sum = $sorted[0]+$sorted[1]+$sorted[2];
	return 0 if $sum == 0;
	return Round($sorted[2]/$sum,3);
}

# Position of maximum value in array.
sub IndexOfMax {
	my $index = 0;
	for my $i (0 .. $#_) { $index = $i if  $_[$index] < $_[$i]; }
	return $index;
}

__END__

=head1 NAME

concur.pl - Using Getopt::Long and Pod::Usage

=head1 SYNOPSIS

concur.pl [options] -i alignment.bam -g genome -o output_name

 Mandatory arguments:
   -i --input       Input bam file
   -g --genome      Reference genome [hg38, hg19, mm10, mm9]
   -o --out         Output file name
 Additional options:
   -h --help        Help message & quit
   -m --man         Help message & quit
   -v --version     Version message & quit

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<CONCUR> will read the given input file(s) and do someting
useful with the contents thereof.

=cut
