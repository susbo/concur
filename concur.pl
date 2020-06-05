#!/bin/perl -w
print `date`;
#    CONCUR: quick and robust calculation of codon usage from ribosome profiling data
#    Copyright (C) 	2020	Susanne Bornelöv
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published
#    by the Free Software Foundation, version 3
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

use strict;
use Getopt::Long;
use Pod::Usage;

my $version = "0.9";

GetOptions("version" => sub { VersionMessage() }
          , "help" => \my $help
          , "man" => \my $man
          , 'input=s' => \my $input
          , 'genome=s' => \my $genome
          , 'out=s' => \my $out
          , 'name=s' => \my $name
          , 'size=s' => \my $size
          , 'noR' => \my $noR
          ) or pod2usage(2);

pod2usage(2) if $help or $man;
unless ($input && $genome && $out) { pod2usage("ERROR: Missing mandatory arguments.\n"); pod2usage(2); }

if (!$name) { # Set name to input file name by default.
	my @tmp = split "/",$input;
	$name = $tmp[$#tmp]; $name =~ s/.bam//;
}

# Set size range to 20-50 nt by default.
if (!$size) { $size = "20-50" }
my $min_size; my $max_size;
if ($size =~ /^(\d+)-(\d+)$/) {
	($min_size,$max_size) = ($1,$2);
} else {
	pod2usage("ERROR: Unknown format used for fragment size range.\n"); pod2usage(2);
}

#### Info #####
# Print welcome message and information.

print "Running CONCUR v$version\n\n";
print "### Parameters ###\n";
print "Input file: $input\n";
print "Genome: $genome\n";
print "Fragment size range tested: $size\n";
print "Output folder: $out\n";
print "Output file name: $name\n";
print "Run without R: "; if ($noR) { print "TRUE"; } else { print "FALSE"; } print "\n";
print "##################\n";

#### Running CONCUR ####
my $tmp = "$out/tmp";
mkdir $out;
mkdir $tmp;

# Run each subroutine
GenomeToTranscript();
PeriodicityAtTSS();
PredictFrame();
CodonFrequency();
CodonFrequency2();
PlotCodonFrequency("codon_frequency/$name","correlation") unless $noR;
CorrelateToBest();
PlotCodonFrequency("codon_frequency/$name.selected","correlation.selected") unless $noR;
FinalCodonFrequency();
Validate() unless $noR;

print `date`;

sub GenomeToTranscript {
	print "[Step 1/10] Mapping genomic reads to transcripts...\n";
	mkdir "$tmp/genome_to_transcript";

	# Extract reads with any overlap to exons; consider start for reads on + strand, and end for reads on - strand.
	SystemBash("paste <(bedtools bamtobed -i $input | cut -f1-4,6) <(samtools view $input | cut -f10) | awk -v OFS=\"\t\" '\$5==\"-\" {\$2=\$3-length(\$6)} \$5==\"+\" {\$3=\$2+length(\$6)} {print \$1,\$2,\$3,\$4,\$6,\$5}' | gzip -c > $tmp/genome_to_transcript/$name.bed.gz");
	SystemBash("bedtools intersect -a <(zcat $tmp/genome_to_transcript/$name.bed.gz) -b <(zcat data/$genome.bed.gz) -wo -split | gzip -c > $tmp/genome_to_transcript/$name.intersect.gz");
}

# Calculate periodicity at TSS per read length
sub PeriodicityAtTSS {
	print "[Step 2/10] Calculating periodicity...\n";
	mkdir "$tmp/periodicity";
	my $file = "$tmp/genome_to_transcript/$name.intersect.gz";
#	foreach my $len (20..50) {
	foreach my $len ($min_size..$max_size) {
     SystemBash("zcat $file | awk -v len=\"$len\" 'length(\$5)==len && \$6==\$12' | awk '{if (\$12==\"+\") shift=\$11+\$2-\$8+12; else shift=\$11+\$9-\$3+12; if (shift<13 && shift>-13) print shift;}' | sort | uniq -c | sort -k2,2 -n > $tmp/periodicity/$name.$len.txt");
	}
}

# Predict reading frame for each read length
# This function creates an initial guess of the best shift for each read length and frame
sub PredictFrame {
	print "[Step 3/10] Predicting frame per read length...\n";
	open OUT,">$out/$name.predicted_frame.txt";
	my $max_ratio = 0;
	my $max_i;
	my $max_N = 0;

	print OUT "Length\tShift\tScores\tDiffs\tReads\tSelected\n";

#	foreach my $i (20..50) {
	foreach my $i ($min_size..$max_size) {
		my @counts;	# Counts per position relative to the TIS
		open IN,"$tmp/periodicity/$name.$i.txt";
		while (my $row = <IN>) {
			$row =~ /\s(\d+) ([0-9-]+)$/;
			next if $2 < -20;
			$counts[$2+20] += $1;
		}
		close IN;
		my @total = (0,0,0);	# Counts per reading frame in TIS region -- READS column
		foreach my $j (-20..99) {
			$total[$j%3] += $counts[$j+20]//0;
			my $b = $j;
			$b = $j+1 if (($counts[$b+20]//0) < ($counts[$j+1+20]//0));
			$b = $j+2 if (($counts[$b+20]//0) < ($counts[$j+2+20]//0));
		}
		print OUT "$i\t"; # Fragment length
		my $N = 0; $N += $_ for @total;	# Total number of reads
		my @shift_score = (0)x3; # SCORES column
		my @shift_diff = (0)x3; # DIFFS column
		my @predicted_shift = ("-","-","-");
		# Make a frame prediction if total count is above 1000 reads.
		if ($N>1000) {
			foreach my $j (-6..6) {	# Find predicted shift as first local maxima position at predicted frame
				my $frame = $j%3;
				if ($total[$frame] > 1000 && ($counts[$j+20]//0) > ($counts[$j+20-3]//0) && ($counts[$j+20]//0) > ($counts[$j+20+3]//0)) { # Only consider local maxima (2.1)
					my $score = Round(($counts[$j+20]//0) / (($counts[$j+20-3]//0) + ($counts[$j+20]//0) + ($counts[$j+20+3]//0)+1),3); # Use for second condition in 2.1.
					my $diff = Round(($counts[$j+20]//0)*2 - ($counts[$j+20-9]//0)/9 - ($counts[$j+20-6]//0)/3 - ($counts[$j+20-3]//0) - ($counts[$j+20+3]//0)/3 - ($counts[$j+20+6]//0)/9 - ($counts[$j+20+9]//0)/9 ,0);
					if ($diff >= $shift_diff[$frame] || $shift_diff[$frame] eq "-") { # Only replace it better or not set yet.
						# Keep track of best score and difference for the predicted shift.
						$shift_score[$frame] = $score;
						$shift_diff[$frame] = $diff;
						$predicted_shift[$frame] = $j;
					}
				}
			}
		}
		print OUT "",join(";",@predicted_shift),"\t",join(";",@shift_score),"\t",join(";",@shift_diff)."\t",join(";",@total)."\t";

		my @score = (0,0,0);
		foreach my $f (0..2) {
			$score[$f] = 1 if $shift_score[$f] > 0.5 && $total[$f]>1000;
		}
		print OUT "",join(";",@score)."\n";
	}
	close OUT;
}

# Calculate codon frequencies
sub CodonFrequency {
	print "[Step 4/10] Calculating codon frequency per read length...\n";
	# Retrieve codon frequencies
	mkdir "$tmp/codon_frequency";

	# Start by reading all predicted shifts
	my @all_shifts;
	open IN,"$out/$name.predicted_frame.txt";
	<IN>; # Skip header line
	while (my $row = <IN>) {
		chomp $row;
		my @line = split "\t",$row;
		my @shift = split ";",$line[1]; # Predicted best shift
		my @selected = split ";",$line[$#line]; # Selected or not?
		foreach my $i (0..2) {
			$shift[$i] = "-" if $selected[$i] == 0;
			$all_shifts[$line[0]][$i] = $shift[$i];
		}
	}
	close IN;

	open IN,"zcat $tmp/genome_to_transcript/$name.intersect.gz |" or die "Cannot open infile\n";
	open OUT,"| gzip -c > $tmp/codon_frequency/$name.tmp.gz" or die "Cannot open outfile\n";

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
		my ($a_chr,$a_start,$a_end,$a_id,$a_seq, $a_strand) = @line[0..5];
		my ($b_chr,$b_start,$b_end,$b_transcript,$b_shift,$b_strand,$b_gene,$b_frame,$b_overlap) = @line[6..$#line];

		next if $a_strand ne $b_strand; # Only consider strand-specific matches.

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
		next if $pred_shift eq "-";

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
	print "[Step 5/10] Calculating codon frequency (step 2) per read length...\n";

	mkdir "$tmp/codon_frequency/$name";
	SystemBash("rm -f $tmp/codon_frequency/$name/$name.codon.*.txt");

	foreach my $i (1..9) {
   	my @lengths=`zcat $tmp/codon_frequency/$name.tmp.gz | cut -f10 | sort | uniq`;
   	foreach my $length (@lengths) {
			chomp $length;
			next if $length > 50 or $length < 20;
			foreach my $frame (0..2) {
	      	SystemBash("zcat $tmp/codon_frequency/$name.tmp.gz | awk -v len=$length -v frame=$frame '\$10==len && \$11==frame' | cut -f$i | grep -v N | awk 'length>2' | sort | uniq -c | sort -k2,2 > $tmp/codon_frequency/$name/$name.codon.$i.$length-$frame.txt");
			}
		}
	}
}

sub PlotCodonFrequency {
	my ($folder,$out) = @_;
	print "[Step 6/10] Plotting codon per read length correlation...\n" if $out !~ "selected";
	print "[Step 8/10] Plotting codon per read length correlation...\n" if $out =~ "selected";
	mkdir "$tmp/$out";
	SystemBash("Rscript --vanilla scripts/correlation.R $tmp/$folder $name $genome $tmp/$out");
}

sub CorrelateToBest {
	print "[Step 7/10] Calculating correlation between read lengths...\n";
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
			$score++ if $key =~ /\.$pos$/;
		}
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
		if ($cor_score{$read} < 0) { $exclude{$read}++; $flag{$read} .= "BG,"; } # Filter 2.2.1
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
			if (!defined($scores{$pos}{$read})) { # This should not happen...
				$flag{$read} .= "MISSING$pos,";
				$exclude{$read}++;
				next;
			}
			my @tmp = split "\\.",$best_read[$pos];
			my @cor_to_best;
			foreach my $i (1..9) {
				$cor_to_best[$i] = $cor{"$read.$i"}{"$best.$pos"}//0;
			}
			my @order = sort { $cor_to_best[$b] <=> $cor_to_best[$a] } (1..9);

			my @cor_to_best_opp; # Check that is does not correlate better to another position from best
			foreach my $i (1..9) {
				$cor_to_best_opp[$i] = $cor{"$read.$pos"}{"$best.$i"}//0;
			}
			my @order_opp = sort { $cor_to_best_opp[$b] <=> $cor_to_best_opp[$a] } (1..9);

			print OUT "\t[$best]$pos==$order[0],$order[1],$order[2]";
			print OUT "\t[$read]$order_opp[0],$order_opp[1],$order_opp[2]==$pos";
			print OUT "\t$scores{$pos}{$read}/$max_score[$pos]";
			$scoreSum += $scores{$pos}{$read};
			$maxSum += $max_score[$pos];
			if ($scores{$pos}{$read} < $max_score[$pos]*0.5) { # Filter 2.2.3
				$flag{$read} .= "SCORE$pos,";
				$exclude{$read}++;
			}
			if ($cor{"$read.$pos"}{"BG"} > $cor{"$read.$pos"}{$best_read[$pos]}) { # Not used currently
				$flag{$read} .= "COR$pos,";
#				$exclude{$read}++;
			}
			if ($order[0] != $pos) { # Filter 2.2.2
				$flag{$read} .= "ORDER$pos,";
				$exclude{$read}++;
			}
			if ($order_opp[0] != $pos) { # Similar to 2.2.2, not used currently
				$flag{$read} .= "OPP$pos,";
#				$exclude{$read}++;
			}
		}
		print OUT "\t$scoreSum/$maxSum"; # Not used currently (combines the score for P and A sites)
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
	print "[Step 9/10] Calculate final codon frequency...\n";

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
	print "[Step 10/10] Make final figures for validation...\n";
	SystemBash("Rscript --vanilla scripts/validate.R $out $name $genome");
}

# This function makes sure that Bash is used on Ubuntu systems.
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

__END__

=head1 NAME

concur.pl - Using Getopt::Long and Pod::Usage

=head1 SYNOPSIS

concur.pl -i alignment.bam -g genome -o output_folder [-n output_name] [-s min_size-max_size] [--noR]

 Mandatory arguments:
   -i --input       Input bam file
   -g --genome      Reference genome [hg38, hg19, mm10, mm9, sc3, rn9]
   -o --out         Output folder name
 Optional arguments:
   -n --name        Output file name [default: input file name]
	-s --size        Fragment size range to use [default: 20-50]
	--noR            Don't run any function that require R [default: FALSE]
 Additional options:
   -h --help        Print help message and exit
   -m --man         Print help message and exit
   -v --version     Print version and exit

=head1 OPTIONS

=over 8

=item B<--input>

Input bam file

=item B<--genome>

Reference genome [hg38, hg19, mm10, mm9, sc3, rn9]

=item B<--out>

Output folder name

=item B<--name>

Output file name [default: input file name]

=item B<--size>

Fragment size range to use for analysis [default: 20-50]

=item <--noR>

Don't run any function that require R [default: FALSE]

=item B<--help>

Print a brief help message and exit

=item B<--man>

Print a brief help message and exit

=item B<--version>

Print version message and exit

=back

=head1 DESCRIPTION

B<CONCUR> will read the given input bam file, detect the best frame
for identifying the A-site, and return a table with codon usage.

Please see https://github.com/susbo/concur for more details.

=cut
