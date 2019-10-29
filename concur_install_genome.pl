#!/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;

my $version = "1.0";

GetOptions("version" => sub { VersionMessage() }
          , "help" => \my $help
          , "man" => \my $man
          , 'gtf=s' => \my $gtf
          , 'fasta=s' => \my $fasta
          , 'short=s' => \my $short
          , 'recalculate' => \my $recalculate
          , 'pcg' => \my $pcg
          ) or pod2usage(2);

pod2usage(2) if $help or $man;
unless ($fasta && $gtf && $short) { pod2usage("Missing mandatory arguments."); pod2usage(2); }
unless ($recalculate) { $recalculate = 0; }

#### Info #####
# Print welcome message and information.

print "Running CONCUR v$version\n\n";
print "### Parameters ###\n";
print "GTF file: $gtf\n";
print "Fasta file: $fasta\n";
print "Short genome name: $short\n";
print "Recalculate frame: ".flag($recalculate)."\n";
print "Protein-coding genes (Gencode) format instead of cds: ".flag($pcg)."\n";
print "##################\n";

## Perform the installation

# Step 1: Convert gtf/gff file to bed
# Re-sort the gtf/gff since it may not be sorted correctly.
mkdir "data";
if ($gtf =~ /.gz$/) {
	open IN,"zcat $gtf | grep -v \"^#\" | awk '\$3==\"CDS\"' | awk -v OFS=\"\\t\" '{if (\$7==\"-\"){\$4=-\$4; \$5=-\$5;}};{print};' | sort -k1,1 -k4,4n | awk -v OFS=\"\\t\" '{if (\$7==\"-\"){\$4=-\$4; \$5=-\$5;}}; {print};' |" or die "Cannot open $gtf\n";
} else { 
	open IN,"cat $gtf | grep -v \"^#\" | awk '\$3==\"CDS\"' | awk -v OFS=\"\\t\" '{if (\$7==\"-\"){\$4=-\$4; \$5=-\$5;}};{print};' | sort -k1,1 -k4,4n | awk -v OFS=\"\\t\" '{if (\$7==\"-\"){\$4=-\$4; \$5=-\$5;}}; {print};' |" or die "Cannot open $gtf\n";
}
#open OUT,'>:gzip',"data/$short.bed.gz" or die "Cannot write to data/$short.bed.gz\n";
open OUT,"| gzip -c > data/$short.bed.gz" or die "Cannot write to data/$short.bed.gz\n";
my %offset;
while (my $row = <IN>) {
   next if $row =~ /^#/; # Not needed atm.
   chomp $row;
   my @line = split "\t",$row;
	$line[8] = join(" ",@line[8..$#line]); # If needed, replace tabs with space in INTO field.
   next if $line[2] ne "CDS"; # Not needed atm.
   my ($start,$end) = ($line[3]-1,$line[4]);
   print OUT "$line[0]\t$start\t$end\t";

   my $transcript_id;
	if ($line[8] =~ /transcript_id[:=\s]([^;]+);/) {
		$transcript_id = $1;
	} else {
	   $line[8] =~ /transcript[:=\s]([^;]+);/;
		$transcript_id = $1;
	}
	$transcript_id =~ s/"//g;
   print OUT "$transcript_id\t";

	# Offset
   my $offset = $offset{$transcript_id}//0;
   print OUT "$offset\t";
   $offset{$transcript_id} += $end-$start;

	# Strand
   print OUT "$line[6]\t";

   my $gene_id;
	if ($line[8] =~ /gene_id[:= ]([^;]+);/) {
	   $gene_id = $1;
	} elsif ($line[8] =~ /gene[:= ]([^;]+);/) {
	   $gene_id = $1;
	} elsif ($line[8] =~ /CDS[:= ]([^;]+);/) {
	   $gene_id = $1;
	}
	$gene_id =~ s/"//g;

	# Frame
	my $frame = $line[7];
	if ($recalculate) {
		$frame = (3-($offset{$transcript_id}-($end-$start))%3)%3;
	}


   print OUT "$gene_id\t$frame"; # Gene_id and frame_pos

   print OUT "\n";
}
close IN;
close OUT;

# Step 2: Calculate codon frequencies
if ($fasta =~ /.gz$/) {
	open IN,"zcat $fasta |" or die "Cannot open $fasta\n";
} else {
	open IN,"$fasta" or die "Cannot open $fasta\n";
}
open OUT,">data/$short.bg.txt" or die "Cannot write open/$short.txt";

my %counts;
my @bases = ("A","C","G","T");
my @codons;
my ($used,$excluded) = (0,0);
foreach my $a (@bases) { foreach my $b (@bases) { foreach my $c (@bases) { push @codons, "$a$b$c"; }}}

if (!$pcg) { # Default action; only cds in fasta
	while (my $row = <IN>) {
   	if ($row =~ /^>/) {
			next;
   	} else {
	      chomp $row;
   	   if (length($row)%3 != 0) {
				$excluded++;
         	next;
	      } else {
				$used++;
	   	   my @line = $row =~ /.{3}/g;
   	   	foreach my $e (@line) { 
	   	      $counts{$e}++;
				}
   	   }
	   }
	}
} else { # Gencode format; whole pcg in fasta
	my ($start,$end);
	my $seq = "";
	while (my $row = <IN>) {
   	chomp $row;
	   if ($row =~ /CDS/) {
      	if (length($seq) > 0) { # Do nothing for first row; otherwise do calcs
	         $seq = substr $seq,$start-1,$end-$start+1;
#print length($seq)." [".(length($seq)%3)."] ";
   	      if (length($seq)%3 != 0) {
					$excluded++;
#print "$start-$end";
				} else {
					$used++;
	            my @line = $seq =~ /.{3}/g;
		         foreach my $e (@line) {
   		         $counts{$e}++;
					}
      	   }
	      }
   	   # Next seq
      	$row =~ /CDS:(\d+)-(\d+)/;
	      ($start,$end) = ($1,$2);
      	$seq = "";
#print "\n";
	   }
   	else {
      	$seq = $seq.$row;
#print length($seq)." [".(length($seq)%3)."] ";
	   }
	}
}

##############
foreach my $k (@codons) { print OUT "$k\t$counts{$k}\n"; }
close IN;
close OUT;

# Print flags as true or false
sub flag {
	return "yes" if $_[0];
	return "no";
}

# Print version message and then exit.
sub VersionMessage {
	print "Concur v$version\n";
	exit();
}

print "\nFinished.\n";
print "$used transcripts were used and $excluded were excluded (length not a multiple of 3).\n";
print "$short.bed.gz and $short.bg.txt have been created in the directory \"data\".\n";
print "Make sure that those are present in the folder where you are running CONCUR.\n";

__END__

=head1 NAME

concur_install_genome.pl - Using Getopt::Long and Pod::Usage

=head1 SYNOPSIS

perl  concur_install_genome.pl -f cds.fa -g annotation.gtf -s xx1

 Mandatory arguments:
   -f --fasta          Fasta file with CDS sequences
   -g --gtf            GTF/GFF file
   -s --short          Short genome name [such as hg38, hg19, mm10, mm9, sc3]
 Advanced options:
   -r --recalculate    Recalculate frame instead of using column 8 in GTF/GFF
 Additional options:
   -h --help           Help message & quit
   -m --man            Help message & quit
   -v --version        Version message & quit

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
