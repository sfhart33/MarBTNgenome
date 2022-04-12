#!C:\Perl64\bin\perl -w

####################################################################################################
#		Michael J. Metzger
#		Pacific Northwest Research Institute
#		2020-06-13
#
#		Title: two_fasta_prefix_compile_firstX.pl
#
#		Project: Integration site mapping
#	
#		Input: Two FASTA files
#
#		Output: Single FASTA file with combined data and prefixes on names
#
####################################################################################################

use strict;
use warnings;

open (my $log, '>', "two_fasta_prefix_compile_firstX.log");

my $usage = "perl two_fasta_prefix_compile.pl sequences1.fasta sequences2.fasta prefix1 prefix2 outfile.fasta numberofsequences\n";

###################
### INPUT FILES ###
###################

my $seq1file = shift(@ARGV) or die $usage;
my $seq2file = shift(@ARGV) or die $usage;
my $seq1pre = shift(@ARGV) or die $usage;
my $seq2pre = shift(@ARGV) or die $usage;
my $outfile = shift(@ARGV) or die $usage;
my $seqnum =  shift(@ARGV) or die $usage;

#Take list of contig names and generate a FASTA file from reference
#open each file and write to a single new file with prefix	

open (my $seq1fasta, '<', $seq1file) or die "Could not open file '$seq1file'";
open (my $seqout, '>', $outfile) or die "Could not open output file";

my $seq1count = 0;
while (my $seqline = <$seq1fasta>) {
	chomp $seqline;
	if ($seqline =~ m/^>/) {
		last if ($seq1count >= $seqnum);
		$seqline =~ s/^>/>$seq1pre/;
		print $seqout "$seqline\n";
		$seq1count = ($seq1count+1);
		for(\*STDOUT, $log) {print $_ "$seqline\n";}
	} else {
	print $seqout "$seqline\n";
	}
}

close $seq1fasta;

open (my $seq2fasta, '<', $seq2file) or die "Could not open file '$seq2file'";

my $seq2count = 0;
while (my $seqline = <$seq2fasta>) {
	chomp $seqline;
	if ($seqline =~ m/^>/) {
		last if ($seq2count >= $seqnum);
		$seqline =~ s/^>/>$seq2pre/;
		print $seqout "$seqline\n";
		$seq2count = ($seq2count+1);
		for(\*STDOUT, $log) {print $_ "$seqline\n";}
	} else {
	print $seqout "$seqline\n";
	}
}

for(\*STDOUT, $log) {print $_ "$seq1count sequences from $seq1file with prefix $seq1pre\n";}
for(\*STDOUT, $log) {print $_ "$seq2count sequences from $seq2file with prefix $seq2pre\n";}

close $seq2fasta;
close $seqout;

exit;


