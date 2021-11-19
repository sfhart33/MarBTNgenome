#!C:\Perl64\bin\perl -w

####################################################################################################
#		Michael J. Metzger
#		Pacific Northwest Research Institute
#		2020-06-28
#
#		Title: intflanktorefSAMtoBEDplus3.pl
#
#		Project: MappingTEs
#	
#		Input: SAM files of sequence flanking TE mapped to genome
#		mapped_upflanks.sam
#		mapped_downflanks.sam
#		output_prefix (optional)
#			
#
#		Output: BED-like file of primary reads with integration locations 
#		(intreads.bed) and another with supplementary matches (intreads.sup.bed)
#			
#		Dependency: 
#
####################################################################################################

use strict;
use warnings;


my $usage = "perl intflanktorefSAMtoBEDplus2.pl mapped_upflanks.sam mapped_downflanks.sam output_prefix\n";

###################
### INPUT FILES ###
###################

my $samfileup = shift(@ARGV) or die $usage;
my $samfiledown = shift(@ARGV) or die $usage;
my $outpre = shift(@ARGV);
my $outfile;
my $outfilesup; #sup files include reads with 0x800 flag
if (not defined $outpre) {
	$outfile = "intreads.bed";
	$outfilesup = "intreads.sup.bed"; 
} else {
	$outfile = $outpre . ".intreads.bed";
	$outfilesup = $outpre . ".intreads.sup.bed";
}
#
open (my $samup, '<', $samfileup) or die "Could not open file '$samfileup'";
open (my $samdown, '<', $samfiledown) or die "Could not open file '$samfiledown'";
open (my $intBED, '>', $outfile) or die "Could not open output file";
open (my $intBEDsup, '>', $outfilesup) or die "Could not open output file";
my @samarray;
my $samarray;
my @CIGARcount;
my $CIGARcount;
my $matchlen = 0;
my $flankmatch;
my $pos1;
my $start;
my $end;
my $CIGARlet;
my $CIGARnum = 0;
while (my $sam_line = <$samdown>) { #mapping of downstream reads
	chomp $sam_line;
	if ($sam_line !~ m/^@/) { #skip header
		@samarray = split (/\t/, $sam_line);
		if ($samarray[1] & 0x10) { #downstream reads that map reverse to genome
			@CIGARcount = split(/(M|S|D|N|I|H)/, $samarray[5]);
			$CIGARnum = shift @CIGARcount;
			$CIGARlet = shift @CIGARcount;
			$matchlen = $CIGARnum; #initial M length
			while (scalar(@CIGARcount) != 0) { #loop for complete match length
				$CIGARnum = shift @CIGARcount;
				$CIGARlet = shift @CIGARcount;
				if ($CIGARlet =~ m/(M|D|N)/) {
					$matchlen = $matchlen + $CIGARnum;
				}
			}
			if ($CIGARlet =~ m/M/) { #confirm match goes to TE junction
				$flankmatch = $samarray[3] + $matchlen - 1; #1-indexed first base of flanking reads
				$pos1 = $samarray[3] + $matchlen - 1; #1-indexed start site of int
				$start = $samarray[3] + $matchlen - 1 - 4 - 1; #0-indexed bed-style start of TSD
				$end = $samarray[3] + $matchlen - 1 + 1 - 1; #0-indexed bed-style end of TSD
				if ($samarray[1] & 0x800) {
					print $intBEDsup "$samarray[2]\t$start\t$end\t$samarray[0]\t$samarray[4]\t\-\tdown\t$flankmatch\t$samarray[2]_$pos1\_R\t\n";
				} else {
					print $intBED "$samarray[2]\t$start\t$end\t$samarray[0]\t$samarray[4]\t\-\tdown\t$flankmatch\t$samarray[2]_$pos1\_R\t\n";
				}
			}
		} else { #downstream reads that map forward relative to genome
			@CIGARcount = split(/(M|S|D|N|I|H)/, $samarray[5]);
			$CIGARnum = shift @CIGARcount;
			$CIGARlet = shift @CIGARcount;
			if ($CIGARlet =~ m/M/) { #confirm match starts at TE junction
				$flankmatch = $samarray[3]; #1-indexed first base of flanking reads
				$pos1 = $samarray[3]; #1-indexed start site of int 
				$start = $samarray[3] - 1; #0-indexed bed-style start of TSD
				$end = $samarray[3] + 5 - 1; #0-indexed bed-style end of TSD
				if ($samarray[1] & 0x800) {
					print $intBEDsup "$samarray[2]\t$start\t$end\t$samarray[0]\t$samarray[4]\t\+\tdown\t$flankmatch\t$samarray[2]_$pos1\_F\t\n";
				} else {
					print $intBED "$samarray[2]\t$start\t$end\t$samarray[0]\t$samarray[4]\t\+\tdown\t$flankmatch\t$samarray[2]_$pos1\_F\t\n";
				}
			}
		}
	}	
}
close $samdown;

#do it all again with upstream reads	
while (my $sam_line = <$samup>) { #mapping of upstream reads
	chomp $sam_line;
	if ($sam_line !~ m/^@/) {
		@samarray = split (/\t/, $sam_line);
		if ($samarray[1] & 0x10) { #upstream reads that map reverse to genome
			@CIGARcount = split(/(M|S|D|N|I|H)/, $samarray[5]);
			$CIGARnum = shift @CIGARcount;
			$CIGARlet = shift @CIGARcount;
			if ($CIGARlet =~ m/M/) { #confirm match starts at TE junction
				$flankmatch = $samarray[3]; #1-indexed first base of flanking reads
				$pos1 = $samarray[3] + 4; #1-indexed start site of int
				$start = $samarray[3] - 1; #0-indexed bed-style start of TSD
				$end = $samarray[3] + 5 - 1; #0-indexed bed-style end of TSD
				if ($samarray[1] & 0x800) {
					print $intBEDsup "$samarray[2]\t$start\t$end\t$samarray[0]\t$samarray[4]\t\-\tup\t$flankmatch\t$samarray[2]_$pos1\_R\t\n";
				} else {
					print $intBED "$samarray[2]\t$start\t$end\t$samarray[0]\t$samarray[4]\t\-\tup\t$flankmatch\t$samarray[2]_$pos1\_R\t\n";
				}
			}
		} else { #upstram reads that map forward relative to genome
			@CIGARcount = split(/(M|S|D|N|I|H)/, $samarray[5]);
			$CIGARnum = shift @CIGARcount;
			$CIGARlet = shift @CIGARcount;
			$matchlen = $CIGARnum; #initial M length
			while (scalar(@CIGARcount) != 0) {
				$CIGARnum = shift @CIGARcount;
				$CIGARlet = shift @CIGARcount;
				if ($CIGARlet =~ m/(M|D|N)/) {
					$matchlen = $matchlen + $CIGARnum;
				}
			} 
			if ($CIGARlet =~ m/M/) { #confirm match goes to TE junction
				$flankmatch = $samarray[3] + $matchlen - 1; #1-indexed first base of flanking reads
				$pos1 = $samarray[3] + $matchlen - 1 - 4; #1-indexed start site of int 
				$start = $samarray[3] + $matchlen - 1 - 4 - 1; #0-indexed bed-style start of TSD
				$end = $samarray[3] + $matchlen - 1 + 1 - 1; #0-indexed bed-style end of TSD
				if ($samarray[1] & 0x800) {
					print $intBEDsup "$samarray[2]\t$start\t$end\t$samarray[0]\t$samarray[4]\t\+\tup\t$flankmatch\t$samarray[2]_$pos1\_F\t\n";
				} else {
					print $intBED "$samarray[2]\t$start\t$end\t$samarray[0]\t$samarray[4]\t\+\tup\t$flankmatch\t$samarray[2]_$pos1\_F\t\n";
				}
			}
		}
	}	
}
close $samup;
close $intBED;
close $intBEDsup;

exit;


