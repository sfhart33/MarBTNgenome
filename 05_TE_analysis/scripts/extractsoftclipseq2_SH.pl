#!C:\Perl64\bin\perl -w

####################################################################################################
#		Michael J. Metzger
#		Pacific Northwest Research Institute
#		2019-12-16
#
#		Title: extractsoftclipseq.pl
#
#		Project: MappingTEs
#	
#		Input: SAM file of illumina seq mapped to LTR
#			need to input TElength into script below
#
#		Output: fastq file of soft-clipped portions of sam, and error file (S>1)
#			
#		Dependency: 
#
####################################################################################################

use strict;
use warnings;


my $usage = "perl extractsoftclipseq.pl mappedfile.sam\n";

###################
### INPUT FILES ###
###################

my $samfile = shift(@ARGV) or die $usage;

#
open (my $sam, '<', $samfile) or die "Could not open file '$samfile'";
open (my $fastqup, '>', "softclips_up.fastq") or die "Could not open output file";
open (my $fastqdown, '>', "softclips_down.fastq") or die "Could not open output file";
open (my $internalmap, '>', "internalmap.sam") or die "Could not open output file";
my @samarray;
my $samarray;
my @CIGARsplit;
my $CIGARsplit;
my @CIGARcount;
my $CIGARcount;
my $seq;
my $qscore;
my $flanklen;
my $seqname;
my $TElength = 177; #set based on length of TE reference sequence
print "TE length is $TElength bp.\n";
print "If this is incorrect, input the correct length  of \$TElength in the script.\n";
my $matchlen;
my $CIGARlet;
my $CIGARnum;
my $readnum;
while (my $sam_line = <$sam>) {
	chomp $sam_line;
	if ($sam_line !~ m/^@/) {
		@samarray = split (/\t/, $sam_line);
		if ($samarray[1] & 0x40) { # first read in a pair
			$readnum = 1;
			$seqname = "$samarray[0]*$readnum";
		} else {
			if ($samarray[1] & 0x80) { # second read in a pair
				$readnum = 2;
				$seqname = "$samarray[0]*$readnum";
			} else { # no paired data
			$seqname = $samarray[0];
			}
		}
		if ($samarray[4] != 0) { # continues with all mapped reads
			if ($samarray[5] =~ m/S{1}/) { # continues with all reads with 1 soft clip
				if ($samarray[5] =~ m/^\d+S/) {
					# soft clip is in front, save $CIGARsplit[0] number of bases from front
					if ($samarray[3] == 1) { 
						# confirms matching to beg of TE
						@CIGARsplit = split(/S/, $samarray[5]);
						$flanklen = ($CIGARsplit[0]+0);
						$seq = substr($samarray[9], 0, $flanklen);
						$qscore = substr($samarray[10], 0, $flanklen);
						print $fastqup "\@$seqname\n$seq\n\+\n$qscore\n";
					} else {
						print $internalmap "$sam_line\n";
					}
				} else {
				if ($samarray[5] =~ m/\d+S$/) { 
					#soft clip is in back, save CIGARsplit[-1] number of bases from back
					@CIGARcount = split(/(M|S|D|N|I)/, $samarray[5]);
					$matchlen = 0;
					while (scalar(@CIGARcount) != 0) {
						$CIGARnum = shift @CIGARcount;
						$CIGARlet = shift @CIGARcount;
					if ($CIGARlet =~ m/(M|D|N)/) {
						$matchlen = $matchlen + $CIGARnum;
						}
					} # generates match length based on CIGAR	
					if ($TElength = $matchlen + $samarray[3] - 1) { 
						# confirms matching to end of TE
						@CIGARsplit = split(/[M|S]/, $samarray[5]);
						$flanklen = ($CIGARsplit[-1]+0);
						$seq = substr($samarray[9], -$flanklen);
						$qscore = substr($samarray[10], -$flanklen);
						print $fastqdown "\@$seqname\n$seq\n\+\n$qscore\n";
					} else {
						print $internalmap "$sam_line\n";
					}
				}
			}	
			}	else {
			if ($samarray[5] =~ m/S{2,}/) {
				#more than one soft clip
				print "more than one soft clip";
				open (my $doubleclip, '>>', "doubleclip.sam") or die "Could not open output file";
				print $doubleclip "$sam_line\n";
				close $doubleclip;
			}	
			}
			}
		}
	}  


close $sam;
close $fastqup;
close $fastqdown;
close $internalmap;

exit;


