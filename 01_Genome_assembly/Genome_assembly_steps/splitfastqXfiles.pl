#!C:\Perl64\bin\perl -w

####################################################################################################
#		Michael J. Metzger
#		Pacific Northwest Research Institute
#		2020-05-06
#
#		Title: splitfastqXfiles.pl
#
#		Project: Integration site mapping
#	
#		Input: FASTQ file and # of subfiles to split it into (2-10)
#
#		Output: FASTQ files (.A.fastq, .B.fastq)
#			Really just reads 4 lines into each file sequentially.
#
####################################################################################################

use strict;
use warnings;

open (my $log, '>', "splitfastqXfiles.log");

my $usage = "perl splitfastqXfiles.pl input.fastq #ofsplits\n";

###################
### INPUT FILES ###
###################

my $input_fastq = shift(@ARGV) or die $usage;
my $num = shift(@ARGV) or die $usage;
my $split_num = int($num);

#Take a fastq files and split into a number of smaller files sequentially (read 1 in file A, read 2 in file B, etc)

my $output1;
my $output2;
my $output3;
my $output4;
my $output5;
my $output6;
my $output7;
my $output8;
my $output9;
my $output10;
	

if (($split_num < 2) | ($split_num > 10)) {
	print "number should be between 2-10";
	exit;
}
open (my $input, '<', $input_fastq) or die "Could not open file '$input_fastq'";

my @input_name = split (/.fastq/, $input_fastq);
open ($output1, '>', "$input_name[0].A.fastq") or die "Could not open output file";
open ($output2, '>', "$input_name[0].B.fastq") or die "Could not open output file";
if ($split_num > 2) {
	open ($output3, '>', "$input_name[0].C.fastq") or die "Could not open output file";
	if ($split_num > 3) {
		open ($output4, '>', "$input_name[0].D.fastq") or die "Could not open output file";
		if ($split_num > 4) {
			open ($output5, '>', "$input_name[0].E.fastq") or die "Could not open output file";
			if ($split_num > 5) {
				open ($output6, '>', "$input_name[0].F.fastq") or die "Could not open output file";
				if ($split_num > 6) {
					open ($output7, '>', "$input_name[0].G.fastq") or die "Could not open output file";
					if ($split_num > 7) {
					open ($output8, '>', "$input_name[0].H.fastq") or die 	"Could not open output file";
						if ($split_num > 8) {
							open ($output9, '>', "$input_name[0].I.fastq") or die "Could not open output file";
							if ($split_num > 9) {
								open ($output10, '>', "$input_name[0].J.fastq") or die "Could not open output file";
							}
						}
					}
				}
			}
		}
	}
}

my @matcharray;
my $namematch = 0;
while (my $q_line = <$input>) {
	print $output1 "$q_line";
	$q_line = <$input>;
	print $output1 "$q_line";
	$q_line = <$input>;
	print $output1 "$q_line";
	$q_line = <$input>;
	print $output1 "$q_line";
	$q_line = <$input>;
	last if eof;
	print $output2 "$q_line";
	$q_line = <$input>;
	print $output2 "$q_line";
	$q_line = <$input>;
	print $output2 "$q_line";
	$q_line = <$input>;
	print $output2 "$q_line";
	last if eof;
	if ($split_num > 2) {
		$q_line = <$input>;
		print $output3 "$q_line";
		$q_line = <$input>;
		print $output3 "$q_line";
		$q_line = <$input>;
		print $output3 "$q_line";
		$q_line = <$input>;
		print $output3 "$q_line";
		last if eof;
		if ($split_num > 3) {
			$q_line = <$input>;
			print $output4 "$q_line";
			$q_line = <$input>;
			print $output4 "$q_line";
			$q_line = <$input>;
			print $output4 "$q_line";
			$q_line = <$input>;
			print $output4 "$q_line";
			last if eof;
			if ($split_num > 4) {
				$q_line = <$input>;
				print $output5 "$q_line";
				$q_line = <$input>;
				print $output5 "$q_line";
				$q_line = <$input>;
				print $output5 "$q_line";
				$q_line = <$input>;
				print $output5 "$q_line";
				last if eof;
				if ($split_num > 5) {
					$q_line = <$input>;
					print $output6 "$q_line";
					$q_line = <$input>;
					print $output6 "$q_line";
					$q_line = <$input>;
					print $output6 "$q_line";
					$q_line = <$input>;
					print $output6 "$q_line";
					last if eof;
					if ($split_num > 6) {
						$q_line = <$input>;
						print $output7 "$q_line";
						$q_line = <$input>;
						print $output7 "$q_line";
						$q_line = <$input>;
						print $output7 "$q_line";
						$q_line = <$input>;
						print $output7 "$q_line";
						last if eof;
						if ($split_num > 7) {
							$q_line = <$input>;
							print $output8 "$q_line";
							$q_line = <$input>;
							print $output8 "$q_line";
							$q_line = <$input>;
							print $output8 "$q_line";
							$q_line = <$input>;
							print $output8 "$q_line";
							last if eof;
							if ($split_num > 8) {
								$q_line = <$input>;
								print $output9 "$q_line";
								$q_line = <$input>;
								print $output9 "$q_line";
								$q_line = <$input>;
								print $output9 "$q_line";
								$q_line = <$input>;
								print $output9 "$q_line";
								last if eof;
								if ($split_num > 9) {
									$q_line = <$input>;
									print $output10 "$q_line";
									$q_line = <$input>;
									print $output10 "$q_line";
									$q_line = <$input>;
									print $output10 "$q_line";
									$q_line = <$input>;
									print $output10 "$q_line";
									close $output10;
								}
							}
						}
					}
				}
			}
		}
	}
}
close $input;
close $output1;
close $output2;
close $output3 if defined($output3);
close $output4 if defined($output4);
close $output5 if defined($output5);
close $output6 if defined($output6);
close $output7 if defined($output7);
close $output8 if defined($output8);
close $output9 if defined($output9);
close $output10 if defined($output10);
for(\*STDOUT, $log) {print $_ "$input_fastq has been split into $split_num new fastq files\n";}

exit;


