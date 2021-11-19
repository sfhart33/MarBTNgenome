#!C:\Perl64\bin\perl -w

####################################################################################################
#		Michael J. Metzger
#		Pacific Northwest Research Institute
#		2020-05-06
#
#		Title: fastq_matefinder.pl
#
#		Project: Integration site mapping
#	
#		Input: FASTQ file (LR not CRLF) with list of reads and one FASTQ file with potential mates
#
#		Output: FASTQ file of mate contigs (_mates.fastq)
#				FASTQ file of input reads with mate pairs (_paired.fastq)	
#				FASTQ file of input reads without mate pairs (_unpaired.fastq)	
#
####################################################################################################

use strict;
use warnings;

open (my $log, '>', "fastq_matefinder.log");

my $usage = "perl fastq_matefinder.pl query.fastq potentialmates.fastq\n";

###################
### INPUT FILES ###
###################

my $query_fastq = shift(@ARGV) or die $usage;
my $potentialmates_fastq = shift(@ARGV) or die $usage;


#Take a fastq file and generate a new file with all mates 
	

open (my $query, '<', $query_fastq) or die "Could not open file '$query_fastq'";
open (my $pot_mates, '<', $potentialmates_fastq) or die "Could not open file '$potentialmates_fastq'";
my @query_name = split (/.fastq/, $query_fastq);

open (my $mates, '>', "$query_name[0]_mates.fastq") or die "Could not open output file '$query_name[0]_mates.fastq'";
open (my $query_paired, '>', "$query_name[0]_paired.fastq") or die "Could not open output file '$query_name[0]_paired.fastq'";
open (my $query_unpaired, '>', "$query_name[0]_unpaired.fastq") or die "Could not open output file '$query_name[0]_unpaired.fastq'";


my @q_array;
my @pot_array;
my $namematch = 0;
my $mapcount = 0;
my $unmapcount = 0;
my $q_line = <$query>;
my $lastq_line = $q_line;
my $q_line234;
my $pot_line;
my $pot_line234;

while ($q_line = <$query>) {
	chomp $lastq_line;
	if ($lastq_line =~ m/^@/) {
		@q_array = split (/[ \/]/, $lastq_line);
		open ($pot_mates, '<', $potentialmates_fastq) or die "Could not open output file '$potentialmates_fastq'";
		while ($pot_line = <$pot_mates>) {
			chomp $pot_line;
			if ($pot_line =~ m/^@/) {
				@pot_array = split (/[ \/]/, $pot_line);
				if ($q_array[0] eq $pot_array[0]) {
					if ((($q_array[1] =~ m/^1/) & ($pot_array[1] =~ m/^2/)) | (($q_array[1] =~ m/^2/) & ($pot_array[1] =~ m/^1/))) { 
						print $mates "$pot_line\n";
						$namematch = 1;
						$mapcount = ($mapcount + 1);
						while ($pot_line234 = <$pot_mates>) {
							last if ($pot_line234 =~ m/^@/);
							print $mates "$pot_line234";
						} 
						print $query_paired "$lastq_line\n";
						print $query_paired "$q_line";
						$lastq_line = $q_line;
						while ($q_line = <$query>) {
							last if ($q_line =~ m/^@/);
							$lastq_line = $q_line;
							print $query_paired "$q_line";
						}
						$lastq_line = $q_line;
						last;
					}
				}
			}
		}
	}
	if ($namematch == 0) { 
		print $query_unpaired "$lastq_line\n";
		print $query_unpaired "$q_line";
		$lastq_line = $q_line;
		while ($q_line = <$query>) {
			last if ($q_line =~ m/^@/);
			$lastq_line = $q_line;
			print $query_unpaired "$q_line";
		}
		$lastq_line = $q_line;
		$unmapcount = ($unmapcount + 1);
	}
	$namematch = 0;
}

close $query;
close $pot_mates;
close $mates;
close $query_paired;
close $query_unpaired;

for(\*STDOUT, $log) {print $_ "'$mapcount' reads with correctly matched mates\n'$unmapcount' reads with no match\n";}

for(\*STDOUT, $log) {print $_ "'$query_name[0]_mates.fastq' generated with sequences that are mates of reads in '$query_fastq'\n";}

if ($unmapcount == 0) {
	unlink("$query_name[0]_paired.fastq") or die "Can't delete '$query_name[0]_paired.fastq': $!\n";
	unlink("$query_name[0]_unpaired.fastq") or die "Can't delete '$query_name[0]_unpaired.fastq': $!\n";
	for(\*STDOUT, $log) {print $_ "All reads mapped. Paired and unpaired temp files removed.\n";}
}

exit;


