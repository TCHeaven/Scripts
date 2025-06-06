#!/usr/bin/perl -w

use strict;

# a program to read in every four lines of a file count the number of A,T,G,C, & N's
# in the second line & use these to estimate genome coverage.
# Andrew Armitage 07/03/14




### define usage of the program ###
my ($USAGE) = "count_nucleotides.pl -i seq_data.fastq -i seq_data2.fastq(optional) -g estimated genome size (Mb)";
my ($size);
my ($input);
my ($input2);
my ($total_nuc_count) = 0;
my ($est_cov) = 0;


###################
#	Collect inputs
###################

#	If input values are not entered, or contain whitespace print error messages and exit.
unless(@ARGV) {
	print $USAGE;
	exit;
}

if($ARGV[3] =~ /^\s*$/) {
	print "Entry values undefined, please adhere to: $USAGE\n";
	exit;
}

# Collect input values from the stating commands. Unless multiple inputs are entered
# set $input2 to "empty".
# When collecting estimated genome size, convert this to bp rather than Mbp

$input = $ARGV[1];

if ($ARGV[2] =~ /^-g/) {
		$size = ($ARGV[3]*1000000);
		$input2 = "empty";
}
elsif ($ARGV[2] =~ /^-i/ && $ARGV[4] =~ /^-g/ ) {
		$input2 = $ARGV[3];
		$size = ($ARGV[5]*1000000);
}
else {
	print "Entry values undefined (error in multiple entry values) please adhere to: $USAGE\n";
	exit;
}
		
###############


print "The estimated genome size is: $size bp\n\n";


###############
# Send input files to subroutine to
# count nucleotides & collect results
###############

$total_nuc_count = (count_nucl ($input));

unless ($input2 eq "empty") {
	$total_nuc_count = $total_nuc_count + (count_nucl ($input2));
}

################


			

################
# Estimate genome coverage 
################

$est_cov = ($total_nuc_count/$size);
$est_cov = (sprintf "%.2f", $est_cov);

################



################
# Print results & exit
################

print "\n\nTotal results:\n";
print "\n There are a total of $total_nuc_count nucleotides in this file.\n";
print "\n This equates to an estimated genome coverage of $est_cov .\n";

exit;

################







####################################################
#	Main program - Subroutine for counting nucleotides
####################################################
#
#	Read each line from a fastq file and print it to an outfile. Every 4 lines ($linecount1),
#	increase the count of number of fastq accessions read by 1 ($linecount2). Also count
#	every line read and if the number of total number of lines read so far ($totalline)
#	If the number of lines read equals the desired number ($critvalue) then close the current
#	output file, open the second outfile and divert the remaining output to that file.

# Read each line from a fastq file and count every base on the second line. If this line
# is empty, add 1 to a count of blank lines ($blankline). Every 4 lines ($linecount1) 
# reset this value. Also count every line read ($totalline).
####################################################

sub count_nucl {

	my ($infile) = @_;
	my ($linecount1) = 0;
	my ($currentline) = "";
	my ($nuc_count) = 0;
	my ($totalline) = 0;
	my ($blankline) = 0;
	my ($line_length) = 0;


print "\nThe input file is: $infile\n";


unless ( open(SEQFILE, $infile) ) {
			print "cannot open the file \"$infile\" \n\n";
			exit;
	}

	while ($currentline = <SEQFILE>) {
				$linecount1++;
				if ($linecount1 == 2) {
						$line_length = ($currentline =~ tr/ATGCN//);
						$nuc_count = ($nuc_count + $line_length);
						if ($line_length == 0) {
								$blankline++;
								}
						}
				elsif ($linecount1 == 4) {
						$linecount1 = 0;
						$totalline++;
						}
	
				}

close (SEQFILE);

print "\nResults for: $infile\n";
print " Within this file of $nuc_count bp there were $totalline fastq sequences\n";
print " of these $blankline lines were empty.\n\n";

return $nuc_count;			
}

####################################################



