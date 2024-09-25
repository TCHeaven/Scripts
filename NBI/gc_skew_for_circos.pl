#!/usr/bin/perl 

#A program to analyze GC skew in a (circular) genome with a given window and step size.
#I have modified this to produce output for circos

use warnings;
use strict;
use POSIX;

my $usage = "\nUSAGE: $0 fastaFile windowSize stepSize chromosomeName[forCircos] posColor[forCircos] negColor[forCircos]\n\n";

my $file = shift or die ($usage);
my $windowSize = shift or die ($usage);
my $stepSize = shift or die ($usage);
my $chromosomeName = shift or die ($usage);
my $posColor = shift or die ($usage);
my $negColor = shift or die ($usage);
my $leftSize = ceil(($windowSize - 1)/2);
my $rightSize = floor (($windowSize - 1)/2);


my ($headerRef, $seqRef) = get_fasta_names_and_seqs($file);
my @seqArray = @$seqRef;

my $seq = $seqArray[0];
my $length = length ($seq);

for (my $i = 0; $i < $length; $i += $stepSize){
	my $gcount = 0;
	my $ccount = 0;
	
	my $left = $i - $leftSize;
	my $right = $i + $rightSize;
	
	
	if ($left < 0){
		$left = $length + $left;
	} 
	if ($right >= $length){
		$right = $right - $length;
	}
	

	if ($left < $right){
		my $substr = uc(substr ($seq, $left, $right - $left + 1));
		$ccount += $substr =~ tr/C/C/;
		$gcount += $substr =~ tr/G/G/;		
	}else{
		my $substr1 = uc(substr ($seq, 0, $right + 1));
		my $substr2 = uc(substr ($seq, $left));
		$ccount += $substr1 =~ tr/C/C/;
		$gcount += $substr1 =~ tr/G/G/;		
		$ccount += $substr2 =~ tr/C/C/;
		$gcount += $substr2 =~ tr/G/G/;		
	}	
	print "$chromosomeName\t", $i+1, "\t", $i+1, "\t", ($gcount - $ccount)/($gcount + $ccount), "\tfill_color=";
	if (($gcount - $ccount)/($gcount + $ccount) > 0){
		print $posColor, "\n";
	}else{
		print $negColor, "\n";
	}

}

exit;



sub get_fasta_names_and_seqs {
	use strict;
	use warnings;

	my ($inputfilename) = @_;
	my @fasta_names = ();
	my @fasta_seqs= ();

		   
	unless ( open(FILEDATA, $inputfilename) ) {
		print STDERR "Cannot open file \"$inputfilename\"\n\n"; #print error message
		exit; #exit the program
	}	

	my @filedata = <FILEDATA>; #Read the lines of the file into an array
	close FILEDATA;
	
	my $seq_count = 0; #this will be used to keep track of the number of sequences
	foreach my $line (@filedata){
		chomp $line;
		if ($line =~ /^\s*$/) {next;} #ignore line if it is blank
		
		elsif ($line =~ /^>/) { #if the line is a header line (begins with ">")...
			if ($line =~ /^>.*[\w]+/){
				my $partialLine = substr ($&, 1);
				push (@fasta_names, $partialLine); #add that line to an array of fasta names
				push (@fasta_seqs, ""); #and add a new blank element to an array of sequences
				++$seq_count; #also increment our counter which keeps track of sequence number
			}
		}	
		
		else { #if the line's not blank or a header, add it to the current sequence 
			$fasta_seqs[$seq_count-1] .= $line;
		}
		
		$fasta_seqs[$seq_count-1] =~s/\s//g; #remove all whitespace from the current  sequence
	}
	
	return (\@fasta_names, \@fasta_seqs);

}
