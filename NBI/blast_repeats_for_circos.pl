#!/usr/bin/perl

use strict;
use warnings;
use Bio::SearchIO; 

#print usage statement if blast output filename is not provided on the command line
my $usage = "\nUSAGE: $0 blast_file min_length min_ID circos_chrom_name\n\n";
my $file = shift or die ($usage);
my $min_len = shift or die ($usage);
my $min_id = shift or die ($usage);
my $name = shift or die ($usage);

#Import Blast output file as a BioPerl object
my $SearchIO_obj = new Bio::SearchIO(-format => 'blast',-file => $file);
	
#loop through the blast file, going through each hsp and print output that can be read by Circos.	
while( my $result_obj = $SearchIO_obj->next_result ) {
	if( my $hit_obj = $result_obj->next_hit ) {
		my $count = 0;
		while (my $hsp_obj = $hit_obj->next_hsp){
			unless ($count){
				$count = 1;
				next;
			}
			if ($hsp_obj->percent_identity >= $min_id && $hsp_obj->length('total') >= $min_len){
				print "repeat$count\t$name\t",$hsp_obj->start('query'),"\t",$hsp_obj->start('query'),"\n";
				if ($hsp_obj->strand('hit') == 1){
					print "repeat$count\t$name\t",$hsp_obj->start('hit'),"\t",$hsp_obj->start('hit'),"\n";
				}else{
					print "repeat$count\t$name\t",$hsp_obj->end('hit'),"\t",$hsp_obj->end('hit'),"\n";
				}
				++$count;
			}
		}
	}
}
