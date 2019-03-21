#! /usr/bin/perl -w
use strict;

#my @files = <*bam>;

#foreach my $file (@files){

	my $file = shift();
	
	my @temp = split("_acc", $file);
	my $id = $temp[0];

	my $output_dir = $id."_cufflinks";
    	print "Output directory: $output_dir\n";
	print "Running cufflinks issuing the command:";
	my $cmd = "cufflinks -p 8 --library-type fr-firststrand -G /scratch/cgsb/ercan/annot/c_elegans.WS220.annotations_forCufflinks.gff3 -o $output_dir $file";
	print "$cmd\n";
	my $start = time();
	my $res = `cufflinks -p 8 --library-type fr-firststrand -G /scratch/cgsb/ercan/annot/c_elegans.WS220.annotations_forCufflinks.gff3 -o $output_dir $file`;
	my $time = int((time() - $start)/60);
	print "Finished running CUFFLINKS ($time min)\n\n\n";
#}

