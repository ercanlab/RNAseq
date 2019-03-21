#! /usr/bin/perl -w
use strict;

#my @files = <*bam>;

#foreach my $file (@files){

	my $file = shift();
	my @temp = split(".bam", $file);
        my $id = $temp[0];


	my $output_file = $id."_counts.txt";
	#my $output_file = $id."_counts_newGTF.txt";
	my $sam_file = $id.".sam";
        print "Converting bam to sam file... ";
        my $sam = `samtools view -h -o $sam_file $file`;
    	print "Output file: $output_file\n";
	print "Running htseq-count issuing the command:";
	my $cmd = "htseq-count -s reverse $sam_file /scratch/cgsb/ercan/annot/c_elegans.WS220.annotations_new.gtf > $output_file";
	#my $cmd = "htseq-count -t exon -s reverse $sam_file /home/sea283/Annotations/Caenorhabditis_elegans.WBcel235.84.gtf > $output_file";
        print "$cmd\n";
	print "$cmd\n";
	my $start = time();
	my $res = `htseq-count -t exon -s reverse $sam_file /scratch/cgsb/ercan/annot/c_elegans.WS220.annotations_new.gtf > $output_file` ;
	#my $res = `htseq-count -t exon -s reverse $sam_file /home/sea283/Annotations/Caenorhabditis_elegans.WBcel235.84.gtf > $output_file` ;
	my $time = int((time() - $start)/60);
	print "Finished running htseq-count ($time min)\n\n\n";
	system("rm $sam_file");
#}

