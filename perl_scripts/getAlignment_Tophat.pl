#! /usr/bin/perl -w
use strict;

#my @files = <*fastq*>;

#foreach my $file (@files){

my $file = shift();

	my @temp = split(/\.fastq/, $file);
	my $id = $temp[0];

	my $log2 = "Read_Alignment_module_log_".$id.".txt";
	open(LOG2,">", $log2);
	print LOG2 "Input Sequenced Reads .fastq file: $file\n\n\n";
	       
	if ($file =~ /gz/) {
		my $cmd = "gzip -d $file";
		chomp ($cmd);
		system($cmd);
		$file =~ s/\.gz//g;
		print LOG2 "Unzipped file: $file\n";	
	}
	       
	my $indexfile = "/scratch/cgsb/ercan/annot/forBowtie/c_elegans.WS220"; #p_pacificus.WS228"; #change me 
	print LOG2 "Index file is: $indexfile\n";
    	my $output_dir = $id."_ce10_tophat";
    	print LOG2 "Output directory: $output_dir\n";
	print LOG2 "Running TOPHAT issuing the command:";
	my $cmd = "tophat -p 8 --library-type fr-firststrand -o $output_dir $indexfile $file";
	print LOG2 "$cmd\n";
	my $start = time();
	my $tophat = `tophat -p 8 --library-type fr-firststrand -o $output_dir $indexfile $file`;
	my $time = int((time() - $start)/60);
	print LOG2 "Finished running TOPHAT ($time min)\n\n\n";
	print LOG2 "Counting reads in fastq file:\n";
	my $n_fastq = `wc -l $file`;
    	$n_fastq =~ s/\s+//g;
    	$n_fastq =~ s/$file//g;
    	$n_fastq /= 4;
    	print LOG2 "$n_fastq\n";
    	print LOG2 "Counting aligned reads in bam file:\n";
    	my $n_bam = `samtools view $output_dir/accepted_hits.bam | cut -f 1 | sort | uniq | wc -l`;
    	print LOG2 "$n_bam\n";
    	my $percent = ($n_bam / $n_fastq) * 100;
    	print LOG2 sprintf("%.2f", $percent)."% of reads aligned.\n";
	rename("$output_dir/accepted_hits.bam","$output_dir/".$id."_accepted_hits.bam");
	rename("$output_dir/unmapped.bam","$output_dir/".$id."_unmapped.bam");
	close LOG2;
#}

