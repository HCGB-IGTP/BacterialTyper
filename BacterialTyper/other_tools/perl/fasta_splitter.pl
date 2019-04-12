#!/usr/bin/perl
use strict;
use warnings;

my $file = $ARGV[0]; 
my $parts = $ARGV[1];
my $folder = $ARGV[2];

if (scalar @ARGV == 0) {
	print "Usage:\nperl $0 fasta_file parts\n\n";
	print "This script splits fasta in as many parts as stated...\n";
	exit();	
}
my $size = -s $file; ## get file size
print "Stats for file: $file\n";
print "Chars: $size\n";
my $block = int($size/$parts);

open (FH, "<$file") or die "Could not open source file. $!";
print "\t+ Splitting file into blocks of $block characters aprox ...\n";
my $j = 0; my @files;
while (1) {
	my $chunk;
	my $block_file = $folder."/tmp_part-".$j.".fna";
	push (@files, $block_file);
	open(OUT, ">$block_file") or die "Could not open destination file";
	if (!eof(FH)) { read(FH, $chunk,$block);  
		if ($j > 0) { $chunk = ">".$chunk; }
		print OUT $chunk;
	} ## Print the amount of chars	
	if (!eof(FH)) { $chunk = <FH>; print OUT $chunk; } ## print the whole line if it is broken	
	if (!eof(FH)) { 
		$/ = ">"; ## Telling perl where a new line starts
		$chunk = <FH>; chop $chunk; print OUT $chunk; 
		$/ = "\n";
	} ## print the sequence if it is broken
	$j++; close(OUT); last if eof(FH);
}
close(FH);
