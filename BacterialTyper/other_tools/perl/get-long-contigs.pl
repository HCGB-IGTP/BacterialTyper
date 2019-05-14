#!/usr/bin/perl
use strict;
use warnings;

## get contigs longer than 5kb
my $file = $ARGV[0];
my $len = $ARGV[1];

if (!$file || !$len) {
	print "No inputs provided\nPlease provide a contig file and a size\n";
	print "perl $0 contig_file.fasta length\n";
	print "This script retrieves sequences longer than length stated (len_seq > user_len)\n";
	exit();
}

my $user_len = int($len);
my %hash;
open(FILE, $file) || die "Could not open the $file ...\n";
$/ = ">"; ## Telling perl where a new line starts
while (<FILE>) {		
	next if /^#/ || /^\s*$/;
	chomp;
	my ($titleline, $sequence) = split("\n",$_,2);
	next unless ($sequence && $titleline);
    $sequence =~ s/\n//g;
	my $len_seq = length($sequence);
	if ($len_seq > $user_len) {
		print ">$titleline\n$sequence\n";
}}
close(FILE);
$/ = "\n";
