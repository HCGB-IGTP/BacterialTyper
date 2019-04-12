#!/usr/bin/perl
use strict;
use warnings;

my $ids = $ARGV[0];
my $file = $ARGV[1];

if (!@ARGV) {
	print "Usage: perl $0 ids fasta\n";
	exit();
}

## get ids
open (ID,"$ids");
my @IDS = (<ID>);
close (ID);

chomp @IDS;

#read proteins
my %hash;
open(FILE, $file) || die "Could not open the $file ...\n";
$/ = ">"; ## Telling perl where a new line starts
while (<FILE>) {
	next if /^#/ || /^\s*$/;
        chomp;
        my ($titleline, $sequence) = split(/\n/,$_,2);
        next unless ($sequence && $titleline);
	my @id = split(" ", $titleline); 
        chop $sequence;
	$sequence =~ s/\n//g;
	if (grep /$id[0]$/, @IDS) {
	        print ">".$titleline."\n".$sequence."\n";
	}
}
close(FILE); $/ = "\n";

