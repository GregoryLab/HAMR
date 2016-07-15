#! /usr/bin/perl
use warnings;
use strict;

###version 2

### parses coverageBed -s -d output (per base coverage) to generate per-base histogram (R is super-slow with these operations) Assumes integer data!

unless (scalar @ARGV == 0) {
print "usage: perl coverageBed_depth_histogram.pl <coverageBed_depth_input.bed >output_file\n";
exit;		
}

#my $coverageBed_file = shift or die;

# ----------------------------------------------
# generate histogram from each line of coverageBed output
# ----------------------------------------------

##input example
# chr start stop transcript id score strand relative_pos coverage
#chr1    16774864        16777247        AT1G44110.1     #AT1G44110;lots_of_exons        -       1       0
my ($chr, $start, $stop, $transcript, $id, $score, $strand, $relative_pos, $coverage);
my %hist;

#open (COVBED, "<$coverageBed_file") or die "Unable to initialize file $coverageBed_file for writing: $!";
while (my $line = <STDIN>) {
	chomp $line;
	my @line = split(/\t/, $line);
	$coverage = pop(@line);
    $hist{$coverage}++; #actual binning step.
}


print STDOUT "$_\t$hist{$_}\n" for (sort {$a <=> $b} keys %hist);
