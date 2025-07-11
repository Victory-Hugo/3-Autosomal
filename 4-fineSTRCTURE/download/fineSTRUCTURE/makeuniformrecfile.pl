#!/usr/bin/perl
#
use strict;
use warnings;

my $snplocs;
my @snps;

sub trim($){  # remove whitespace from beginning and end of the argument
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


if($#ARGV+1 !=2) {
print "Usage ./makeuniformrecfile.pl <phasefile> <outputfile>\n";
print "     <phasefile> is a valid chromopainter inputfile ending in .phase (in ChromoPainter v1 or v2 format) \n";
print("     <outputfile> will be a recombination file usable with <phasefile> in ChromoPainter, nominally in Morgans/base.\n");
print("The recombination rate is scaled to be approximately that in humans (0.1 Morgans/Mb). Because of this, it will NOT be usable directly and should only ever be used in conjunction with EM parameter estimation, which corrects for the global amount of recombination. If you are working on non-humans or simulated data, you may experience problems with EM estimation. The parameter may get stuck at a local mode where there is effectively infinite (or no) recombination. In this case, you should specify the initial conditions of ChromoPainter to have a much smaller or larger Ne (-n) value.\n");
die "\n";
}
my $phasefile = $ARGV[0];#'myres-1-1.phase';	
my $outputfile = $ARGV[1];#'test.recombfile';

## Set up the outputfile
open OUTPUTFILE, ">", $outputfile or die $!;
print OUTPUTFILE "start.pos recom.rate.perbp\n";

## Read in the SNPS
open PHASEFILE, $phasefile or die $!;
for (my $i = 0; $i < 4; ++$i){
  $snplocs=trim(<PHASEFILE>);
  if(substr($snplocs,0,1) eq "P") {last;}
}
close PHASEFILE;

@snps = split(/ /, $snplocs);
my $numsnps=scalar(@snps)-1;

print "Found $numsnps SNPS\n";
my $finalSNP=$snps[$numsnps];
## Loop over all SNPs to calculate their recombination position
for (my $i = 0; $i < $numsnps-1; ++$i){
  my $snploc=$snps[$i+1];
  my $overallrecrate = 1.0 /10000000;
  if($i<$numsnps-1){if($snps[$i+2]<$snps[$i+1]){
	$overallrecrate = -9;
  }}
  printf OUTPUTFILE "%i %.20f\n", $snploc, $overallrecrate;
}
printf OUTPUTFILE "%i %.20f\n", $finalSNP, 0;

close OUTPUTFILE;

print "Note: No absolute recombination rate.\nYou must perform EM estimation of Ne (e.g. -in -i 10) using chromopainter to use this recombination file\n";

