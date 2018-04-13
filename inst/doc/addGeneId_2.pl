#!/usr/bin/perl -w

#
# Given probe (per line) with probeset identification, add column with gene/group id
# - takes 3 arguments, the probe file, the gene file, the output file
# matches the pbsetid number in each column (column 1 in both files); if pbset not in genefile, then dropped.
# Returns .flat file

# use Getopt::Std;
# use strict;

if (scalar(@ARGV) == 0){
 system("perldoc $0");
 exit(0);
}

my $PSR;
my $EXONPOS;
my $probeFile = $ARGV[0];
my $geneFile=$ARGV[1];
my $ADFfile = $ARGV[2];
my $line;
my %probeset;
my $linecount=0;
print "Processing $geneFile\n";
open(GENEFILE,"<$geneFile");
while($line = <GENEFILE>) {
  chop $line;
  my @lineSplit = split("\t",$line); #PSID\tProbeSet\tTCID
  # print @lineSplit;
  # print $lineSplit[0];
  # print $lineSplit[1]; 
  # print $lineSplit[2];
  $probeset{$lineSplit[0]}=join("\t",$lineSplit[2]); 
}
close(GENEFILE);

 print "Processing $probeFile, creating tab deliminated .flat file\n";

open(ADFFILE,">$ADFfile");
open(FILE,"<$probeFile");
# so matches .bgp file
 my $header="TCID";
 $linecount=0;
 while($line = <FILE>) {
   chop $line; 
   if($linecount==0){print ADFFILE $line."\t".$header."\n";}
   my @lineSplit = split(/\t/,$line);  
   if($linecount==1){print $line."\n";}
   if(exists($probeset{$lineSplit[0]})){print ADFFILE $line."\t".$probeset{$lineSplit[0]}."\n";}
   $linecount++;
   if($linecount % 100000 ==0) {print $linecount."...";}
   if($linecount % 1000000 ==0){ print "\n";}
 }
close(FILE);
close(ADFFILE);



