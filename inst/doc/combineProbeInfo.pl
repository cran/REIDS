#!/usr/bin/perl -w

#
# Based on Code Written by Mark Robinson, WEHI/Garvan
# mrobinson@wehi.edu.au
# Edited to make flat file, but minus the gene/group identification (i.e. no mps file). 
# - takes 1 argument which is the prefix, then takes
#   <prefix>.pgf, <prefix>.clf and creates
#   <prefix>.probeflat 
#Also now keeps all of the info in pgf and clf, and makes output with one row per probe (like a .bgf) . 
#Also create separate file of just the probeset names (unique)

use Getopt::Std;
use strict;

if (scalar(@ARGV) == 0){
 system("perldoc $0");
 exit(0);
}

my $PSR;
my $EXONPOS;
my $prefix = $ARGV[0];

my $PGFfile = "$prefix.pgf";
my $CLFfile = "$prefix.clf";
my $ADFfile = "$prefix.probeflat";
my $PSRfile = "$prefix.psr";
my $MISSfile = $prefix."_missingxy.probeflat";
my $line;
my %probes;
##print "Processing MPS $MPSfile\n";
##open(FILE,"<$MPSfile");
##while($line = <FILE>) {
##  if ($line =~ /^#/) { next; }
##  chop $line;
##  my @lineSplit = split("\t",$line);
##  my @psSplit = split(" ",$lineSplit[3]);
##  foreach (@psSplit) {
##    $probesets{$_}=$lineSplit[0];
##  }
##}
##close(FILE);

my $linecount=0;
my %psrCount;
my %psr;
my %psrid;
print "Processing CLF $CLFfile\n\tNumber of probes finished:\n";
open(FILE,"<$CLFfile");
while($line = <FILE>) {
  if ($line =~ /^#/) { next; }
  chop $line;
  my @lineSplit = split("\t",$line); #probeid \t x \t y
  $probes{$lineSplit[0]}=($lineSplit[1]) . "\t" . ($lineSplit[2]); #EP fixed here to get rid of +1
  $linecount++;
##  if($linecount % 100000 ==0) {print $linecount."...";}
##  if($linecount % 1000000 ==0){ print "\n";}

}
close(FILE);

print "Processing PGF $PGFfile, creating tab deliminated $ADFfile\n";
open(ADFFILE,">$ADFfile");
open(MISSFILE,">$MISSfile");
open(FILE,"<$PGFfile");

#so matches .bgp file
my $header="probeset_id\ttype\tprobeset_name\tatom_id\tprobe_id\ttype\tgc_count\tprobe_length\tinterrogation_position\tprobe_sequence";
print MISSFILE $header."\n";
print ADFFILE $header."\tx\ty\n";
$linecount=0;
while($line = <FILE>) {
  if ($line =~ /^#/) { 
    ##if($line=~/^#%header2(.+)/{ #grab the info here and make it header:
##        print ADFFILE $header.$1."\n";
##    }
    next; 
  }
  chop $line;
  my @lineSplit = split("\t",$line);
  my $probeInfo;
  my @psrLetters;
  my $psrKey;
  if ($lineSplit[0] eq "") {  # tab at start
    shift(@lineSplit); #remove starting tab
    if ($lineSplit[0] eq "") {  # two tabs -- last step
        shift(@lineSplit); #remove starting tab
        if(exists($probes{$lineSplit[0]})){
            $probeInfo=join("\t",@lineSplit);
            print ADFFILE $PSR."\t".$EXONPOS."\t".$probeInfo."\t".$probes{$lineSplit[0]}."\n";
        } else {
            print MISSFILE $PSR."\t".$EXONPOS."\t".$probeInfo."\n";
        }
    } else {                  # one tab -- second step
      $EXONPOS = join("\t",@lineSplit);
    }
  } else { # no tabs -- first step grab the PSR id for next iterations.
    if(!exists($psr{$lineSplit[2]})){ $psr{$lineSplit[2]}=$lineSplit[0];}
    @psrLetters=split("",$lineSplit[2]);
    $psrKey=join("",@psrLetters[0..2]);
    if(exists($psrCount{$psrKey})){ $psrCount{$psrKey}++;}
    else {$psrCount{$psrKey}=1;}
    $PSR = join("\t",@lineSplit);
  }
  $linecount++;
##  if($linecount % 100000 ==0) {print $linecount."...";}
##  if($linecount % 1000000 ==0){ print "\n";}

}
close(FILE);
close(ADFFILE);

print "\n";
foreach my $key (keys %psrCount){
    print $key."\t".$psrCount{$key}."\n";
    
}
open(PSRFILE,">$PSRfile");
print PSRFILE "probesetName\tPSID\n";
foreach my $key (keys %psr){
    print PSRFILE $key."\t".$psr{$key}."\n";    
}
