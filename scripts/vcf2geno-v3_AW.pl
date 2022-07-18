#!/usr/bin/perl
use strict;
use warnings;

use List::Util qw( min max sum );

#Script to take Sanger's VCF files and make calls based off of the read calls. Assumes that all vcf filters have been applied as desired. Also adds in the maf at the position for that population

## test file: test4-5000

## v2 - change requiring the dominant allele to have 10 reads to requiring a total of 10 reads required, regardless of alleles

## v3 takes into account the listed presence of more than 2 possible alleles in the data even if only 2 alleles are present

open(IN, "$ARGV[0]");
open(OUTALL, ">$ARGV[0]-maj-all-v3");
open(OUTMISS, ">$ARGV[0]-maj-invar-v3");
open(OUTKEEP, ">$ARGV[0]-maj-keep-v3");
my @line;
my $ref;
my $alt;
my $reffreq;
my $altfreq;
my $missing;
my @outline;

while(<IN>){
  @outline=();
  $reffreq = 0;
  $altfreq = 0;
  $missing = 0;
  chomp;
  next if ($_ =~m/##/);
  if ($_ =~m/#/ ){
      my $header="${_}\tnumref\tnumalt\tmissing\n";
      $header=~ s/\n/ /g;
       $header=~ s/INFO//g;
      $header=~ s/FORMAT//g;
     $header=~ s/(\s)+/$1/g;
    push(@outline,"$header\n");
    print OUTALL join("\t",@outline);
    print OUTMISS join("\t",@outline);
    print OUTKEEP join("\t",@outline);
    next;}
  my @line = split(/\t/, $_);
  $ref = $line[3];
  $alt = $line[4];
  my @popalts=split(/,/,$alt);
  my $linelength = @line;
  for (my $i = 0; $i < $linelength; $i++){
    if ($i < 9){
        if($i < 7){push(@outline, "$line[$i]");}
    }
    else {
      my @geno = split(/:/,$line[$i]);
        if(defined($geno[1])){
            $geno[1]=~ s/\./0/g;
            my @depth = split(/,/,$geno[1]);
            if(!defined($depth[0])){$depth[0]=0;}
            if(!defined($depth[1])){$depth[1]=0;}
            if(sum(@depth)<5){push(@outline, "N"); $missing++;}
            elsif ($depth[0] == $depth[1]){push(@outline, "1"); $altfreq++;}
            elsif ($depth[0]< $depth[1]){push(@outline, "1"); $altfreq++;}
            elsif ($depth[0]>$depth[1]){push(@outline, "0"); $reffreq++;}
            else {push(@outline, "N");$missing++;}
        }
    }
  }
  push(@outline, "$reffreq\t$altfreq\t$missing\n");
  print OUTALL join("\t",@outline);
  if ($reffreq == 0 or $altfreq == 0){print OUTMISS join("\t",@outline);}
  else {print OUTKEEP join("\t",@outline);}; #getc;
}


