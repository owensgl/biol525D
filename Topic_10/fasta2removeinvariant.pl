#!/bin/perl
use warnings;
use strict;
#This script removes invariant sites from a fasta file. 
my $sample_counter = 0;
my %data_hash;
my %name_hash;
my %base_hash;
my $n_bases;
while (<STDIN>){
  chomp;
  my $line = $_;
  if ($line =~ /^\>/){
    $sample_counter++;
    $name_hash{$sample_counter} = $line;
  }else{
    my @bases=split(//,$line);
    $n_bases = $#bases;
    foreach my $n (0..$#bases){
      $data_hash{$sample_counter}{$n} = $bases[$n];
      if ($bases[$n] ne "N"){
        $base_hash{$n}{$bases[$n]}++;
      }
    }
  }
}
#For the first sample_counter
print "$name_hash{1}\n";
foreach my $i (0..$n_bases){
  my $bases = keys %{ $base_hash{$i} };
  if ($bases > 1){
    print "$data_hash{1}{$i}";
  }
}

#For the rest
foreach my $n (2..$sample_counter){
  print "\n$name_hash{$n}\n";
  foreach my $i (0..$n_bases){
    my $bases = keys %{ $base_hash{$i} };
    if ($bases > 1){
      print "$data_hash{$n}{$i}";
    }
  }
}
