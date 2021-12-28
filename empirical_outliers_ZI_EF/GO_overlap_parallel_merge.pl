#!/usr/bin/perl -w
use strict;

#merges output files from GO_overlap_parallel into one output files



#my $InputFileKey = 'KF_new_PBE_regions_01_GO';
my $InputFileKey = 'EFht_ZI_winFST_regions_01_GO';
my $OutputFile = $InputFileKey . '_merged.txt';

my $i = 0;
my $j = 0;
my $PValue = 0;
my @line = ();


#get sim peaks from input file(s)
my @OutputAoA = ();
my @AllFiles = ();
my @InputFiles = ();
opendir DIR, "." or die "couldn't open directory\n";
@AllFiles = readdir(DIR);
closedir DIR;
for ($i = 0; $i < @AllFiles; $i++){
  if (($AllFiles[$i] =~ m/$InputFileKey/) && ($AllFiles[$i] =~ m/txt/)){
    unless ($AllFiles[$i] =~ m/merged/){
      push @InputFiles, $AllFiles[$i];
    }
  }
}
$i = @InputFiles;
print "Found $i input files\n";

open I, "<$InputFiles[0]" or die;
while (<I>){
  chomp;
  last if m/^$/;
  @line = split("\t",$_);
  push @OutputAoA, [ @line ];
}
close I;

$i = @OutputAoA;
print "Found $i GO categories in first permutation file\n";


for ($i = 1; $i < @InputFiles; $i++){
  open I, "<$InputFiles[$i]" or die;
  $j = 0;
  while (<I>){
    chomp;
    last if m/^$/;
    @line = split("\t",$_);
    $PValue = $line[5];
    $OutputAoA[$j][5] = $OutputAoA[$j][5] + $PValue;
    $j++;
  }
  close I;
}

for ($i = 0; $i < @OutputAoA; $i++){
  $OutputAoA[$i][5] = $OutputAoA[$i][5] / @InputFiles;
}

open O, ">$OutputFile";
for ($i = 0; $i < @OutputAoA; $i++){
  for ($j = 0; $j < @{$OutputAoA[$i]}; $j++){
    print O $OutputAoA[$i][$j];
    if ($j == (@{$OutputAoA[$i]} - 1)){
      print O "\n";
    }
    else{
      print O "\t";
    }
  }
}
close O;


