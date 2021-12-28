#!/usr/bin/perl -w
use strict;
use Cwd;

# coded for Complete sweeps, window size of 5000bp.
# modified by Tiago Ribeiro (tribeiro@wisc.edu, tiaaagosr@gmail.com)

my $minValue=int($ARGV[0]); # flexible to define incomplete sweep boundaries, add random numbers for neutral ou complete sweeps
my $maxValue=int($ARGV[1]);

my $InputFile = $ARGV[2];

my $MinTract = 500;  #threshold bp length to analyze identical pairs of haplotypes

my $MinCount = 2;  #if fewer than this number of individuals carry the less common allele at a site, that polymorphism is ignored

my $neutral = 0;  #set to 1 for a neutral locus, 0 for selection
my $RequireComplete = 1;  #set to 1 to ensure complete sweeps only, 0 otherwise
my $InitialSampleSize1=0;
if($RequireComplete==1 && $neutral==0){
	$InitialSampleSize1 = 52;  #simulate population 1 with enough extra alleles to allow sampling of FinalSampleSize alleles carrying an adaptive allele
}else{
	$InitialSampleSize1=50;
}

my $FinalSampleSize = 50;  #the number of alleles to use in the analysis, for both populations
my $InitialSampleSize2 = 50;  #simulate extra pop2 alleles so we can select those without the adaptive allele
#Note - for neutral case, just simulate with final sample size (no extras)

my $IncSweep = 0;  #Set to 1 if you only want to consider incomplete sweeps with final adaptive allele frequencies within a certain range (otherwise set to 0)
my $MaxCount = $maxValue;  #The maximum number of sampled chromosomes in population 1 that carry the adaptive allele (if it's higher, skip this replicate).  Only matters if $IncSweep = 1.
my $MinAlleleCount = $minValue;  #The minimum number of sampled chromosomes in population 1 that carry the adaptive allele (if it's lower, skip this replicate).  Only matters if $IncSweep = 1.

my $XPEHHreps = 3;  #Only evaluate XPEHH for this many replicates

my $i = 0;
my $msmsOutputKeyword = "txt";  #only the msms output file should have this in the file name

#Declare some stuff

my $j = 0;
my $s = 0;
my $TotalPops = 0;
my $TotalSites = 0;
my $reps = 0;
my $median = 0;
my $denominator = 0;
my $NumberToTrim1 = $InitialSampleSize1 - $FinalSampleSize;
my $NumberToTrim2 = $InitialSampleSize2 - $FinalSampleSize;
my @InputAoA = ();
my @InputInd = ();
my @PopSizes = ();
my @LengthSumsA = ();
my @LengthSumsB = ();
my @denominators = ();
my @NonAdaptive = ();


#unpack the zip file containing msms and a few other files
my $pwd = cwd();


#READ THE INPUT DATA

open I, "<$InputFile"
  or die "Can't open input file \"$InputFile\" ($!)";

#get info from command line
$_ = (<I>);
chomp;
my @line = split;
if ($line[0] =~ m/msms/){
  shift @line;
}
my $TotalInds = $line[1];
my $TotalReps = $line[2];
for ($i = 2; $i < @line; $i++){
  if ($line[$i] =~ m/\-r/){
    $TotalSites = $line[$i+2];
    next;
  }
  if ($line[$i] =~ m/\-I/){
    $TotalPops = $line[$i+1];
    for ($j = 0; $j < $TotalPops; $j++){
      $s = $i + $j + 2;
      push @PopSizes, $line[$s];
    }
    last;
  }
}
print "From input file $InputFile, expecting data for $TotalInds individuals from $TotalPops populations, for $TotalReps replicates of a locus with $TotalSites sites.\n\nReplicates finished: ";

if (($PopSizes[0] != $InitialSampleSize1) || ($PopSizes[1] != $InitialSampleSize2)){
  die "Found simulated sample sizes $PopSizes[0] and $PopSizes[1].  Expected $InitialSampleSize1 and $FinalSampleSize.\n";
}
$PopSizes[0] = $FinalSampleSize;
$PopSizes[1] = $FinalSampleSize;

if ($TotalPops != 2){
  die "This program only works with two populations.  Found $TotalPops\n";
}
my @PopStarts = ();
my @PopStops = ();
my $start = 0;
my $stop = -1;
push @PopStarts, 0;
for ($i = 0; $i < @PopSizes; $i++){
  $stop = $stop + $PopSizes[$i];
  push @PopStops, $stop;
  last if ($i == (@PopSizes - 1));
  $start = $start + $PopSizes[$i];
  push @PopStarts, $start;
}
my @PopComps = ();
my $comps = 0;
for ($i = 0; $i < @PopSizes; $i++){
  $comps = 0;
  for ($j = $PopSizes[$i] - 1; $j > 0; $j--){
    $comps += $j;
  }
  push @PopComps, $comps;
}
my $BtwnComps = $PopSizes[0] * $PopSizes[1];

#Declare some more stuff
my $TargetSite = 0;
my $TargetLoc = 0.50000;
my $d = 0;
my $p = 0;
my $count = 0;
my $length = 0;
my $LengthSum = 0;
my $ratio = 0;
my $pi = 0;
my $diffs = 0;
my $pos = 0;
my $Dxy = 0;
my $FST = 0;
my $XPEHHmax = 0;
my $skip = 1;
my @positions = ();
my @DiffSites = ();
my @CHIASums = ();
my @CHIARatios = ();
my @MedMinRatios = ();
my @TargetFreqs = ();
my @pis = ();
my @FSTs = ();
my @XPEHHs = ();
my @AdaptiveAlleles = ();
my @NumberAdaptiveAlleles = ();

#Input loop
scalar (<I>);
while (<I>){
  chomp;
  $_ =~ s/\s+$//;
  next if (m%/%);
  next if (m/segsites/);
  if (m/positions/){
    @positions = split;
    shift @positions;
    if (@positions == 0){
      push @CHIARatios, 1;
    }
    next;
  }
  if ((m/0/) || (m/1/)){
    @line = split //, $_;
    push @InputAoA, [ @line ];

#Each time a replicate's data matrix is complete (all individuals present), calculate the CHIA ratios...
    if (@InputAoA == $TotalInds){

###
      for ($i = 1; $i < @InputAoA; $i++){
	if (@{$InputAoA[$i]} != @{$InputAoA[$i-1]}){
	  open E, ">error_output.txt";
	  for ($i = 0; $i < @positions; $i++){
	    print E "$positions[$i]\t";
	  }
	  print E "\n";
	  for ($i = 0; $i < @InputAoA; $i++){
	    for ($j = 0; $j < @{$InputAoA[$i]}; $j++){
	      print E "$InputAoA[$i][$j]\t";
	    }
	    print "\n";
	  }
	  die;
	  close E;
	}
      }
###
      
#Find the advantageous allele and report its frequency in population 
      if ($neutral == 0){
	@AdaptiveAlleles = ();
	$TargetSite = -1;
	for ($s = 0; $s < @positions; $s++){
	  if ($positions[$s] == $TargetLoc){
	    $TargetSite = $s;
	    last;
	  }
	}
	if ($TargetSite >= 0){
	  $count = 0;
	  @NonAdaptive = ();
	  for ($i = 0; $i < $InitialSampleSize1; $i++){
	    if ($InputAoA[$i][$s] ne '0'){
	      $count++;
	      push @AdaptiveAlleles, $InputAoA[$i][$s];
	      if ($InputAoA[$i][$s] ne '1'){
		$InputAoA[$i][$s] = 1
	      }
	    }
	    else{
	      push @NonAdaptive, $i;
	    }
	  }
	  push @TargetFreqs, $count;
	}
	else{
	  push @TargetFreqs, 0;
	  print "failed to find adaptive allele for a replicate\n";
	}
#Find the number of unique adaptive alleles in the population
	@AdaptiveAlleles = sort(@AdaptiveAlleles);
	for ($i = 1; $i < @AdaptiveAlleles; $i++){
	  if ($AdaptiveAlleles[$i] eq $AdaptiveAlleles[$i-1]){
	    splice @AdaptiveAlleles, $i, 1;
	    $i--;
	  }
	}
	$i = @AdaptiveAlleles;
	push @NumberAdaptiveAlleles, $i;

#If simulating incomplete sweeps with a specific range of final frequencies, skip replicates that fail this criterion
    if ($IncSweep == 1){
      if ( ($count > $MaxCount) || ($count < $MinAlleleCount) ){
        print "Skipped a replicate because $count individuals had an adaptive allele\n";
        @InputAoA = ();
        pop @TargetFreqs;
        pop @NumberAdaptiveAlleles;
        next;
      }
    }


#Filter replicates and individuals to obtain completely swept cases in pop 1 (if we're simulating complete sweeps)
	elsif ($count < ($FinalSampleSize * $RequireComplete)){
	  @InputAoA = ();
	  print "Skipped a replicate because not enough population 1 individuals had adaptive allele: $count\n";
	  pop @TargetFreqs;
	  pop @NumberAdaptiveAlleles;
	  next;
	}
	for ($i = 0; $i < $NumberToTrim1; $i++){
	  if ($i <= (@NonAdaptive - 1)){
	    $j = $NonAdaptive[$i];
	  }
	  else{
	    $j = $FinalSampleSize;
	  }
	  splice @InputAoA, $j, 1;
	}

      }


    
#evalute CHI ratio
#first remove singletons from InputAoA and positions
      for ($s = 0; $s < @positions; $s++){
	$count = 0;
	for ($i = 0; $i < ($FinalSampleSize * 2); $i++){
	  $count += $InputAoA[$i][$s];
	}
	if ( ($count < ($MinCount - 0.5)) || ($count > (@InputAoA - ($MinCount-0.5))) ){
	  splice @positions, $s, 1;
	  for ($i = 0; $i < @InputAoA; $i++){
	    splice @{$InputAoA[$i]}, $s, 1;
	  }
	  $s--;
	}
      }

#For each population, cycle through each pair of sequences
      @CHIASums = ();
      for ($p = 0; $p < @PopSizes; $p++){
	$LengthSum = 0;
	for ($i = $PopStarts[$p]; $i < $PopStops[$p]; $i++){
	  for ($j = $i + 1; $j <= $PopStops[$p]; $j++){

#Record sites that differ between this pair of sequences
	    @DiffSites = ();
	    push @DiffSites, 0;
	    for ($s = 0; $s < @positions; $s++){
	      if ($InputAoA[$i][$s] != $InputAoA[$j][$s]){
		push @DiffSites, $positions[$s];
	      }
	    }
	    push @DiffSites, 1;

#Identify tracts of identical sequence longer than the threshold
	    for ($d = 0; $d < @DiffSites - 1; $d++){
	      $length = (($DiffSites[$d+1] - $DiffSites[$d]) * $TotalSites) - 1;
	      if ($length >= $MinTract){
		$LengthSum += $length;
	      }
	    }
	  }
	}
	if ($p == 0){
	  push @LengthSumsA, $LengthSum;
	}
	else{
	  push @LengthSumsB, $LengthSum;
	}
      }
      if ($LengthSumsB[-1] > 0){
	$ratio = ($LengthSumsA[-1] / $PopComps[0]) / ($LengthSumsB[-1] / $PopComps[1]);
      }
      else{
	$ratio = ($LengthSumsA[-1] / $PopComps[0]) / (($MinTract / 2) / $PopComps[1]);
      }
      @InputAoA = ();
      $reps++;
      print "replicate $reps, Ben. Allele Freq. $TargetFreqs[-1], CHIA Ratio $ratio.\n";
    }
  }
}

#evaluate simple ratios
#if ($UseMedianDenom == 0){
  for ($i = 0; $i < @LengthSumsA; $i++){
    if ($LengthSumsB[$i] > 0){
      $ratio = ($LengthSumsA[$i] / $PopComps[0]) / ($LengthSumsB[$i] / $PopComps[1]);
    }
    else{
      $ratio = ($LengthSumsA[$i] / $PopComps[0]) / (($MinTract / 2) / $PopComps[1]);
    }
    push @CHIARatios, $ratio;
  }
#}

#evaluate ratios using median denominator as minimum denominator 
#else{
  for ($i = 0; $i < @LengthSumsB; $i++){
    $denominator = $LengthSumsB[$i] / $PopComps[1];
    push @denominators, $denominator;
  }
  @denominators = sort { $a <=> $b } @denominators;
  $i = @denominators;
  if ( ($i % 2) == 1){
    $i = (@denominators - 1) / 2;
    $median = $denominators[$i];
  }
  else{
    $i = @denominators / 2;
    $median = ($denominators[$i] + $denominators[$i-1]) / 2;
  }
  if ($median == 0){
    die "median denominator was zero.  try changing your parameters\n";
  }
  for ($i = 0; $i < @LengthSumsA; $i++){
    $denominator = $LengthSumsB[$i] / $PopComps[1];
    if ($denominator > $median){
      $ratio = ($LengthSumsA[$i] / $PopComps[0]) / $denominator;
      push @MedMinRatios, $ratio;
    }
    else{
      $ratio = ($LengthSumsA[$i] / $PopComps[0]) / $median;
      push @MedMinRatios, $ratio;
    }
  }
#}

#@CHIARatios = sort { $a <=> $b } @CHIARatios;

#output to file
chdir ("$pwd");
my $file = 'idhap' . '_' . $InputFile;
open O, ">$file";
for ($i = 0; $i < @CHIARatios; $i++){
  print O "$TargetFreqs[$i]\t$NumberAdaptiveAlleles[$i]\t$CHIARatios[$i]\t$MedMinRatios[$i]\t";

  print O "\n";
}
close O;

chdir("$pwd");

