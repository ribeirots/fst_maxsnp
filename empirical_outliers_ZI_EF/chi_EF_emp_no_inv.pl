#!/usr/bin/perl -w
#use strict;

#Calculate the ratio of identical haplotype sharing (excluding singleton polymorphisms) in previously defined windows

my @chrs = ('ChrX');

my $WinProp = 0.1; #set this to determine the proportion of a window that must be identical between two individuals in order to count that haplotype identity block

my $WinFilePrefix = 'windows_ZI250_';
my $OutputFile = 'chi_EF_ZI.txt';

my @pops = ('EF','ZI');

my @InLines1 = ("EF101N", "EF102N", "EF103N", "EF10N", "EF110N", "EF112N", "EF114N", "EF115N", "EF116N", "EF117N", "EF119N", "EF11N", "EF120N", "EF122N", "EF126N", "EF12N", "EF130N", "EF131N", "EF132N", "EF134N", "EF135N", "EF136N", "EF15N", "EF16N", "EF19N", "EF1N", "EF24N", "EF25N", "EF26N", "EF27N", "EF2", "EF31N", "EF32N", "EF35N", "EF39", "EF3N", "EF43N", "EF46N", "EF48N", "EF52N", "EF53N", "EF54N", "EF59N", "EF64N", "EF65", "EF66N", "EF67N", "EF69N", "EF6N", "EF70N", "EF73N", "EF74N", "EF75N", "EF78", "EF7N", "EF80N", "EF81N", "EF82N", "EF83N", "EF86N", "EF88N", "EF8N", "EF93N", "EF95N", "EF96N", "EF98N", "EF9N");
my @InLines2 = ("ZI10", "ZI103", "ZI104", "ZI114N", "ZI117", "ZI118N", "ZI126", "ZI129", "ZI134N", "ZI136", "ZI138", "ZI152", "ZI157", "ZI161", "ZI164", "ZI165", "ZI167", "ZI170", "ZI172", "ZI173", "ZI176", "ZI177", "ZI178", "ZI179", "ZI181", "ZI182", "ZI184", "ZI185", "ZI188", "ZI190", "ZI191", "ZI193", "ZI194", "ZI196", "ZI197N", "ZI198", "ZI199", "ZI200", "ZI202", "ZI206", "ZI207", "ZI210", "ZI211", "ZI212", "ZI213", "ZI214", "ZI216N", "ZI218", "ZI219", "ZI220", "ZI221", "ZI225", "ZI226", "ZI227", "ZI228", "ZI230", "ZI231", "ZI232", "ZI233", "ZI235", "ZI237", "ZI239", "ZI240", "ZI241", "ZI250", "ZI251N", "ZI252", "ZI253", "ZI254N", "ZI255", "ZI26", "ZI261", "ZI263", "ZI264", "ZI265", "ZI267", "ZI268", "ZI269", "ZI27", "ZI271", "ZI273N", "ZI276", "ZI279", "ZI281", "ZI284", "ZI286", "ZI292", "ZI293", "ZI295", "ZI296", "ZI303", "ZI311N", "ZI313", "ZI314", "ZI316", "ZI317", "ZI319", "ZI31N", "ZI320", "ZI321", "ZI324", "ZI329", "ZI33", "ZI332", "ZI333", "ZI335", "ZI336", "ZI339", "ZI341", "ZI342", "ZI344", "ZI348", "ZI351", "ZI352", "ZI353", "ZI357N", "ZI358", "ZI359", "ZI362", "ZI364", "ZI365", "ZI368", "ZI370", "ZI373", "ZI374", "ZI377", "ZI378", "ZI379", "ZI380", "ZI381", "ZI384", "ZI386", "ZI388", "ZI392", "ZI394N", "ZI395", "ZI396", "ZI397N", "ZI398", "ZI400", "ZI402", "ZI405", "ZI413", "ZI418N", "ZI420", "ZI421", "ZI429", "ZI431", "ZI433", "ZI437", "ZI443", "ZI444", "ZI445", "ZI446", "ZI448", "ZI453", "ZI455N", "ZI456", "ZI457", "ZI458", "ZI460", "ZI466", "ZI467", "ZI468", "ZI471", "ZI472", "ZI476", "ZI477", "ZI486", "ZI488", "ZI490", "ZI491", "ZI504", "ZI505", "ZI508", "ZI50N", "ZI514N", "ZI517", "ZI523", "ZI524", "ZI527", "ZI530", "ZI56", "ZI59", "ZI61", "ZI68", "ZI76", "ZI81", "ZI85", "ZI86", "ZI90", "ZI91", "ZI99");


# Chr2L
# my @InLines1 = ('EF101N', 'EF102N', 'EF103N', 'EF108N', 'EF10N', 'EF110N', 'EF112N', 'EF114N', 'EF115N', 'EF116N', 'EF117N', 'EF119N', 'EF11N', 'EF120', 'EF120N', 'EF122N', 'EF126N', 'EF12N', 'EF130N', 'EF131N', 'EF132N', 'EF134N', 'EF135N', 'EF136N', 'EF15N', 'EF16N', 'EF19N', 'EF1N', 'EF2', 'EF24N', 'EF25N', 'EF26N', 'EF27N', 'EF31N', 'EF32N', 'EF35N', 'EF39', 'EF3N', 'EF43N', 'EF46N', 'EF48N', 'EF52N', 'EF53N', 'EF54N', 'EF59N',  'EF64N', 'EF65', 'EF66N', 'EF67N', 'EF69N', 'EF6N', 'EF70N', 'EF73N', 'EF74N', 'EF75N', 'EF78', 'EF7N', 'EF80N', 'EF81N', 'EF82N', 'EF83N', 'EF86N', 'EF88N', 'EF8N', 'EF93N', 'EF95N', 'EF96N', 'EF98N');
# my @InLines2 = ('ZI10', 'ZI103', 'ZI104', 'ZI126', 'ZI129', 'ZI134N', 'ZI138', 'ZI152', 'ZI157', 'ZI164', 'ZI165', 'ZI167', 'ZI170', 'ZI172', 'ZI173', 'ZI176', 'ZI177', 'ZI178', 'ZI179', 'ZI181', 'ZI182', 'ZI184', 'ZI185', 'ZI188', 'ZI190', 'ZI191', 'ZI193', 'ZI194', 'ZI196', 'ZI197N', 'ZI199', 'ZI200', 'ZI202', 'ZI207', 'ZI210', 'ZI211', 'ZI212', 'ZI213', 'ZI214', 'ZI216N', 'ZI218', 'ZI219', 'ZI220', 'ZI221', 'ZI225', 'ZI226', 'ZI228', 'ZI230', 'ZI231', 'ZI233', 'ZI237', 'ZI239', 'ZI241', 'ZI250', 'ZI251N', 'ZI252', 'ZI254N', 'ZI255', 'ZI26', 'ZI261', 'ZI264', 'ZI265', 'ZI267', 'ZI268', 'ZI271', 'ZI273N', 'ZI279', 'ZI28', 'ZI281', 'ZI286', 'ZI291', 'ZI292', 'ZI293', 'ZI295', 'ZI296', 'ZI303', 'ZI311N', 'ZI314', 'ZI317', 'ZI319', 'ZI31N', 'ZI320', 'ZI324', 'ZI329', 'ZI33', 'ZI332', 'ZI333', 'ZI335', 'ZI336', 'ZI341', 'ZI344', 'ZI348', 'ZI351', 'ZI353', 'ZI358', 'ZI359', 'ZI362', 'ZI364', 'ZI365', 'ZI368', 'ZI370', 'ZI373', 'ZI377', 'ZI378', 'ZI379', 'ZI380', 'ZI381', 'ZI382', 'ZI384', 'ZI386', 'ZI388', 'ZI392', 'ZI398', 'ZI400', 'ZI402', 'ZI413', 'ZI418N', 'ZI420', 'ZI421', 'ZI429', 'ZI431', 'ZI433', 'ZI437', 'ZI443', 'ZI447', 'ZI448', 'ZI453', 'ZI455N', 'ZI457', 'ZI458', 'ZI460', 'ZI466', 'ZI468', 'ZI471', 'ZI472', 'ZI476', 'ZI477', 'ZI490', 'ZI491', 'ZI504', 'ZI508', 'ZI50N', 'ZI514N', 'ZI517', 'ZI523', 'ZI56', 'ZI59', 'ZI68', 'ZI85', 'ZI86', 'ZI90', 'ZI91', 'ZI99');

# Chr2R
# my @InLines1 = ("EF101N", "EF102N", "EF103N", "EF108N", "EF10N", "EF110N", "EF112N", "EF114N", "EF115N", "EF116N", "EF117N", "EF119N", "EF11N", "EF120N", "EF122N", "EF126N", "EF12N", "EF130N", "EF131N", "EF132N", "EF134N", "EF135N", "EF136N", "EF15N", "EF16N", "EF19N", "EF1N", "EF24N", "EF25N", "EF26N", "EF27N", "EF2", "EF31N", "EF32N", "EF35N", "EF39", "EF3N", "EF43N", "EF46N", "EF48N", "EF52N", "EF53N", "EF54N", "EF59N", "EF64N", "EF65", "EF66N", "EF67N", "EF69N", "EF6N", "EF70N", "EF73N", "EF74N", "EF75N", "EF78", "EF7N", "EF80N", "EF81N", "EF82N", "EF83N", "EF86N", "EF88N", "EF8N", "EF93N", "EF95N", "EF96N", "EF98N", "EF9N");
# my @InLines2 = ("ZI103", "ZI104", "ZI114N", "ZI117", "ZI118N", "ZI126", "ZI129", "ZI134N", "ZI138", "ZI152", "ZI157", "ZI161", "ZI164", "ZI165", "ZI167", "ZI170", "ZI172", "ZI173", "ZI176", "ZI177", "ZI178", "ZI181", "ZI184", "ZI185", "ZI188", "ZI190", "ZI191", "ZI193", "ZI194", "ZI196", "ZI197N", "ZI198", "ZI199", "ZI200", "ZI202", "ZI206", "ZI207", "ZI210", "ZI211", "ZI212", "ZI214", "ZI218", "ZI219", "ZI220", "ZI221", "ZI225", "ZI227", "ZI228", "ZI230", "ZI231", "ZI232", "ZI233", "ZI239", "ZI240", "ZI241", "ZI250", "ZI251N", "ZI252", "ZI253", "ZI254N", "ZI26", "ZI261", "ZI263", "ZI264", "ZI265", "ZI267", "ZI268", "ZI271", "ZI273N", "ZI279", "ZI28", "ZI281", "ZI284", "ZI286", "ZI291", "ZI292", "ZI293", "ZI296", "ZI303", "ZI311N", "ZI313", "ZI314", "ZI316", "ZI317", "ZI319", "ZI31N", "ZI320", "ZI324", "ZI329", "ZI33", "ZI332", "ZI333", "ZI335", "ZI336", "ZI339", "ZI341", "ZI342", "ZI344", "ZI352", "ZI359", "ZI362", "ZI368", "ZI370", "ZI373", "ZI374", "ZI377", "ZI378", "ZI379", "ZI380", "ZI381", "ZI382", "ZI386", "ZI392", "ZI394N", "ZI395", "ZI396", "ZI397N", "ZI398", "ZI400", "ZI402", "ZI405", "ZI413", "ZI418N", "ZI420", "ZI421", "ZI429", "ZI431", "ZI433", "ZI437", "ZI443", "ZI444", "ZI445", "ZI446", "ZI447", "ZI453", "ZI455N", "ZI456", "ZI457", "ZI458", "ZI460", "ZI468", "ZI471", "ZI472", "ZI476", "ZI477", "ZI486", "ZI488", "ZI490", "ZI491", "ZI504", "ZI505", "ZI508", "ZI50N", "ZI514N", "ZI517", "ZI523", "ZI59", "ZI61", "ZI68", "ZI76", "ZI81", "ZI85", "ZI86", "ZI90", "ZI91", "ZI99");

# Chr3L
# my @InLines1 = ("EF101N", "EF102N", "EF103N", "EF108N", "EF10N", "EF110N", "EF114N", "EF115N", "EF116N", "EF117N", "EF119N", "EF11N", "EF120N", "EF122N", "EF126N", "EF12N", "EF130N", "EF131N", "EF132N", "EF135N", "EF136N", "EF15N", "EF16N", "EF19N", "EF1N", "EF24N", "EF25N", "EF26N", "EF27N", "EF2", "EF31N", "EF32N", "EF35N", "EF39", "EF3N", "EF43N", "EF46N", "EF48N", "EF52N", "EF53N", "EF54N", "EF59N", "EF64N", "EF65", "EF66N", "EF67N", "EF69N", "EF6N", "EF70N", "EF73N", "EF74N", "EF75N", "EF78", "EF7N", "EF80N", "EF81N", "EF82N", "EF83N", "EF86N", "EF88N", "EF8N", "EF93N", "EF95N", "EF96N", "EF98N", "EF9N");
# my @InLines2 = ("ZI10", "ZI103", "ZI104", "ZI117", "ZI118N", "ZI126", "ZI134N", "ZI138", "ZI152", "ZI157", "ZI164", "ZI167", "ZI170", "ZI172", "ZI176", "ZI177", "ZI179", "ZI181", "ZI184", "ZI185", "ZI188", "ZI190", "ZI191", "ZI193", "ZI196", "ZI197N", "ZI198", "ZI199", "ZI200", "ZI202", "ZI206", "ZI207", "ZI211", "ZI212", "ZI214", "ZI216N", "ZI218", "ZI220", "ZI221", "ZI225", "ZI226", "ZI227", "ZI228", "ZI230", "ZI232", "ZI233", "ZI235", "ZI237", "ZI240", "ZI241", "ZI250", "ZI251N", "ZI252", "ZI253", "ZI254N", "ZI255", "ZI261", "ZI263", "ZI264", "ZI265", "ZI267", "ZI268", "ZI269", "ZI27", "ZI271", "ZI273N", "ZI276", "ZI279", "ZI28", "ZI284", "ZI286", "ZI291", "ZI295", "ZI296", "ZI311N", "ZI313", "ZI316", "ZI317", "ZI319", "ZI31N", "ZI320", "ZI321", "ZI324", "ZI329", "ZI33", "ZI332", "ZI333", "ZI335", "ZI336", "ZI342", "ZI344", "ZI348", "ZI351", "ZI352", "ZI353", "ZI359", "ZI362", "ZI368", "ZI370", "ZI373", "ZI374", "ZI377", "ZI378", "ZI379", "ZI380", "ZI381", "ZI382", "ZI386", "ZI392", "ZI394N", "ZI395", "ZI396", "ZI397N", "ZI398", "ZI400", "ZI402", "ZI405", "ZI413", "ZI418N", "ZI420", "ZI421", "ZI429", "ZI431", "ZI433", "ZI437", "ZI443", "ZI444", "ZI445", "ZI446", "ZI447", "ZI453", "ZI455N", "ZI456", "ZI457", "ZI467", "ZI468", "ZI472", "ZI476", "ZI477", "ZI488", "ZI490", "ZI491", "ZI504", "ZI505", "ZI508", "ZI50N", "ZI514N", "ZI517", "ZI523", "ZI524", "ZI527", "ZI56", "ZI59", "ZI61", "ZI68", "ZI76", "ZI81", "ZI85", "ZI86", "ZI90", "ZI91");

# Chr3R
# my @InLines1 = ('EF101N', 'EF102N', 'EF103N', 'EF108N', 'EF10N', 'EF110N', 'EF112N', 'EF114N', 'EF115N', 'EF116N', 'EF117N', 'EF119N', 'EF120N', 'EF122N', 'EF126N', 'EF12N', 'EF130N', 'EF131N', 'EF132N', 'EF135N', 'EF136N', 'EF15N', 'EF16N', 'EF19N', 'EF1N', 'EF2', 'EF24N', 'EF25N', 'EF26N', 'EF27N', 'EF31N', 'EF32N', 'EF35N', 'EF39', 'EF3N', 'EF43N', 'EF46N', 'EF48N', 'EF52N', 'EF53N', 'EF54N', 'EF59N', 'EF64N', 'EF65', 'EF66N', 'EF67N', 'EF69N', 'EF6N', 'EF70N', 'EF73N', 'EF74N', 'EF75N', 'EF78', 'EF7N', 'EF80N', 'EF81N', 'EF82N', 'EF83N', 'EF86N', 'EF88N', 'EF8N', 'EF93N', 'EF95N', 'EF96N', 'EF98N', 'EF9N');
# my @InLines2 = ('ZI10', 'ZI104', 'ZI114N', 'ZI117', 'ZI118N', 'ZI126', 'ZI129', 'ZI134N', 'ZI136', 'ZI138', 'ZI152', 'ZI157', 'ZI161', 'ZI164', 'ZI165', 'ZI167', 'ZI170', 'ZI172', 'ZI173', 'ZI176', 'ZI178', 'ZI179', 'ZI181', 'ZI182', 'ZI184', 'ZI185', 'ZI191', 'ZI193', 'ZI194', 'ZI196', 'ZI198', 'ZI200', 'ZI206', 'ZI207', 'ZI210', 'ZI211', 'ZI212', 'ZI213', 'ZI216N', 'ZI219', 'ZI225', 'ZI226', 'ZI231', 'ZI232', 'ZI235', 'ZI239', 'ZI240', 'ZI250', 'ZI251N', 'ZI252', 'ZI253', 'ZI254N', 'ZI255', 'ZI26', 'ZI261', 'ZI263', 'ZI265', 'ZI268', 'ZI269', 'ZI27', 'ZI271', 'ZI276', 'ZI279', 'ZI28', 'ZI281', 'ZI286', 'ZI295', 'ZI303', 'ZI311N', 'ZI313', 'ZI314', 'ZI316', 'ZI31N', 'ZI320', 'ZI321', 'ZI324', 'ZI329', 'ZI33', 'ZI332', 'ZI335', 'ZI336', 'ZI339', 'ZI342', 'ZI344', 'ZI348', 'ZI357N', 'ZI358', 'ZI359', 'ZI364', 'ZI365', 'ZI370', 'ZI373', 'ZI374', 'ZI377', 'ZI378', 'ZI386', 'ZI388', 'ZI392', 'ZI394N', 'ZI395', 'ZI398', 'ZI402', 'ZI413', 'ZI418N', 'ZI420', 'ZI421', 'ZI431', 'ZI433', 'ZI444', 'ZI445', 'ZI447', 'ZI453', 'ZI455N', 'ZI456', 'ZI457', 'ZI458', 'ZI460', 'ZI466', 'ZI467', 'ZI468', 'ZI471', 'ZI472', 'ZI477', 'ZI488', 'ZI490', 'ZI491', 'ZI504', 'ZI50N', 'ZI514N', 'ZI517', 'ZI524', 'ZI527', 'ZI530', 'ZI56', 'ZI61', 'ZI68', 'ZI76', 'ZI85', 'ZI86', 'ZI91', 'ZI99');


push @InLines, @InLines1;
push @InLines, @InLines2;
my @PopStarts = ();
my @PopStops = ();
my $i = 0;
push @PopStarts, 0;
$i += @InLines1;
push @PopStarts, $i;
$i = @InLines1 - 1;
push @PopStops, $i;
$i += @InLines2;
push @PopStops, $i;

if (@InLines != ($PopStops[-1]+1)){
  die "problem with pop starts/stops\t";
}

die if (@PopStarts != 2);
die if (@pops != 2);

my @InHandles = @InLines;
for ($i = 0; $i < @InHandles; $i++){
  $InHandles[$i] =~ s/1/A/;
  $InHandles[$i] =~ s/2/B/;
  $InHandles[$i] =~ s/3/C/;
  $InHandles[$i] =~ s/4/D/;
  $InHandles[$i] =~ s/5/E/;
  $InHandles[$i] =~ s/6/F/;
  $InHandles[$i] =~ s/7/G/;
  $InHandles[$i] =~ s/8/H/;
  $InHandles[$i] =~ s/9/I/;
  $InHandles[$i] =~ s/0/J/;
}

my $p = 0;
my $c = 0;
my $f = 0;
my $w = 0;
my $j = 0;
my $s = 0;
my $file = '';
my $file1 = '';
my $sum = 0;
my $RelDataProp = 0;
my $diffs = 0;
my $comps = 0;
my $SiteComps = 0;
my $SiteCov = 0;
my $kb = -1;
my $window = 0;
my $SplicePos = 0;
my $SpliceNum = 0;
my $pi = 0;
my $threshold = 0;
my $NSites = 0;
my $IndA = 0;
my $IndB = 0;
my $site = 0;
my $pos = 0;
my $length = 0;
my $LengthSum = 0;
my $ratio = 0;
my $MinComps = 0;
my $denominator = 0;
my $median = 0;

my @line = 0;
my @ChrWinOffset = ();
my @PiAoA = ();
my @PopPi = ();
my @RGSiteCov = ();
my @InputAoA = ();
my @NextAoA = ();
my @WinStarts = ();
my @WinStops = ();
my @LengthSums = ();
my @ratios = ();
my @TotalComps = ();
my @LengthSumsA = ();
my @SiteCompsA = ();
my @LengthSumsB = ();
my @SiteCompsB = ();
my @NSites = ();
my @denominators = ();
my @SortedValues = ();
my @RatioQuantiles = ();

#Get window boundaries
push @ChrWinOffset, 0;
for ($c = 0; $c < @chrs; $c++){
  $file = $WinFilePrefix . $chrs[$c] . ".txt";
  open W, "<$file" or die "Can't open window file $file\n";
  while (<W>){
    chomp;
    last if m/^$/;
    @line = split;
    push @WinStarts, $line[0];
    push @WinStops, $line[1];
  }
  close W;
  if ($c < (@chrs - 1) ){
    $sum = @WinStarts;
    push @ChrWinOffset, $sum;
  }
}

#go through each chromosome...
for ($c = 0; $c < @chrs; $c++){
  $kb = -1;
  $window = $ChrWinOffset[$c];
  
#Open all the data files
  @line = ();
  for ($f = 0; $f < @InLines; $f++){
    $file = $InLines[$f] . "_" . $chrs[$c] . "\.fas1k";
    open $InHandles[$f], "<../$chrs[$c]/$file" or print "Can't open ../$chrs[$c]/$file\n";
    push @{$InputAoA[0]}, @line;
  }
  $file1 = $InHandles[0];
  while (<$file1>){
#    @InputAoA = ();
    $kb++;
    chomp;
    last if m/^$/;
    @line = split(//, $_);
#    push @InputAoA, [ @line ];
    push @{$InputAoA[0]}, @line;
    for ($f = 1; $f < @InLines; $f++){
      $file = $InHandles[$f];
      $_ = (<$file>);
      chomp;
      @line = split(//, $_);
#      push @InputAoA, [ @line ];
      push @{$InputAoA[$f]}, @line;
    }

#Check whether we've got our full data - if so, transfer the relevant part to LastAoA
    if ( ($WinStops[$window]) <= (($kb*1000)+999) ){
      $SpliceNum = 999 - ($WinStops[$window] % 1000);
      $SplicePos = @{$InputAoA[0]} - $SpliceNum;
      @NextAoA = ();
      for ($i = 0; $i < @InputAoA; $i++){
	@line = splice( @{$InputAoA[$i]}, $SplicePos, $SpliceNum );
	push @NextAoA, [ @line ];
      }
      $window++;
      print "got data for $chrs[$c], window $window\n";

#Decide how long an identical tract must be to be counted
      $NSites = @{$InputAoA[0]};
      $threshold = $NSites * $WinProp;
#      $MinComps = $NSites * 6;

#Ignore (recode) singleton polymorphisms
      for ($s = 0; $s < $NSites; $s++){
	$As = 0;
	$Cs = 0;
	$Gs = 0;
	$Ts = 0;
	$Ns = 0;
	for ($i = 0; $i < @InputAoA; $i++){
	  if ($InputAoA[$i][$s] eq 'A'){
	    $As++;
	  }
	  elsif ($InputAoA[$i][$s] eq 'C'){
	    $Cs++;
	  }
	  elsif ($InputAoA[$i][$s] eq 'G'){
	    $Gs++;
	  }
	  elsif ($InputAoA[$i][$s] eq 'T'){
	    $Ts++;
	  }
	  else{
	    $Ns++;
	  }
	}
	$MajorAllele = 'A';
	$MajorAlleleCount = $As;
	if ($Cs > $MajorAlleleCount){
	  $MajorAllele = 'C';
	  $MajorAlleleCount = $Cs;
	}
	if ($Gs > $MajorAlleleCount){
	  $MajorAllele = 'G';
	  $MajorAlleleCount = $Gs;
	}
	if ($Ts > $MajorAlleleCount){
	  $MajorAllele = 'T';
	  $MajorAlleleCount = $Ts;
	}
	last if ($MajorAlleleCount < 1.5);
	if ($As == 1){
	  for ($j = 0; $j < @InputAoA; $j++){
	    if ($InputAoA[$j][$s] eq 'A'){
	      $InputAoA[$j][$s] = $MajorAllele;
	      last;
	    }
	  }
	}
	if ($Cs == 1){
	  for ($j = 0; $j < @InputAoA; $j++){
	    if ($InputAoA[$j][$s] eq 'C'){
	      $InputAoA[$j][$s] = $MajorAllele;
	      last;
	    }
	  }
	}
	if ($Gs == 1){
	  for ($j = 0; $j < @InputAoA; $j++){
	    if ($InputAoA[$j][$s] eq 'G'){
	      $InputAoA[$j][$s] = $MajorAllele;
	      last;
	    }
	  }
	}
	if ($Ts == 1){
	  for ($j = 0; $j < @InputAoA; $j++){
	    if ($InputAoA[$j][$s] eq 'T'){
	      $InputAoA[$j][$s] = $MajorAllele;
	      last;
	    }
	  }
	}
      }

#Within each population, cycle through each pair of sequences
      @LengthSums = ();
      @TotalComps = ();
      for ($p = 0; $p < @pops; $p++){
	$SiteComps = 0;
	$LengthSum = 0;
	for ($IndA = $PopStarts[$p]; $IndA < $PopStops[$p]; $IndA++){
	  for ($IndB = $IndA + 1; $IndB <= $PopStops[$p]; $IndB++){

#Record sites that differ between this pair of sequences
	    $IdSites = 0;
	    for ($site = 0; $site < $NSites; $site++){
	      next if (($InputAoA[$IndA][$site] eq 'N') || ($InputAoA[$IndB][$site] eq 'N'));
	      $SiteComps++;
	      if ($InputAoA[$IndA][$site] eq $InputAoA[$IndB][$site]){
		$IdSites++;
	      }
	      else{
		if ($IdSites > $threshold){
		  $LengthSum += $IdSites;
		}
		$IdSites = 0;
	      }
	    }
	    if ($IdSites > $threshold){
	      $LengthSum += $IdSites;
	    }
	  }
	}
	push @LengthSums, $LengthSum;
	push @TotalComps, $SiteComps;
	if ($p == 0){
	  push @LengthSumsA, $LengthSum;
	  push @SiteCompsA, $SiteComps;
	  push @NSites, $NSites;
	}
	else{
	  push @LengthSumsB, $LengthSum;
	  push @SiteCompsB, $SiteComps;
	}
      }
#      if ($LengthSums[1] == 0){
#	$LengthSums[1] = $threshold / 2;
#      }
#      if ($LengthSums[0] == 0){
#	$LengthSums[0] = $threshold / 2;
#      }
#      if (($TotalComps[0] < $MinComps) || ($TotalComps[1] < $MinComps)){
#	$ratio = -999;
#	push @ratios, $ratio;
#      }
#      else{
#	$ratio = log (($LengthSums[0] / $TotalComps[0]) /  ($LengthSums[1] / $TotalComps[1])) ;
#	push @ratios, $ratio;
#      }
#      print " - ratio = $ratio\n";
      @InputAoA = @NextAoA;
    }
  }
}

#calculate ratio, using max of local and global denominators
#to avoid the lowest pop 2 haplotype-sharing values driving outliers
for ($i = 0; $i < @LengthSumsB; $i++){
  $MinComps = $NSites[$i] * 6;
  if (($SiteCompsA[$i] > $MinComps) && ($SiteCompsB[$i] > $MinComps)){
    $denominator = $LengthSumsB[$i] / $SiteCompsB[$i];
    push @denominators, $denominator;
  }
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
for ($i = 0; $i < @LengthSumsA; $i++){
  $MinComps = $NSites[$i] * 6;
  if (($SiteCompsA[$i] > $MinComps) && ($SiteCompsB[$i] > $MinComps)){
    $denominator = $LengthSumsB[$i] / $SiteCompsB[$i];
    if ($denominator > $median){
      $ratio = ($LengthSumsA[$i] / $SiteCompsA[$i]) / $denominator;
      push @ratios, $ratio;
    }
    else{
      $ratio = ($LengthSumsA[$i] / $SiteCompsA[$i]) / $median;
      push @ratios, $ratio;
    }
  }
  else{
    $ratio = -999;
    push @ratios, $ratio;
  }
}

@SortedValues = @ratios;
@SortedValues = sort { $b <=> $a } @SortedValues;
for ($i = (@SortedValues - 1); $i >= 0; $i--){
  if ($SortedValues[$i] < -0.9){
    splice @SortedValues, $i, 1;
  }
  else{
    last;
  }
}
for ($i = 0; $i < @ratios; $i++){
  if ($ratios[$i] < -0.9){
    push @RatioQuantiles, -1;
    next;
  }
  for ($j = 0; $j < @ratios; $j++){
    if ($ratios[$i] >= $SortedValues[$j]){
      $quantile = ($j / @SortedValues);
      push @RatioQuantiles, $quantile;
      last;
    }
  }
}



#Output ratio for each window/pop (transposing PiAoA)
open O, ">$OutputFile";
$c = 0;
print O "chr\tWinStart\tWinStop\tNSites\tPop1SiteComps\tPop1LengthSum\tPop2SiteComps\tPop2LengthSum\tPopLogRatio\tRatioQuantile\n";
for ($w = 0; $w < @WinStarts; $w++){
#  $window = $w - $chrs[$c];
  print O "$chrs[$c]\t$WinStarts[$w]\t$WinStops[$w]";
#  if ( ($c < (@chrs-1)) && ( ($w+1) == $chrs[$c+1]) ){
  if ( (defined($WinStarts[$w+1])) && ($WinStarts[$w+1] < $WinStops[$w]) ){
    $c++;
  }
  print O "\t$NSites[$w]\t$SiteCompsA[$w]\t$LengthSumsA[$w]\t$SiteCompsB[$w]\t$LengthSumsB[$w]\t$ratios[$w]\t$RatioQuantiles[$w]\n";
}
