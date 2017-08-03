#!/usr/bin/perl
#
###@StandardDocs 
#  This perl script takes input clc run directory with subdir contains 
#  bam file(s), a panel specific configure file and call dlmp_cvgFbamAnnotCont.pl per 
#  sample amd the summarize data and output coverage_summary.xls file
#  
#  INPUT: 
#
#  OUTPUT
#    file with multiple sessions with coverage per sample, per lane, pergene, and per gene and exon 
#
# HISOTY
#  -- 20150926 Zheng Jin Tu:   add dup and mq0
#  -- 20151014 Zheng Jin Tu:   add contamination check portion
#  -- 20151018 Zheng Jin Tu:   add picard duplication rate portion
#  -- 20170127 Zheng Jin Tu:   add code to deal overlap within bed and intro/exon repeat names
#  -- 20170214 Zheng Jin Tu:   add code to read merged.cvg1 file due to cava style target.bed with closeby frag
#  -- 20170316 Zheng Jin Tu:   correct for loop with i to t change  else { for ( $t = 0; $t < $cntchar; $t++ ){ $outstr .= "\t0"; }}, old way will hang in for loop if any folder without any fastq count
#  -- 20170330 Zheng Jin Tu:   use cava/<sample>_cmb.vcf file if it exists
#
my $Prg = "sum_cvgFbam.pl";

sub usage { print "\nUsage: $Prg ";
            print "\n\t\t -id indir          # REQUIRED input top directory";
            print "\n\t\t  -c configfile     # REQUIRED panel specific configure file";
            print "\n\t\t -od outdir         # REQIURED outdir where output write to";
            print "\n\t\t -ot output         # REQUIRED output name suggest runbatch name as output name";
            print "\n\t\t -td tmpdir         # REQUIRED outdir name where multiple cvg files located";
            print "\n\t\t[-ss samplesheet    # OPTION Illumina SampleSheet.csv or Ion Torrent SampleSheet.txt file";
            print "\n\t\t  -t technology     # CHOICE: illumina, iontorrent; DEFAULT illumina";
            print "\n\t\t  -l label          # Torrent sample name (I006-12345678) has label CHOICE: yes no, DEFAULT yes";
            print "\n\t\t -dg debug          # for debug use DEFAULT no]\n\n"; }

use File::Basename;
use Getopt::Long;
use Data::Dumper;

my ( $indir, $configfile, $outdir, $output, $tmpdir, $samplesheet, $technology, $label, $debug );

if ($#ARGV < 3) { warn "number of commandline args: $#ARGV\n" if $debug; &usage; exit 1; }

my @Options =
(
    "id=s"  => \$indir,
    "c=s"   => \$configfile,
    "od=s"  => \$outdir,
    "td=s"  => \$tmpdir,
    "ot=s"  => \$output,
    "ss=s"  => \$samplesheet,
    "t=s"   => \$technology,
    "l=s"   => \$label,
    "dg=s"  => \$debug
);
if (!&GetOptions(@Options) ){ &usage; exit 1;                                                     }
if (!$indir       ){ warn "\nMust specify an input directory   \n"; &usage; exit 1;               }
if (!$tmpdir      ){ warn "\nMust specify a tmpdir for cvg file\n"; &usage; exit 1;               }
if (!$configfile  ){ warn "\nMust specify panel configfile file\n"; &usage; exit 1;               }
if (!$outdir      ){ warn "\nMust specify an outdir            \n"; &usage; exit 1;               }
if (!$output      ){ $output      = "$Prg"; $output =~ s/\.pl//;                                  }
if (!$technology  ){ $technology  = "illumina";                                                   }
if (!$label       ){ $label       = "yes";                                                        }
if (!$debug       ){ $debug       = 0;                                                            }
  
# --------------------------------------------------------------------------------
# first read in configure file
# --------------------------------------------------------------------------------
my ( $cvgTargetBed, $cvgGenomeBed, $zerobasebed, $genelist, $desiredCoverageIllumina, $desiredCoverageTorrent, $variantmincvg, $allelefrequency );
open ( IN, "< $configfile" ) or print "\nFail to open $configfile\n\n";
while ( <IN> ) {
  if ( !/^\#/ ){
    chomp; my @tp = split ( /\"/, $_ );
    if    ( /^cvgTargetBed=/              ){ $cvgTargetBed               = $tp[1]; }
    elsif ( /^cvgGenomeBed=/              ){ $cvgGenomeBed               = $tp[1]; }
    elsif ( /^genelist=/                  ){ $genelist                   = $tp[1]; }
    elsif ( /^desiredCoverageIllumina=/   ){ $desiredCoverageIllumina    = $tp[1]; }
    elsif ( /^desiredCoverageTorrent=/    ){ $desiredCoverageTorrent     = $tp[1]; }
    elsif ( /^variantmincvg=/             ){ $variantmincvg              = $tp[1]; }
    elsif ( /^allelefrequency=/           ){ $allelefrequency            = $tp[1]; }
  }
}
close ( IN );

my $desiredCoverage;
if    ( $desiredCoverageIllumina && $technology eq "illumina"   ){ $desiredCoverage = $desiredCoverageIllumina; }
elsif ( $desiredCoverageTorrent  && $technology eq "iontorrent" ){ $desiredCoverage = $desiredCoverageTorrent;  }

if ( $debug ){ print "desiredCoverage=$desiredCoverage $technology desiredCoverageIllumina=$desiredCoverageIllumina $desiredCoverageTorrent\n"; }

if ( !$cvgGenomeBed    ){ print "\nERROR: can not identify cvgGenomeBed in $configfile\n"; exit 1; }
if ( !$cvgTargetBed    ){ print "\nERROR: can not identify cvgTargetBed in $configfile\n"; exit 1; }
if ( !$variantmincvg   ){ $variantmincvg   = 10;  }
if ( !$allelefrequency ){ $allelefrequency = 5;   }

if ( !$cvgTargetBed ) { print "\nERROR: can not identify cvgTargetBed or in $configfile\n"; exit 1; }

if ( opendir DH, $outdir ){ closedir DH; } else { `mkdir $outdir`; }

# -------------------------------------------------------------------------------------------
# read $cvgTargetBed and mark region to be checked 
# -------------------------------------------------------------------------------------------
$cntline = 0; $cntpos = 0;
open ( IN, "< $cvgTargetBed" ) or print "\nFail to open $cvgTargetBed\n" and die;
while ( <IN> ) {
  # chr10   88681246        88681482        BMPR1A:NM_004329        9       +
  chomp; @inline = split ( /\t/, $_ );   $chr = $inline[0];   
  @inline2 = split ( /\:/, $inline[3] ); $gene = $inline2[0];
  if ( $inline[1] < $inline[2] ) { $begin = $inline[1]; $end = $inline[2]; } else { $begin = $inline[2]; $end = $inline[1]; }

  if ( $zerobasebed eq "yes" ){ $begin++; }

  # check any gene is too big (>50MB) that is typo error for a few cases in the past
  if ( abs ( $end - $begin ) > 50000000 ) { print "\n$end - $begin \> 50 Mb in line\n$_\n"; exit 1; }

  for ( $i = $begin; $i <= $end; $i++ ){ 
    if ( !exists $tgtdna_hash{$chr}{$i} ){ $tgtdna_hash{$chr}{$i} = $gene; $cntpos++; }
  }
  if ( $cntline % 100 == 1 ) { print "  $cntline\t$chr\t$begin\t$end\t$gene\t$cntpos\n"; }
  $cntline++;
}
close ( IN );
print "$cvgTargetBed has record of $cntline line(s) with $cntpos nucleotides\n";

# -------------------------------------------------------------------------------------------
# if exists, read $bedseq and record ref_hash
# -------------------------------------------------------------------------------------------
$cntlocseq = 0;
if ( open ( IN, "< $bedseq1" )){ 
  while ( <IN> ) {
    chomp; @inline = split( /\t/, $_ ); 
    $ref_hash{$inline[0]}{$inline[1]} = $inline[2]; 
    $cntlocseq++;
  } 
  close( IN );
  print "\n$bedseq1 has $cntlocseq lines ...\n";
}

# -------------------------------------------------------------------------------------------
# if exists, read $dbsnp file
# -------------------------------------------------------------------------------------------
$cntdbsnp = 0; $cntline = 0;
if ( open ( IN, "< $dbsnp"  )) {
  print "\nOpen $dbsnp for reading ...\n";
  while ( <IN> ) {
    # 585  chr1   10144  10145   rs144773400  0   +  A   A  -/A  genomic deletion  unknown 0  0  unknown exact   1  1  BL, 0
    @inline = split ( /\t/, $_ );
    if ( $inline[2] < $inline[3] ) { $start = $inline[2]; $stop = $inline[3]; } else { $start = $inline[3]; $stop = $inline[2]; }
    for ( $i = $start; $i <= $stop; $i++ ) {
      if ( exists $tgtdna_hash{$inline[1]}{$i} ) {
        $dbsnp_hash{$inline[1]}{$i} = $inline[4]; $cntdbsnp++;
        $dbsnpvar_hash{$inline[1]}{$i} = "$inline[8]\>\>$inline[9]";
      }
    }
    if ( $cntline % 10000 == 1 ) { print " $cntline\t$cntdbsnp\t$inline[1] $inline[2] $inline[4]\n"; }
    $cntline++;
  }
  close ( IN );
  print "$dbsnp has cntdbsnp=$cntdbsnp\n";
}

# -------------------------------------------------------------------------------------------
# first round read $indir and creat list of dir_lists
# -------------------------------------------------------------------------------------------
$cntdirlist = 0;
opendir DH, $indir or print "\n\tFail to open $indir\n" and die; @dir_list = grep !/^\./, readdir(DH);
foreach $dir_list ( @dir_list ){ 
  if ( $dir_list =~ /log/ ){} # to avoid log or logs folder hold process
  else { 
    $dirlist_hash{$dir_list} = $dir_list; 
    $dirlistArr[$cntdirlist] = $dir_list; 
    $cntdirlist++; 
  }
}
closedir DH;

# -------------------------------------------------------------------------------------------
# if exists, read $samplesheet and mark samples to be checked
# -------------------------------------------------------------------------------------------
if ( open ( IN, "< $samplesheet" ) ) {
  $cntsample = 0;
  while ( <IN> ) {
    if ( $technology eq "iontorrent" ) {
      # IonXpress_001    1095240         40,276,721      36,251,699      343,930         117 bp         BAM BAI
      @inline = split ( /\t/, $_ );
      $IXnum = $inline[0]; $IXid = $inline[1]; $IXid =~ s/\s//g;

      if ( $label eq "yes" || $label eq "Yes" ){ if ( $IXnum =~ /IonXpress_/ ){ $IXnum =~ s/IonXpress_//; } $sampleID = "I$IXnum-$IXid"; }
      else { $sampleID = $IXid; }

      if ( $inline[0] eq "No barcode" ) { $sampleID = "noBarcode"; }
      $sample_hash{$sampleID} = $sampleID;
      $sampleArr[$cntsample] = $sampleID;
      $cntsample++;
    }
    elsif ( $technology eq "illumina" ) {
      # 000000000-AA6VC,1,1095240,human,ATCACG,unknown,N,Solid_Tumor_Targeted_Cancer_Panel-valid:2013-08-09-MiSeq,CGSL,CAPN_35_AA6VC_MS
      if ( !/FCID/ ) {
        @inline = split ( /\,/, $_ ); $sampleID = $inline[2];
        $sample_hash{$sampleID} = $sampleID;
        $sampleArr[$cntsample] = $sampleID;
        $cntsample++;
      }
    }
    print "cntsample=$cntsample\tsampleID=$sampleID\n";
  }
  close ( IN );

  if ( $technology eq "illumina" ) {
    for ( $l = 0; $l <= 8; $l++ ) {
      $sampleID = "lane$l";
      if ( exists $dirlist_hash{$sampleID} && !exists $sample_hash{$sampleID} ) { 
       $sample_hash{$sampleID} = $sampleID; $sampleArr[$cntsample] = $sampleID; $cntsample++; 
      }
    }
  }
  elsif ( $technology eq "iontorrent" ) { 
    $sampleID = "noBarcode"; 
    if ( exists $dirlist_hash{$sampleID} && !exists $sample_hash{$sampleID} ) {
      $sample_hash{$sampleID} = $sampleID; $sampleArr[$cntsample] = $sampleID; $cntsample++;
    }
  }
  print "$samplesheet has $cntsample samples\n";
}
else { $cntsample = $cntdirlist; @sampleArr = @dirlistArr; } # %dirlist_hash = (); }

print "\nsampledirlist\t";
$c = 0; foreach $key ( %dirlist_hash ) { if ( $dirlist_hash{$key} ne "" ) { $c++; print "$c $key $dirlist_hash{$key} "; }}
print "\n\nsamplefilelist\t";
for ( $s = 0; $s < $cntsample; $s++ ) { $sp = $s + 1; print " $sp $sampleArr[$s] "; }
print "\n\n";

# -------------------------------------------------------------------------------------------
# read in $indir and record info
# -------------------------------------------------------------------------------------------
$cntbamfile = 0; $cntfqfile = 0;
for ( $s = 0; $s < $cntsample; $s++ ) {
  $dir_list = $sampleArr[$s];  # print "$corename ";
  $sindir = "$indir/$dir_list";  $cntcmbvcf = 0;
  if ( opendir SDH, $sindir ) { 
    $tgtvcfoutstr = ""; $ctmvcfoutstr = "";  $cntrawsnp = 0; $cntrawdip = 0; $cnttgtsnp = 0; $cnttgtdip = 0; $cntctmsnp = 0; $cntctmdip = 0;

    @sdir_list = grep !/^\./, readdir(SDH); 

    if ( $dir_list =~ /clc/ || $dir_list =~ /CLC/ || $dir_list =~ /Clc/ ){} 
    else { 
      foreach $sdir_list ( @sdir_list ) { 
        if ( $sdir_list =~ /extended/ ) {} 
        elsif ( $sdir_list =~ /fastq.gz$/ ) { 
          $lane = ""; $corename = "";
          @tpl = split ( /\./, $sdir_list ); # 800425_CAACT_L001_R2_001 mapping SNP Detection Table.csv
          for ( $t1 = 0; $t1 < scalar(@tpl); $t1++ ){
            @tp2 = split ( /_/, $tpl[$t1] ); if ( $t1 == 0 ) { $corename = $tp2[0]; }
            for ( $t2 = 0; $t2 < scalar(@tp2); $t2++ ) {
              $cf = $tp2[$t2];
              if ( $cf eq "L001" || $cf eq "L002" || $cf eq "L003" || $cf eq "L004" || $cf eq "L005" || $cf eq "L006" || $cf eq "L007" || $cf eq "L008" ) {
                $lane = $cf; 
              }
            }
          }
          if ( $lane ne "" && $corename ne "" ){ $lane_hash{$corename} = $lane; $cntlane_hash{$corename}++; }

          $cntfqfile++; if ( $debug ){ print "  $cntfqfile\t$exec"; }
        }
        elsif ( $sdir_list =~ /bam$/ ){
          $cntbamfile++; if ( $debug ){ print "$cntbamfile\t$bedtools $indir/$dir_list/$sdir_list ...\n"; }
        }

        %variant_hash = ();
        if ( $sdir_list =~ /_cmb.vcf/ ){ # for CLC v6.0
          print "\tbam=$cntbamfile\tfastq=$cntfqfile\t$corename\t$lane\t$sdir_list\n";
          # check cava/<sample>_cmb.vvf first
          $cmbvcfile = "$indir/$dir_list/cava/$sdir_list";
          if ( -e $cmbvcfile ){} else { $cmbvcfile = "$indir/$dir_list/$sdir_list"; }

          if ( open ( IN, "< $cmbvcfile" ) ) { 
            $cntline = 0; $cntrawTS = 0; $cnttgtTS = 0;
            $cntcmbvcf++; if ( $cntcmbvcf > 1 ){ print "\nERROR  there is more than one cmb_vcf files $dir_list\n"; exit 1; } 
            while ( <IN> ) {
              if ( /^\#/ ) { $vcfoutstr .= $_; }
              else {
                chomp;
                @inline  = split ( /\t/, $_ );           # chr1  120612040  .  T  TCCTCCGCCG  .  .  .  GT:AD:DP  0/1:1690,598:2382
                @inline1 = split ( /\:/, $inline[9] );   # 0/1:1690,598:2382
                @inline2 = split ( /\//, $inline1[0] );  # 0/1
                @inline3 = split ( /\,/, $inline1[1] );  # 1690,598
                $chr = $inline[0]; $refpos = $inline[1]; $ref = $inline[3]; $reflen = length($ref); $var = $inline[4]; $varlen = length($var);
                $variant = "$chr:g.$refpos$ref\>$var";
                $depth = $inline1[2]; $vardepth = $inline3[1]; 
   
                # check variant whether meet the filters
                $passflag = 0;
                if ( $depth > $variantmincvg ){                 # must greater than $variantmincvg
                  if ( $inline2[0] == 1 ) { $passflag = 1; }    # 1/1 a change of nucleotide
                  else { 
                    $allelefrq = 100 * $vardepth / $depth;
                    if ( $inline1[1] >= $variantmincvg && $allelefrq >= $allelefrequency ){
                      $passflag = 1;
                    }
                  } 
                }
                
                if ( exists $variant_hash{$variant} ){} # use hash to remove any repeat/overlap 
                else { 
                  $variant_hash{$variant} = 1;
                  if ( $reflen == 1 && $varlen == 1 ){ # SNP case
                    $cntrawsnp++;

                    # for transition and transversion
                    if    (( $ref eq "A" || $ref eq "a" ) && ( $var eq "G" || $var eq "g" )){ $cntrawTS++; }
                    elsif (( $ref eq "G" || $ref eq "g" ) && ( $var eq "A" || $var eq "a" )){ $cntrawTS++; }
                    elsif (( $ref eq "C" || $ref eq "c" ) && ( $var eq "T" || $var eq "t" )){ $cntrawTS++; }
                    elsif (( $ref eq "T" || $ref eq "t" ) && ( $var eq "C" || $var eq "c" )){ $cntrawTS++; }

                    if ( $passflag && exists $tgtdna_hash{$chr}{$refpos} ){ 
                      $tgtsnp_hash{$dir_list} .= $_; $tgtvcfoutstr .= "$_\n"; $cnttgtsnp++; 
                      if    (( $ref eq "A" || $ref eq "a" ) && ( $var eq "G" || $var eq "g" )){ $cnttgtTS++; }
                      elsif (( $ref eq "G" || $ref eq "g" ) && ( $var eq "A" || $var eq "a" )){ $cnttgtTS++; }
                      elsif (( $ref eq "C" || $ref eq "c" ) && ( $var eq "T" || $var eq "t" )){ $cnttgtTS++; }
                      elsif (( $ref eq "T" || $ref eq "t" ) && ( $var eq "C" || $var eq "c" )){ $cnttgtTS++; }
                    }
                  }
                  else { # indel and polymorphic case 
                    $cntrawdip++;
                    if ( $passflag ){ 
                      if ( exists $tgtdna_hash{$chr}{$refpos} ){ $tgtdip_hash{$dir_list} .= $_; $tgtvcfoutstr .= "$_\n"; $cnttgtdip++; }
                    } 
                  }
                }
                $cntline++;
              }
            } 
            close ( IN ); 
            $cntrawsnp_hash{$dir_list} = $cntrawsnp;  $cnttgtsnp_hash{$dir_list} = $cnttgtsnp;
            $cntrawdip_hash{$dir_list} = $cntrawdip;  $cnttgtdip_hash{$dir_list} = $cnttgtdip;
            $cntrawTS_hash{$dir_list}  = $cntrawTS;   $cnttgtTS_hash{$dir_list}  = $cnttgtTS;

            # print "cntrawTS=$cntrawTS\n";
            if ( $tgtvcfoutstr ne "" ) { 
              open ( OT, "> $outdir/$dir_list.vcf" ) or print "\nFail to open $outdir/$dir_list.vcf\n\n" and die; print OT $tgtvcfoutstr; close OT;

              # create annotated vcf file for better reading
              $annotvcfstr = "";
              @vcfline = split ( /\n/, $tgtvcfoutstr );
              for ( $v = 0; $v < scalar(@vcfline); $v++ ) { 
                if ( $vcfline[$v] =~ /^\#/ ) { $annotvcfstr .= "$vcfline[$v]\n"; }
                else { 
                  @vcfline2 = split ( /\t/, $vcfline[$v] );
                  $vcfchr = $vcfline2[0]; $vcfloc = $vcfline2[1]; $vcfp = $vcfloc + 1;

                  $nstr = ""; 
                  for ( $n = $vcfloc-2; $n <= $vcfloc+2; $n++ ){
                    if ( exists $ref_hash{$vcfchr}{$n} ){ $nstr .= $ref_hash{$vcfchr}{$n}; } else { $nstr .= "-"; }
                  }
                  if ( exists $dbsnp_hash{$vcfchr}{$vcfp} ){ $snpstr  = $dbsnp_hash{$vcfchr}{$vcfp}; } else { $snpstr  = "-"; }
                  if ( exists $tgtdna_hash{$vcfchr}{$vcfloc} ){ $vcfgene = $tgtdna_hash{$vcfchr}{$vcfloc}; } else { $vcfgene = "-"; } 
                  $annotvcfstr .= "$vcfline[$v]\t$vcfgene\t$snpstr\t$nstr\n";
                }
              }
              open ( OT, "> $outdir/$dir_list.annot.vcf" ) or print "\nFail to open $outdir/$dir_list.annot.vcf\n\n" and die; 
              print OT $annotvcfstr; close ( OT ); 
            }
          }
        } # ENDOF  if ( $sdir_list =~ /_cmb.vcf/ ){ # for CLC v6.0
      } # foreach $sdir_list ( @sdir_list ) {
    } # ENDOF  if ( $dir_list =~ /clc/ || $dir_list =~ /CLC/ || $dir_list =~ /Clc/ ){} else {
  } # ENDOF if ( opendir SDH, $sindir ) {
  closedir SDH;
} # ENDOF for ( $s = 0; $s < $cntsample; $s++ ) {
 
# -----------------------------------------------------------------
# open $tmpdir to summarize data
# -----------------------------------------------------------------
$flag = 0; $totfastqnum = 0; 
for ( $s = 0; $s < $cntsample; $s++ ){ 
  $corename = $sampleArr[$s];
  $filename = "$tmpdir/$corename.zcat";
  if ( open ( IN, "< $filename" )){
    $sum = 0; while (<IN>) { chomp; $sum += $_; } close( IN ); $fastqsum_hash{$corename} = $sum / 4; $totfastqnum += $sum / 4;
    if ( $sum ) { $flag = 1; }
  } 

  $filename = "$tmpdir/$corename.dup.mq0";
  if ( open ( IN, "< $filename" )){
    while (<IN>){ chomp; @tp = split ( /\t/, $_ ); $tgtdup_hash{$corename} = "$tp[0]\t$tp[1]\t$tp[2]\t$tp[3]"; $mq0_hash{$corename} = $tp[4]; } close( IN );
  }

  $filename = "$tmpdir/$corename.picard.matrics";
  if ( open ( IN, "< $filename" )){
    $cntline = 0; $dupflag = 0;
    while ( <IN> ){ 
      chomp; @tp = split ( /\t/, $_ ); 
      if ( $dupflag ){ $pctdup_hash{$corename} = sprintf ( "%.2f", $tp[7] ); } 
      if ( /PERCENT_DUPLICATION/ ){ $dupflag = 1; } else { $dupflag = 0; } 
      $cntline++;
    }
    close ( OT );
  }
}

$cntAvgLow = 0; # used to monitor the whole exon or amplicon drop out
for ( $s = 0; $s < $cntsample; $s++ ) {
  $corename = $sampleArr[$s]; print "$corename ";
  $filename = "$tmpdir/$corename.cvg1";

  $cntline = 0;
  if ( open ( IN, "< $filename" )){
    while (<IN>) { 
      # chr4    1803551 1803653 FGFR3   1       +       9288    102     102     1.0000000
      chomp; @inline = split( /\t/, $_ ); 
      $tmpcvg = $inline[scalar(@inline)-4]; 
      $exonum = "$inline[4] $inline[0]:$inline[1]-$inline[2]"; 
      @tp = split ( /\:/, $inline[3] );

      if ( $cntline == 0 ) { 
        $tmpsum = 0; $generec = $tp[0]; # APC:NM_000038.5
        if ( $tmpcvg ) { $exonCvg[$cntexon] = $tmpcvg; $exonNumArr[$cntexon] = $exonum; $tmpsum += $tmpcvg; $cntexon++; }
      }
      elsif ( $generec eq $tp[0] && $tmpcvg ) { 
        $exonCvg[$cntexon] = $tmpcvg; $tmpsum += $tmpcvg; $exonNumArr[$cntexon] = $exonum; $cntexon++;
      }
      else { 
        if ( $cntexon && $tmpsum ) { 
          $tmpavg = $tmpsum / $cntexon;  
          for ( $e = 0; $e < $cntexon; $e++ ) { 
            if ( $tmpavg && $exonCvg[$e] ) { 
              $tmpratio = $exonCvg[$e] / $tmpavg; 
              $exongene = "$exonNumArr[$e]\t$generec";
              if ( $debug ){ print "corename=$corename exongene=$exongene tmpratio=$tmpratio\n"; }
              if ( $tmpratio < 0.1 ) { $cvgexongene_hash{$corename}{$exongene} = sprintf( "%.2f\tAvgLow", $tmpratio); }
              else { $cvgexongene_hash{$corename}{$exongene} = sprintf( "%.2f\t", $tmpratio); }
            }
            else { print "Error1 $Prg abort.. filename=$filename\n\n"; exit 1; }
          }
        }
        @exonCvg = "";  $cntexon = 0; $tmpsum = 0;
        if ( $tmpcvg ) { $exonCvg[$cntexon] = $tmpcvg; $exonNumArr[$cntexon] = $exonum; $tmpsum += $tmpcvg; $cntexon++; }
        $generec = $tp[0]; 
      }
      $cntline++;
    } 
    # -----------------------------------------------------------------------------------
    # for last gene
    #
        if ( $cntexon && $tmpsum ) { # avoid single exon and at the end 
          $tmpavg = $tmpsum / $cntexon;
          for ( $e = 0; $e < $cntexon; $e++ ) {
            if ( $tmpavg && $exonCvg[$e] ) {
              $tmpratio = $exonCvg[$e] / $tmpavg;
              $exongene = "$exonNumArr[$e]\t$generec";
              if ( $tmpratio < 0.1 ) { $cvgexongene_hash{$corename}{$exongene} = sprintf( "%.2f\tAvgLow", $tmpratio); }
              else { $cvgexongene_hash{$corename}{$exongene} = sprintf( "%.2f\t", $tmpratio); }
            }
            else { print "Error1 $Prg abort..filename=$filename\n\n"; exit 1; }
          }
          $sum += $tmpsum;
        }
    #   
    # ENDOF for last gene
    # --------------------------------------------------------------------------------
    close ( IN );
  }

  $filename = "$tmpdir/$corename.merged.cvg1";
  $cntline = 0; $sum = 0;
  if ( open ( IN, "< $filename" )) {
    while (<IN>) {
      # chr4    1803551 1803653  9288    102     102     1.0000000
      chomp; @inline = split( /\t/, $_ );
      $tmpcvg = $inline[scalar(@inline)-4]; $sum += $tmpcvg;
    }
    close ( IN );
  }
  $cvg1sum_hash{$corename} = $sum;

  $filename = "$tmpdir/$corename.cvg2";
  if ( open ( IN, "< $filename" )) {
    $sum = 0; $phixsum = 0;
    while (<IN>){ 
      chomp; @inline = split( /\t/, $_ ); 
      if ( $inline[0] eq "phix" ){ $phixsum += $inline[6]; } 
      else { $sum += $inline[scalar(@inline)-4]; }
    } close( IN ); 
    $cvg2sum_hash{$corename} = $sum; $phixsum_hash{$corename} = $phixsum; $totsumphix += $phixsum;
    if ( $sum ){ $flag = 1; }
  }
}

$cnttotN = 0; $sumtotcvg = 0; $cnttotdeptmore = 0;
for ( $s = 0; $s < $cntsample; $s++ ) {
  $corename = $sampleArr[$s]; %chrloc_hash = ();
  $filename = "$tmpdir/$corename.cvg";
  $cntN = 0; $cntdeptmore = 0; $smplmax = 0; $smplmin = 1000000; $smplavgcvg = 0; $sumcvg = 0; $pctmore = 0; $cntrepeat = 0;
  if ( open ( IN, "< $filename" )) {
    $smpl = $corename;
    while ( <IN> ){
      # chr1    156106702       156106829       LMNA    8       +       3       2151
      chomp; @inline = split( /\t/, $_ ); $chr = $inline[0]; $loc = $inline[1] + $inline[scalar(@inline)-2];
      if ( $inline[0] ne "phix" && exists $tgtdna_hash{$chr}{$loc} ){
        if ( exists $chrloc_hash{$chr}{$loc} ){ $cntrepeat++; } # avoid repeat count if there is overlap region within bed file
        else { 
          $chrloc_hash{$chr}{$loc} = 1;      
          $tmpcvg = $inline[scalar(@inine)-1];
          $sumcvg += $tmpcvg; 
          $cntN++;
          if ( $tmpcvg >  $smplmax          ){ $smplmax = $tmpcvg; }
          if ( $tmpcvg <  $smplmin          ){ $smplmin = $tmpcvg; }
          if ( $tmpcvg >= $desiredCoverage ){ $cntdeptmore++;     }
        }
      }
    } 
    close(IN);
    print "Read in $filename and cntN=$cntN\tcntrepeat=$cntrepeat\n";

    if ( $cntN ) { $smplavgcvg = $sumcvg / $cntN; $pctmore = 100 * $cntdeptmore / $cntN; $pctless = 100 - $pctmore; $cntdeptless = $cntN - $cntdeptmore; }
    $cnttotN += $cntN; $sumtotcvg += $sumcvg; $cnttotdeptmore += $cntdeptmore;
    $smpl_hash{$smpl} = sprintf( "$smplmax\t%.1f\t$smplmin\t$cntN\t$cntdeptmore\t%.1f\t$cntdeptless\t%.1f",  $smplavgcvg, $pctmore, $pctless );
  }
} 

# -----------------------------------------------------------------
# write to ouptput
# -----------------------------------------------------------------
open ( OT, "> $outdir/$output.coverage_summary.xls" ) or print "\nFail to open $outdir/$output.coverage_summary.xls\n" and die; 
 
$outstrh  = "SAMPLE-BASED MAPPING, VARIANCE COUNT, COVERAGE, AND PHIX CONTROL SUMMARY TABLE\n";
$outstrh .= "Samples\tReadsNum\tReadsMap2Ref\tPct2ReadsNum\tReadsMap2Target\tPct2ReadsNum\ttotSNP\tSNP2Target\ttotDIP\tDIP2Target\t";
$outstrh .= "MaxCVG\tAvgCVG\tMinCVG\tTargetBasesWOphix\tBasesGrEqCvg$desiredCoverage\tGrEqCvgPct\tBasesLessCvg$desiredCoverage\tLessCvgPct\t";
$outstrh .= "PhixNum\tPct2ReadNum\tDupRate\tDupNum25%perMreads\tDupNum50%perMreads\tDupNum75%perMreads\tDupNum100%perMreads\tMapQuality0Pct\t";
$outstrh .= "TsNum\tTvNum\tTsTvRatio\tTsNumOnTarget\tTvNumOnTarget\tTsTvRatioOnTarget";

if ( $flag ){ 
  # special for NGS20-NONTP sample and sample-mm case
  $cntmm = 0;
  for ( $s = 0; $s < $cntsample; $s++ ){ 
    $key = $sampleArr[$s]; $keymm = "$key-mm";
    if ( exists $fastqsum_hash{$key} && exists $fastqsum_hash{$keymm} ){ 
      $ratio = $fastqsum_hash{$key} / ( $fastqsum_hash{$key} + $fastqsum_hash{$keymm} );
      $samplemmratio_hash{$key} = $ratio; $samplemmratio_hash{$keymm} = 1 - $ratio;
      $cntmm++;
      if ( exists $cvg1sum_hash{$key} && exists $cvg1sum_hash{$keymm} ){
        $ratio = $cvg1sum_hash{$key} / ( $fastqsum_hash{$key} + $fastqsum_hash{$keymm} );
        $samplemmratiotgt_hash{$key} = $ratio;
        $ratio = $cvg1sum_hash{$keymm} / ( $fastqsum_hash{$key} + $fastqsum_hash{$keymm} );
        $samplemmratiotgt_hash{$keymm} = $ratio;
      }
    }
  } 

  if ( $cntmm ){ $outstrh .= "\tsubSmplRatio\tsubSmplRatioTarget"; }

  print OT "$outstrh\n$outstr"; print $outstr; $cntchar = 0;
  for ( $s = 0; $s < $cntsample; $s++ ){ 
    $key = $sampleArr[$s]; $cntc = 0; $outstr = $key;
    if ( exists $fastqsum_hash{$key} && $fastqsum_hash{$key} != 0 ){ 
      $outstr .= "\t$fastqsum_hash{$key}"; $cntc += 1;

      if ( exists $cvg2sum_hash{$key} ){ 
        if ( $debug ) { $tmp1 = $fastqsum_hash{$key}; print "tmp1=$tmp1 - key=$key\n"; }
        if ( $cvg2sum_hash{$key} ){ $pct = 100 * $cvg2sum_hash{$key} / $fastqsum_hash{$key}; } else { $pct = 0; }
        $outstr .= sprintf ("\t$cvg2sum_hash{$key}\t%.1f", $pct );
      } else { $outstr .= "\t0\t0\t"; }  $cntc += 2;

      if ( exists $cvg1sum_hash{$key} ){ 
        if ( $cvg1sum_hash{$key} ) { $pct = 100 * $cvg1sum_hash{$key} / $fastqsum_hash{$key}; } else { $pct = 0; }
        # adjust read and pct if pct is over 100 due to double count reads due to close by bed regions
        if ( $pct > 100 ){ $pct = 100; $cvg1sum_hash{$key} = $cvg2sum_hash{$key}; }
        $outstr .= sprintf ("\t$cvg1sum_hash{$key}\t%.1f", $pct );
      } else { $outstr .= "\t0\t0"; } $cntc += 2;

      if ( exists $cntrawsnp_hash{$key} ){ $outstr .= "\t$cntrawsnp_hash{$key}"; } else { $outstr .= "\t0"; } $cntc += 1;
      if ( exists $cnttgtsnp_hash{$key} ){ $outstr .= "\t$cnttgtsnp_hash{$key}"; } else { $outstr .= "\t0"; } $cntc += 1;
      if ( exists $cntrawdip_hash{$key} ){ $outstr .= "\t$cntrawdip_hash{$key}"; } else { $outstr .= "\t0"; } $cntc += 1;
      if ( exists $cnttgtdip_hash{$key} ){ $outstr .= "\t$cnttgtdip_hash{$key}"; } else { $outstr .= "\t0"; } $cntc += 1;
      if ( exists $smpl_hash{$key}      ){ $outstr .= "\t$smpl_hash{$key}";      } else { $outstr .= "\t0\t0\t0\t0\t0\t0\t0\t0"; } $cntc += 8;
      
      if ( exists $phixsum_hash{$key} && $fastqsum_hash{$key} ){ 
        $phixpct = 100 * $phixsum_hash{$key}/$fastqsum_hash{$key}; 
        $outstr .= sprintf( "\t$phixsum_hash{$key}\t%.2f", $phixpct ); 
      } else { $outstr .= "\t0\t0"; } $cntc += 2;

      # for DuplicatRate
      if ( exists $pctdup_hash{$key} ){ $outstr .= "\t$pctdup_hash{$key}"; } else { $outstr.= "\t0";          } $cntc += 1;
      if ( exists $tgtdup_hash{$key} ){ $outstr .= "\t$tgtdup_hash{$key}"; } else { $outstr.= "\t0\t0\t0\t0"; } $cntc += 4;
      if ( exists $mq0_hash{$key}    ){ $outstr .= "\t$mq0_hash{$key}";    } else { $outstr.= "\t0";          } $cntc += 1;

      # for transition and transversion ratio
      if ( $cntrawsnp_hash{$key} > 0 ){ 
        $ts = $cntrawTS_hash{$key}; $tv = $cntrawsnp_hash{$key} - $ts; $ratio = 0; if ( $tv > 0 ){ $ratio = $ts/$tv; }
        $outstr .= sprintf( "\t$ts\t$tv\t%.1f", $ratio ); 
        $ts = $cnttgtTS_hash{$key}; $tv = $cnttgtsnp_hash{$key} - $ts; $ratio = 0; if ( $tv > 0 ){ $ratio = $ts/$tv; }
        $outstr .= sprintf( "\t$ts\t$tv\t%.1f", $ratio );  
      }
      else { $outstr .= "\t0\t0\t0\t0\t0\t0"; } $cntc += 6;

      # for mm sample if there is one
      if ( $cntmm ){ 
        if ( exists $samplemmratio_hash{$key}    ){ $outstr .= sprintf( "\t%.3f", $samplemmratio_hash{$key} );    } else { $outstr .= "\t0"; }
        if ( exists $samplemmratiotgt_hash{$key} ){ $outstr .= sprintf( "\t%.3f", $samplemmratiotgt_hash{$key} ); } else { $outstr .= "\t0"; }
      }
    }
    else { for ( $t = 0; $t < $cntchar; $t++ ){ $outstr .= "\t0"; }}
    print OT "$outstr\n"; print "$outstr\n";
    if ( $cntc > 5 ){ $cntchar = $cntc; }
  }
  $totavgcvg = 0; $totpassratio = 0; $cnttotdeptless = $cnttotN - $cnttotdeptmore;
  if ( $cnttotN ){ $totavgcvg = $sumtotcvg / $cnttotN; $totpasspct = 100 * $cnttotdeptmore / $cnttotN; $totlesspct = 100 - $totpasspct; } 
  $phixpct = 0; if ( $totfastqnum ){ $phixpct   = 100 * $totsumphix /$totfastqnum; } 
  $outstr  = "RunBatchTotal\t$totfastqnum\t\t\t\t\t\t\t\t\t\t"; 
  $outstr .= sprintf( "%.0f\t\t$cnttotN\t$cnttotdeptmore\t%.1f\t$cnttotdeptless\t%.1f\t$totsumphix\t%.1f\n", $totavgcvg, $totpasspct, $totlesspct, $phixpct );
  print OT $outstr; print $outstr;
}

# -----------------------------------------------------------------
# open $outdir/out/xxxx.cvg to summarize data on per gene base
# -----------------------------------------------------------------
if ( $technology eq "illumina" ) { 
  $cntspl = 0; print "\nSample:Lane --\>  ";
  foreach $sample (%lane_hash){if($lane_hash{$sample} ne ""){print "$lane_hash{$sample}:$sample "; $cntspl++;}} 
}

@inline = split ( /\s+/, $genelist ); $cntgene = 0; $outstr = "\nGene list --\> ";
for ( $i = 0; $i < scalar(@inline); $i++ ){
  if ( $inline[$i] ne "" ){
    $genelist_hash{$inline[$i]} = $inline[$i]; 
    $genearr[$cntgene] = $inline[$i]; 
    $cntgene++; 
    $outstr .= "$i:$inline[$i] ";
  }
} 
print " outstr = $outstr = \n"; $outstr .= "\n\n"; 

$outstrexonhead  = "\nSAMPLE AND GENEorAMPLICON-BASED SUMMARY TABLE\n";
$outstrexonhead .= "Sample\tGene\tFrag\tchr:start-stop Size\tMaxCVG\tAvgCVG\tMinCVG\tNote\tExonCvg2GeneCvg";

for ( $s = 0; $s < $cntsample; $s++ ) { 
  $corename = $sampleArr[$s]; 
  $filename = "$tmpdir/$corename.cvg";

  if ( $corename =~ /noBarcode/ )  {} 
  elsif ( open ( IN, "< $filename" ) ) {
    %gene_hash = (); $lane = $lane_hash{$corename}; $smpl = $corename; 
    $cntgeneN = 0; $cntdeptmore = 0; $smplmax = 0; $smplmin = 1000000;

    print "open $filename for reading corename=$corename\tlane=$lane\n "; 
    while ( <IN> ){ 
      # chr1    979403  979413  AGRN:NM_198576.3        In10/11 +       INTRON  9       2166
      # chr1    979403  979413  AGRN:NM_198576.3        In10/11 +       INTRON  10      2171
      # chr1    979478  979488  AGRN:NM_198576.3        In10/11 +       INTRON  1       2265
      # chr1    979478  979488  AGRN:NM_198576.3        In10/11 +       INTRON  2       2279
      @tp  = split(/\t/, $_); 
      @tp2 = split(/\:/, $tp[3]); 
      if ( exists $genelist_hash{$tp2[0]} ){ $gene_hash{$tp2[0]} .= $_; }
    } 
    close(IN);

    foreach $gene ( %gene_hash ) {
      if ( $gene_hash{$gene} ne "" ) {
        @inline = split ( /\n/, $gene_hash{$gene} );
        $cntexon = 0; %cvg_hash = (); # $max = 0; $min = 1000000; $sumcvgexon = 0;
        for ( $i = 0; $i < scalar(@inline); $i++ ) {
          @inline2 = split ( /\t/, $inline[$i] ); @inline3 = split ( /\:/, $inline2[3] ); # chr1  45805861  45805956  MutYH:NM_001128425  1  -   1  1083
          $exonc = $inline2[4]; $exon  = "$inline2[4] $inline2[0]:$inline2[1]-$inline2[2]"; 
          $cvg = $inline2[scalar(@inline2)-1]; # $cvg = $inline2[7]; # $allcvg += $cvg;
          if ( !exists $cvg_hash{$exon} ){ 
            $exonnum[$cntexon]  = $exon; 
            $exoncnum[$cntexon] = $exonc; 
            $cvg_hash{$exon} = $cvg; 
            $cntexon++; 
          } 
          else { 
            $cvg_hash{$exon} .= " $cvg"; 
          }
        }
        for ( $e = 0; $e < $cntexon; $e++ ) {
          $max = 0; $min = 1000000; $sumcvgexon = 0; $cntexonN = 0;
          @cvgline = split ( /\s/, $cvg_hash{$exonnum[$e]} );
          for ( $c = 0; $c < scalar(@cvgline); $c++ ) {
            if ( $cvgline[$c] ne "" ) {
              $sumcvgexon += $cvgline[$c]; $cntexonN++;
              if ( $max < $cvgline[$c] ) { $max = $cvgline[$c]; } if ( $min > $cvgline[$c] ) { $min = $cvgline[$c]; }
              if ( $cvgline[$c] >= $desiredCoverage && $inline3[0] ne "PHIX" ) { $cntdeptmore++; }
            }
          }
          $lanecvg_hash{$lane} += $sumcvgexon; $lanebase_hash{$lane} += $cntexonN;
          $genecvg_hash{$gene} += $sumcvgexon; $genebase_hash{$gene} += $cntexonN; $cntgeneN += $cntexonN;
          if ( $inline3[0] ne "PHIX" ) { 
            $smplcvg_hash{$smpl}  += $sumcvgexon; if ( $smplmax < $max ) { $smplmax = $max; } 
            $smplbase_hash{$smpl} += $cntexonN;   if ( $smplmin > $min ) { $smplmin = $min; }
            if ( $cntexonN ) { $exonavg = $sumcvgexon / $cntexonN; } else { $exonavg = 0; }
            if ( $min < 100 ) { $min .= "\tLess"; } else { $min .= "\t"; }
            $exongene = "$exonnum[$e]\t$gene"; 
            if ( $debug ){ print "corename=$corename exongene=$exongene end\n"; }
            $cvgpergene = ""; if ( exists $cvgexongene_hash{$corename}{$exongene} ){ $cvgpergene = $cvgexongene_hash{$corename}{$exongene}; }
            if ( $cvgpergene =~ /AvgLow/ ) { $cntAvgLow++; }
            @tpe = split ( /\s/, $exonnum[$e] );
            $outstrexon .= sprintf ( "$corename\t$gene\t$exoncnum[$e]\t$tpe[1] $cntexonN\t$max\t%.0f\t$min\t$cvgpergene\n", $exonavg );
            $sumgenecvg += $sumcvgexon;
          }
        }
        $geneexonnum_hash{$gene} = $cntexon; $genelen_hash{$gene} = $cntgeneN;
      }
    }
    $smpllessdept_hash{$smpl} = $cntdeptmore;
    $smpldeptmoreratio = 0; $smplAvgCvg = 0;
    if ( $smplbase_hash{$smpl} ) { $smpldeptmoreratio = 100 * $cntdeptmore / $smplbase_hash{$smpl}; $smplAvgCvg = $smplcvg_hash{$smpl} / $smplbase_hash{$smpl}; }
    $cntdeptless = $smplbase_hash{$smpl} - $cntdeptmore; $cntdeptlessratio = 100 - $smpldeptmoreratio; 
  }
}

$outstrg1 = "\nGENE-BASED COVERAGE TABLE\nGene"; $outstrg2 = "FlagNum"; $outstrg3 = "CVG";
for ( $g = 0; $g < $cntgene; $g++ ) {
  $gene = $genearr[$g]; $outstrg1 .= "\t$gene";
  $outstrg2 .= "\t$geneexonnum_hash{$gene}";
  if ( $genebase_hash{$gene} ) { $genecvg = $genecvg_hash{$gene} / $genebase_hash{$gene}; } else { $genecvg = 0; }
  $outstrg3 .= sprintf ("\t%.0f", $genecvg );
}

$lanetot = 0; $crosslaneflag = 1; 
for ( $l = 1; $l <= 8; $l++ ) { 
  $lane = "L00$l"; $lanetot += $lanecvg_hash{$lane}; 
  if ( exists $cntlane_hash{$corename} && $cntlane_hash{$corename} > 2 ) { $crosslaneflag = 0; }
}
if ( $lanetot ) {  
  $outstrlane  = "LANE-BASED COVERAGE SUMMARY TABLE\nLANE"; for ( $l = 1; $l <= 8; $l++ ) { $outstrlane .= "\tlane$l"; }
  $outstrlane .= "\nCVG";    
  if ( $crosslaneflag ) { 
    for ( $l = 1; $l <= 8; $l++ ) {
      $lane = "L00$l";
      if ( $lanebase_hash{$lane} ) { $lanecvg = int ($lanecvg_hash{$lane} / $lanebase_hash{$lane}); } else { $lanecvg = 0; }
      $outstrlane .= sprintf ( "\t%.0f", $lanecvg );
    }
    $outstrlane .= "\n"; 
  }
  else { $outstrlane .= "\tNA\n"; }
}

print "\n$outstrlane\n"; print "$outstrg1\n$outstrg2\n$outstrg3\n"; 
if ( $outstrlane ne "" ) { print OT "\n$outstrlane\n"; }
print OT "$outstrg1\n$outstrg2\n$outstrg3\n";
$outstrexonhead .= "\tcntAvgLow=$cntAvgLow";
print OT "$outstrexonhead\n$outstrexon";
close ( OT );
