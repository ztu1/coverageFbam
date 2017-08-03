#!/usr/bin/perl
#
###@StandardDOcs
#  This perl script takes inputs: bam and bedfile
#  then return mapping quality zero ratio and duplicate
#  read number
# 
#  INPUT/OUTPUT
#    -i input         # REQUIRED input sorted bam file
#    -b bedfile       # bed file to define target region to be looked
#    -o output        # output name
#    -s samtools      # samtolls program defualt /usr/local/biotools/samtools/0.1.18/samtools
#   -el normalizeread # normalizeread read number default 1000000
#   -dg debug]        # default 0
#
# HISTORY:
#  20150926	Zheng Jin Tu	initialize this script
#  20151108	Zheng Jin Tu	re-code couting, allow whole genome, and fix normalization issue

$Prg = "dupNum_mq0_bam.pl";
 
sub usage { print "\nUsage: $Prg ";
            print "\n\t\t  -i input         # REQUIRED input sorted bam file";
            print "\n\t\t  -b bedfile       # bed file to define target region to be looked";
            print "\n\t\t  -o output        # output name" ;
            print "\n\t\t  -s samtools      # samtools default /usr/local/biotools/samtools/0.1.18/samtools";
            print "\n\t\t -nr normalizeread # extended length over bed region default 150";
            print "\n\t\t -dg debug]        # default 0\n"; }

use File::Basename;
use Getopt::Long;
#use strict;

if ($#ARGV < 1){ warn "number of commandline args: $#ARGV\n" if $debug; &usage; exit 1; }

@Options =
(
    "i=s"   => \$input,
    "b=s"   => \$bedfile,
    "o=s"   => \$output,
    "s=s"   => \$samtools,
    "nr=s"  => \$normread,
    "dg=s"  => \$debug
);

if (!&GetOptions(@Options) ){ &usage; exit 1;                                        }
if (!$input       ){ warn "\nMust specify an input file \n"; &usage; exit 1;         } $inputc = &lastName($input);
#f (!$bedfile     ){ warn "\nMust specify a bed filea   \n"; &usage; exit 1;         }
if (!$output      ){ $output   = &lastName($input); $output =~ s/bam$/info/;         }
if (!$samtools    ){ $samtools = "/usr/local/biotools/samtools/0.1.18/samtools";     }
if (!$normread    ){ $normread = 1000000;                                            }
if (!$debug       ){ $debug    = 0;                                                  }

# -------------------------------------------------------------------------------------------
# read $bedfile and mark region to be checked
# -------------------------------------------------------------------------------------------
$cntline = 0; $cntpos = 0;
if ( open ( IN, "< $bedfile" )){  
  while ( <IN> ) {
    chomp; @tp = split ( /\t/, $_ ); $chr = $tp[0]; 
    @tp1 = split ( /\:/, $tp[3] ); $gene = $tp1[0];
    if ( $tp[1] < $tp[2] ){ $start = $tp[1]; $stop = $tp[2]; } else { $start = $tp[2]; $stop = $tp[1]; }
    for ( $i = $start; $i <= $stop; $i++ ){ 
      $dna_hash{$chr}{$i} = $tp[3]; $posArr[$cntpos] = "$chr:$i"; $cntpos++;
    }
    $cntline++;
  }
  close ( IN );
  print  "\n$bedfile cntline=$cntline\tcntpos=$cntpos\n";
}

# -------------------------------------------------------------------------------------------
# scan $input bam file 
# -------------------------------------------------------------------------------------------
open IN, "$samtools view $input |" or die "Can't open $samtools pipe: $!\n"; 
print "Read $input data ...\n";
$cntread = 0; $cntreadontarget = 0; $cntmq0 = 0; $totlen = 0; $recloc = 0; $recchr = ""; $maxdup = 0;

if ( $bedfile && $cntpos ){ # if bedfile identified, only looking for target region
  while ( <IN> ) {
    # R0239954:8:000000000-A4EHL:1:2105:18276:3467_1:N:0:TGACCA 99 chrM 8662 60 107M1I5M38S = 9074 563 CAAT  /GFA< RG:Z:15d3ec13 NH:i: 1
    chomp; @tp = split( /\t/, $_ ); $chr = $tp[2]; $loc = $tp[3]; $cigar = $tp[5];
    if ( !/^\@/ && !/\^\#/ && $chr ne "chrPhix" ){ 
      if ( exists $dna_hash{$chr}{$loc} ){ 
        if ( $recloc != $loc || $chr ne $recchr ){ 
          foreach $key ( %dup_hash ){ 
            if ( $dup_hash{$key} ne "" ){ $dupArr[$dup_hash{$key}]++; }
            if ( $maxdup < $dup_hash{$key} ){ $maxdup = $dup_hash{$key}; $maxchrloc = "$recchr:$recloc"; }
          }
          # reset variables
          %dup_hash = (); $recchr = $chr; $recloc = $loc;
        }
        if ( $cigar ne "*" ){ $dup_hash{$cigar}++; }
        $totlen += length( $tp[9] );
        $cntreadontarget++;
      }
      if ( $tp[4] == 0 ){ $cntmq0++; }
      if ( $cntread % 100000 == 1 ){ print "  read=$cntread\tontgt=$cntreadontarget\tmq0=$cntmq0\tmaxdup=$maxdup\tmaxchrloc=$maxchrloc\n"; } 
      $cntread++;
    }
  }    
  if ( $cntreadontarget ){ $avgReadLen = int( $totlen / $cntreadontarget ); }
  $readused = $cntreadontarget;
  print "$input read=$cntread\tcntreadontarget=$cntreadontarget mq0=$cntmq0\treadused=$readused\tavgReadLen=$avgReadLen\n\n"; 
}
else { # for the whole genome
  while ( <IN> ) {
    # R0239954:8:000000000-A4EHL:1:2105:18276:3467_1:N:0:TGACCA 99 chrM 8662 60 107M1I5M38S = 9074 563 CAAT  /GFA< RG:Z:15d3ec13 NH:i: 1
    chomp; @tp = split( /\t/, $_ ); $chr = $tp[2]; $loc = $tp[3]; $cigar = $tp[5];
    if ( !/^\@/ && !/\^\#/ & $chr ne "chrPhix" ){
      if ( $recloc != $loc ){
        foreach $key ( %dup_hash ){
          if ( $dup_hash{$key} ne "" ){
            $dupArr[$dup_hash{$key}]++; 
            if ( $maxdup < $dup_hash{$key} ){ $maxdup = $dup_hash{$key}; $maxchrloc = "$recchr:$recloc"; }
          }
        }
        # reset variables
        %dup_hash = (); $recchr = $chr; $recloc = $loc;
      }
      if ( $cigar ne "*" ){ $dup_hash{$cigar}++; }
      $totlen += length( $tp[9] );
      if ( $tp[4] == 0 ){ $cntmq0++; }
      if ( $cntread % 100000 == 1 ){ print "  read=$cntread\tmq0=$cntmq0\tmaxdup=$maxdup\tmaxchrloc=$maxchrloc\n"; }
      $cntread++;
    }
  }
  if ( $cntread ){ $avgReadLen = int( $totlen / $cntread ); }
  $readused = $cntread;
  print "$input read=$cntread\tmq0=$cntmq0\treadused=$readused\tavgReadLen=$avgReadLen\n\n";
}

close ( IN ); 

# ----------------------------------------------------------------------------------------------
# calcuate numbers and write output
# ----------------------------------------------------------------------------------------------
$normfact = $cntread / $normread;
if ( $readused && $normfact ){
  $mq0ratio = 100 * $cntmq0 / $cntread;
  $sumdupnum = 0; $pct25flag = 1; $dup25 = 0; $pct50flag = 1; $dup50 = 0; $pct75flag = 1; $dup75 =0;
  for ( $i = 0; $i <= $maxdup; $i++ ){
    if ( $dupArr[$i] ){ 
      $sumdupnum += $dupArr[$i] * $i; $ratio = $sumdupnum / $readused;
      if ( $ratio >= 0.25 && $pct25flag ){ $dup25 = $i / $normfact; $pct25flag = 0; } 
      if ( $ratio >= 0.50 && $pct50flag ){ $dup50 = $i / $normfact; $pct50flag = 0; } 
      if ( $ratio >= 0.75 && $pct75flag ){ $dup75 = $i / $normfact; $pct75flag = 0; } 

#     if ( $debug ){ print "  $i\t$dupArr[$i]\t$ratio\n"; }
      elsif ( $i < 10 || $i % 100 == 1 ){ print "  $i\t$dupArr[$i]\t$ratio\n"; }
    }
    if ( $debug ){ print "$i\tdupArr=$dupArr[$i]\tratio=$ratio\tpct25flag=$pct25flag\tpct50flag=$pct50flag\tpct75flag=$pct75flag\n"; }
  }

  if ( $dup50 < $dup25 ){ $dup50 = $dup25; } if ( $dup75 < $dup50 ){ $dup75 = $dup50; }

  $dup100 = $maxdup / $normfact; 
  print "cntread=$cntread\tcntreadontarget=$cntreadontarget\tnormfact=$normfact\n";
  $outstr  = sprintf ( "%.0f\t%.0f\t%.0f\t%.0f\t%.3f\t", $dup25, $dup50, $dup75, $dup100, $mq0ratio );
  $outstr .= "$input\t$bedfile\t$maxdup\t$maxchrloc\n";
}
else { $outstr = "\t\t\t\t\t$input\t$bedfile\t\n"; }

# ----------------------------------------------------------------------------------------------
# write output file
# ----------------------------------------------------------------------------------------------
open ( OT, "> $output" ) or print "\nFail to open $output\n\n" and die;
print OT $outstr; print $outstr;
close ( OT );

# ---------------------------------------------------------------------------------------------
# sub function
# ---------------------------------------------------------------------------------------------
sub lastName {($sub_input) = @_; @sub_inline = split(/\//, $sub_input); $sub_fname = $sub_inline[scalar(@sub_inline)-1]; return $sub_fname;}
