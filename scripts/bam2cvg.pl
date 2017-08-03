#!/usr/bin/perl
#
###@StandardDocs
# This script take indir cotains bam file and a bed file.
# It will call bedtools to generate coverage file
#
# HISTORY
#   -- 20151014: Zheng Jin Tu: initialize
#

$Prg = "dlmp_bam2cvg.pl";

sub usage { print "\nUsage: $Prg ";
            print "\n\t\t -id indir        # REQUIRED input directory";
            print "\n\t\t  -o output       # REQUIRED output";
            print "\n\t\t  -b bedfile      # REQUIRED bed file for target";
            print "\n\t\t -bt bedtools     # REQUIRED path to bedtools";
            print "\n\t\t -bp bedtoolspara # bedtools parameter default \"\" other choice \"-hist -d\";";
            print "\n\t\t -rj runjob       # run job, DEFAULT yes";
            print "\n\t\t -dg debug        # for debug use DEFAULT no]\n\n"; }

use File::Basename;
use Getopt::Long;

if ($#ARGV < 4){ warn "number of commandline args: $#ARGV\n" if $debug; &usage; exit 1; }

@Options =
(
    "id=s"  => \$indir,
    "o=s"   => \$output,
    "b=s"   => \$bedfile,
    "bt=s"  => \$bedtools,
    "bp=s"  => \$bdtlspara,
    "rj=s"  => \$runjob,
    "dg=s"  => \$debug
);

if (!&GetOptions(@Options) ){ &usage; exit 1;                                                   }
if (!$indir      ){ warn "\nMust specify an input directory \n"; &usage; exit 1;                }
if (!$output     ){ warn "\nMust specify an output name     \n"; &usage; exit 1;                }
if (!$bedfile    ){ warn "\nMust specify bedfile            \n"; &usage; exit 1;                }
if (!$bedtools   ){ warn "\nMust specify path/to/bedtools   \n"; &usage; exit 1;                }
if (!$bdtlspara  ){ $bdtlspar   = "";                                                           }
if (!$runjob     ){ $runjob     = "yes";                                                        }
if (!$debug      ){ $debug      = 0;                                                            }
  
# -------------------------------------------------------------------------------------------
# open $indir to find bam file for bedtools
# -------------------------------------------------------------------------------------------
$cntbam = 0;
opendir DH, $indir or print "\nFail to open $indir\n\n" and die;
@dir_list = grep /\.bam$/, readdir(DH); 
foreach $dir_list ( @dir_list ){ 
  if ( $dir_list =~ /extended/ ){}
  else { 
    $exec  = " $bedtools coverage -abam $indir/\"$dir_list\" -b $bedfile $bdtlspara \> $output\n";
    if ( $runjob eq "Yes" || $runjob eq "yes" ){ system( $exec ); } else { print "$cntbam\n$exec"; }
    $cntbam++; 
  }
}
closedir DH;

if ( $cntbam != 1 ){ print "\nThere is either no or more than one bam files detected in $indir\n"; }
