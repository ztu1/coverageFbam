#!/usr/bin/perl
#
###@StandardDocs
# This script take indir cotains bam file and call picard tools for duplication
# call to generate Metrics file
#
# HISTORY
#   -- 20151017: Zheng Jin Tu: initialize
#

$Prg = "dup_via_picard.pl";

sub usage { print "\nUsage: $Prg ";
            print "\n\t\t -id indir        # REQUIRED input directory";
            print "\n\t\t  -o output       # REQUIRED output";
            print "\n\t\t -pp picardpath   # REQUIRED /path/to/picard tools;";
            print "\n\t\t -rj runjob       # run job, DEFAULT no";
            print "\n\t\t -dg debug        # for development use DEFAULT no]\n\n"; }

use File::Basename;
use Getopt::Long;

if ($#ARGV < 2){ warn "number of commandline args: $#ARGV\n" if $debug; &usage; exit 1; }

@Options =
(
    "id=s"  => \$indir,
    "o=s"   => \$output,
    "pp=s"  => \$picardpath,
    "rj=s"  => \$runjob,
    "dg=s"  => \$debug
);

if (!&GetOptions(@Options) ){ &usage; exit 1;                                                   }
if (!$indir      ){ warn "\nMust specify an input directory \n"; &usage; exit 1;                }
if (!$output     ){ warn "\nMust specify an output name     \n"; &usage; exit 1;                }
if (!$picardpath ){ warn "\nMust specify path to picard     \n"; &usage; exit 1;                }
if (!$runjob     ){ $runjob     = "No";                                                         }
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
    $exec  = " java -jar $picardpath/picard.jar SortSam SORT_ORDER=coordinate INPUT=$indir/\"$dir_list\" OUTPUT=$output.sorted.bam\n";
    $exec .= " sleep 1\n";
    $exec .= " java -jar $picardpath/picard.jar MarkDuplicates REMOVE_DUPLICATES=true METRICS_FILE=$output INPUT=$output.sorted.bam OUTPUT=$output.d.bam\n";
    $exec .= " sleep 1\n";
    if ( !$debug ){ $exec .= " rm -fr $output.sorted.bam $output.d.bam\n"; }
    if ( $runjob eq "Yes" || $runjob eq "yes" ){ system( $exec ); } else { print "$cntbam\n$exec"; }
    $cntbam++; 
  }
}
closedir DH;

if ( $cntbam != 1 ){ print "\nThere is either no or more than one bam files detected in $indir\n"; }
