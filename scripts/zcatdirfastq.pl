#!/usr/bin/perl
#
###@StandardDocs
# This script take input run sample directory cotains fastq.gz file(s) 
# and call zcat to count fastq sequence number 
#
# HISTORY
#   -- 20151014: Zheng Jin Tu: initialize  
#

$Prg = "dlmp_zcatdirfastq.pl";

sub usage { print "\nUsage: $Prg ";
            print "\n\t\t -id indir        # REQUIRED input directory";
            print "\n\t\t  -o output       # REQUIRED output name";
            print "\n\t\t -rj runjob       # run job, DEFAULT yes";
            print "\n\t\t -dg debug        # for development use DEFAULT no\n\n"; }

use File::Basename;
use Getopt::Long;

if ($#ARGV < 2){ warn "number of commandline args: $#ARGV\n" if $debug; &usage; exit 1; }

@Options =
(
    "id=s"  => \$indir,
    "o=s"   => \$output,
    "rj=s"  => \$runjob,
    "dg=s"  => \$debug
);

if (!&GetOptions(@Options) ){ &usage; exit 1;                                                   }
if (!$indir      ){ warn "\nMust specify an input directory \n"; &usage; exit 1;                }
if (!$output     ){ warn "\nMust specify an output          \n"; &usage; exit 1;                }
if (!$runjob     ){ $runjob     = "yes";                                                        }
if (!$debug      ){ $debug      = 0;                                                            }
  
# -------------------------------------------------------------------------------------------
# open $indir to find fastq.gz files
# -------------------------------------------------------------------------------------------
$cntfastq = 0; $exec = "rm -fr $output\n";
opendir DH, $indir or print "\nFail to open $indir\n\n" and die;
@dir_list = grep /\.fastq.gz$/, readdir(DH); 
foreach $dir_list ( @dir_list ){ 
  if ( $dir_list =~ /fastq.gz$/ ) {
    $exec .= " zcat $indir/$dir_list | wc -l \>\> $output\n";
    $cntfastq++;
  }
}
closedir DH;

if ( $cntfastq ){ 
  if ( $runjob eq "Yes" || $runjob eq "yes" ){ system( $exec ); } else { print "$cntfastq\n$exec\n"; }
}
else { print "\nThere is no fastq.gz detected in $indir\n"; }
