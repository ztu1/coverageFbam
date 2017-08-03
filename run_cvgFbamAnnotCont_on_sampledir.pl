#!/usr/bin/perl
#
###@StandardDocs 
#  This perl script takes input sample run directory with subdir contains 
#  bam file(s), a panel specific configure file and call BEDTOOLS, zcat, dup scripts 
#  for each sample amd the summarize data and output coverage_summary.xls file
#  
#  INPUT: 
#                -id indir          # REQUIRED input top directory
#                 -c configfile     # REQUIRED panel specific configure file
#                -od outdir         # REQUIRED outdir name
#                -ot output         # REQUIRED output name suggest runbatch name as output name
#               [-ss samplesheet    # OPTION Illumina SampleSheet.csv or Ion Torrent SampleSheet.txt file
#                 -t technology     # CHOICE: illumina, iontorrent; DEFAULT illumina
#                 -l label          # Torrent sample name (I006-12345678) has label CHOICE: yes no, DEFAULT yes
#                -td tmpdir         # tmpdir for job files default outdir
#                -rb runBEDTOOLS    # run BEDTOOLS DEFAULT yes
#                -rd runduplication # run deplication DEFAULT yes
#                -rz runzcat        # run zcat for read num DEFAULT yes
#                -rj runjob         # run job, DEFAULT no
#                 -l logfile        # logfile
#                -dg debug          # for development use DEFAULT no]
#  OUTPUT
#    file with multiple sessions with coverage per sample, per lane, pergene, and per gene and exon 
#
# HISOTY
#  -- 20150926 Zheng Jin Tu:   add dup and mq0
#  -- 20151014 Zheng Jin Tu:   add contamination check portion
#  -- 20160211 Zheng Jin Tu:   rename it from dir2.dlmp_cvgFbamAnnotCont.pl
#  -- 20170214 Zheng Jin Tu:   add mergebed for smpl.merged.cvg1 to deal closeby cava-style target.bed
#                              remove bedfilectm code since it is not use at this package
#
my $Prg = "run_cvgFbamAnnotCont_on_sampledir.pl";

sub usage { print "\nUsage: $Prg ";
            print "\n\t\t -id indir          # REQUIRED input top directory";
            print "\n\t\t  -c configfile     # REQUIRED panel specific configure file";
            print "\n\t\t -od outdir         # REQUIRED outdir name";
            print "\n\t\t -ot output         # REQUIRED output name suggest runbatch name as output name";
            print "\n\t\t[-ss samplesheet    # OPTION Illumina SampleSheet.csv or Ion Torrent SampleSheet.txt file";
            print "\n\t\t  -t technology     # CHOICE: illumina, iontorrent; DEFAULT illumina";
            print "\n\t\t  -l label          # Torrent sample name (I006-12345678) has label CHOICE: yes no, DEFAULT yes";
            print "\n\t\t -td tmpdir         # tmpdir for job files default outdir";
            print "\n\t\t -rb runBEDTOOLS    # run BEDTOOLS DEFAULT yes";
            print "\n\t\t -rd runduplication # run deplication DEFAULT yes";
            print "\n\t\t -rz runzcat        # run zcat for read num DEFAULT yes";
            print "\n\t\t -rj runjob         # run job, DEFAULT no";
            print "\n\t\t -ld logdir         # logdir to hold log file";
            print "\n\t\t -dg debug          # for development use DEFAULT no]\n\n"; }

use File::Basename;
use Getopt::Long;
use Data::Dumper;
use strict;
use warnings;

my ( $indir, $configfile, $outdir, $output, $samplesheet, $technology, $label, $tmpdir, $runbdtl, $runzcat, $rundup );
my ( $logdir, $runjob, $debug );

if ($#ARGV < 4) { warn "number of commandline args: $#ARGV\n" if $debug; &usage; exit 1; }

my @Options =
(
    "id=s"  => \$indir,
    "c=s"   => \$configfile,
    "od=s"  => \$outdir,
    "ot=s"  => \$output,
    "ss=s"  => \$samplesheet,
    "t=s"   => \$technology,
    "l=s"   => \$label,
    "td=s"  => \$tmpdir,
    "rb=s"  => \$runbdtl,
    "rz=s"  => \$runzcat,
    "rd=s"  => \$rundup,
    "ld=s"  => \$logdir,
    "rj=s"  => \$runjob,
    "dg=s"  => \$debug
);

my $SCRIPT_HOME = $ENV{'SCRIPT_HOME'};
my $SAMTOOLS    = $ENV{'SAMTOOLS'};
my $BEDTOOLS    = $ENV{'BEDTOOLS'};
my $PICARD      = $ENV{'PICARD'};
my $JOBLIMIT    = $ENV{'JOBLIMIT'};
my $SLEEPTIME   = $ENV{'SLEEPTIME'};
my $SBATCH      = $ENV{'SBATCH'};
my $SQUEUE      = $ENV{'SQUEUE'};
my $TMPDIR      = $ENV{'TMPDIR'};
my $T           = &gettimestamp();

if (!&GetOptions(@Options)){ &usage; exit 1;                                                      }
if (!$indir       ){ warn "\nMust specify an input directory   \n"; &usage; exit 1;               }
if (!$configfile  ){ warn "\nMust specify panel configfile file\n"; &usage; exit 1;               }
if (!$outdir      ){ warn "\nMust specify outdir name          \n"; &usage; exit 1;               }
if (!$output      ){ $output      = "$Prg"; $output =~ s/\.pl//;                                  }
if (!$tmpdir      ){ $tmpdir      = "$TMPDIR/$Prg.$T";                                         }
if (!$technology  ){ $technology  = "illumina";                                                   }
if (!$label       ){ $label       = "yes";                                                        }
if (!$runbdtl     ){ $runbdtl     = "yes";                                                        }
if (!$runzcat     ){ $runzcat     = "yes";                                                        }
if (!$rundup      ){ $rundup      = "yes";                                                        }
if (!$logdir      ){ $logdir      = "$outdir/log";                                                }
if (!$runjob      ){ $runjob      = "no";                                                         }
if (!$debug       ){ $debug       = 0;                                                            }
  
# --------------------------------------------------------------------------------
# first read in configure file
# --------------------------------------------------------------------------------
my ( $cvgTargetBed, $cvgGenomeBed, $zerobasebed );
open ( IN, "< $configfile" ) or print "\nFail to open $configfile\n\n";
while ( <IN> ) {
  if ( !/^\#/ ){
    chomp; my @tp = split ( /\"/, $_ );
    if    ( /^cvgTargetBed=/              ){ $cvgTargetBed               = $tp[1]; }
    elsif ( /^cvgGenomeBed=/              ){ $cvgGenomeBed               = $tp[1]; }
  }
}
close ( IN );

if ( !$cvgGenomeBed ){ print "\nERROR: can not identify cvgGenomeBed in $configfile\n"; exit 1; }
if ( !$cvgTargetBed ){ print "\nERROR: can not identify cvgTargetBed in $configfile\n"; exit 1; }

my $exec;
if ( -d $outdir ){} else { $exec = " mkdir $outdir\n"; system ( $exec ); print $exec; }
if ( -d $logdir ){} else { $exec = " mkdir $logdir\n"; system ( $exec ); print $exec; }
if ( -d $tmpdir ){} else { $exec = " mkdir $tmpdir\n"; system ( $exec ); print $exec; }

open ( LG, "> $logdir/$Prg.$T.log" ) or print "\nFail to open $logdir/$Prg.$T.log\n\n" and die;
print LG "$T\n";

# merge $cvgTargetBed to eliminate possible overlap inside cvgTargetBed
$exec = "$BEDTOOLS/bedtools sort -i $cvgTargetBed \| $BEDTOOLS/bedtools merge -i stdin \> $tmpdir/merged.sorted.bed\n";
system ( $exec ); print LG $exec;

my ( $JOBHEAD, $jobname, $jobfile, $jobid, $jobnum, @dir_list, $dir_list ); 
my $cntdir = 0;
my $cntjob = 0;
my $jobstr = "";
my $JOBSLEEPTIME = 10 * $SLEEPTIME;

opendir DH, $indir or print "\n\tFail to open $indir\n" and die;
@dir_list = grep !/^\./, readdir(DH);
foreach $dir_list ( @dir_list ){
  $cntdir++;
  if ( $dir_list =~ /log/ ){}
  else {
    # create and submit job: cvg on every base 
    $exec = " $SCRIPT_HOME/scripts/bam2cvg.pl -bt $BEDTOOLS/bedtools -id $indir/$dir_list -o $tmpdir/$dir_list.cvg -b $cvgTargetBed -bp \"-hist -d\"\n";

    $jobname = "cvg.$dir_list"; $jobfile = "$logdir/$jobname.$T.job"; $cntjob++;
    $JOBHEAD = "\#!\/bin\/sh\n\n\#SBATCH -J $jobname\n\#SBATCH -e $tmpdir/$jobname.$T.err\n\#SBATCH -o $tmpdir/$jobname.$T.out\n\n";
    open ( JB, "> $jobfile" ) or print "\nFail to open $jobfile\n" and die; print JB "$JOBHEAD\n$exec"; close ( JB );
    if ( $runbdtl eq "yes" && $runjob eq "yes" ){ 
      do { $jobnum = &checkUserJob(); if ( $jobnum > $JOBLIMIT ){ sleep $JOBSLEEPTIME; print "$jobnum "; }} while ( $jobnum > $JOBLIMIT );
      $jobid = `$SBATCH $jobfile`; chomp( $jobid ); print "jobid=$jobid end\n"; $jobid =~ s/Submitted batch job //; $jobstr .= "$jobid ";
      print " $cntjob submit $jobid\t$jobfile\n"; sleep $SLEEPTIME; 
    } else { print " $cntjob skipped\t$jobfile\n"; } 
    print LG $exec;

    # create and submit job: cvg on bed fragment 
    $exec = " $SCRIPT_HOME/scripts/bam2cvg.pl -bt $BEDTOOLS/bedtools -id $indir/$dir_list -o $tmpdir/$dir_list.cvg1 -b $cvgTargetBed\n";
    $jobname = "cvg1.$dir_list"; $jobfile = "$logdir/$jobname.$T.job"; $cntjob++;
    $JOBHEAD = "\#!\/bin\/sh\n\n\#SBATCH -N 1\n\#SBATCH -J $jobname\n\#SBATCH -e $tmpdir/$jobname.$T.err\n\#SBATCH -o $tmpdir/$jobname.$T.out\n\n";
    open ( JB, "> $jobfile" ) or print "\nFail to open $jobfile\n" and die; print JB "$JOBHEAD\n$exec"; close ( JB );

    if ( $runbdtl eq "yes" && $runjob eq "yes" ){
      do { $jobnum = &checkUserJob(); if ( $jobnum > $JOBLIMIT ){ sleep $JOBSLEEPTIME; print "$jobnum "; }} while ( $jobnum > $JOBLIMIT );
      $jobid = `$SBATCH $jobfile`; chomp( $jobid ); $jobid =~ s/Submitted batch job //; $jobstr .= "$jobid ";
      print " $cntjob submit $jobid\t$jobfile\n"; sleep $SLEEPTIME;
    } else { print " $cntjob skipped\t$jobfile\n"; }
    print LG $exec;

    # create and submit job: cvg on merged.bed fragment 
    $exec = " $SCRIPT_HOME/scripts/bam2cvg.pl -bt $BEDTOOLS/bedtools -id $indir/$dir_list -o $tmpdir/$dir_list.merged.cvg1 -b $tmpdir/merged.sorted.bed\n";
    $jobname = "cvgm.$dir_list"; $jobfile = "$logdir/$jobname.$T.job"; $cntjob++;
    $JOBHEAD = "\#!\/bin\/sh\n\n\#SBATCH -N 1\n\#SBATCH -J $jobname\n\#SBATCH -e $tmpdir/$jobname.$T.err\n\#SBATCH -o $tmpdir/$jobname.$T.out\n\n";
    open ( JB, "> $jobfile" ) or print "\nFail to open $jobfile\n" and die; print JB "$JOBHEAD\n$exec"; close ( JB );
    if ( $runbdtl eq "yes" && $runjob eq "yes" ){ 
      do { $jobnum = &checkUserJob(); if ( $jobnum > $JOBLIMIT ){ sleep $JOBSLEEPTIME; print "$jobnum "; }} while ( $jobnum > $JOBLIMIT );
      $jobid = `$SBATCH $jobfile`; chomp( $jobid ); $jobid =~ s/Submitted batch job //; $jobstr .= "$jobid ";
      print " $cntjob submit $jobid\t$jobfile\n"; sleep $SLEEPTIME;
    } else { print " $cntjob skipped\t$jobfile\n"; }
    print LG $exec;

    # create and submit job: cvg on reference 
    $exec = " $SCRIPT_HOME/scripts/bam2cvg.pl -bt $BEDTOOLS/bedtools -id $indir/$dir_list -o $tmpdir/$dir_list.cvg2 -b $cvgGenomeBed\n";
    $jobname = "cvg2.$dir_list"; $jobfile = "$logdir/$jobname.$T.job"; $cntjob++;
    $JOBHEAD = "\#!\/bin\/sh\n\n\#SBATCH -N 1\n\#SBATCH -J $jobname\n\#SBATCH -e $tmpdir/$jobname.$T.err\n\#SBATCH -o $tmpdir/$jobname.$T.out\n\n";
    open ( JB, "> $jobfile" ) or print "\nFail to open $jobfile\n" and die; print JB "$JOBHEAD\n$exec"; close ( JB );
    if ( $runbdtl eq "yes" && $runjob eq "yes" ){
      do { $jobnum = &checkUserJob(); if ( $jobnum > $JOBLIMIT ){ sleep $JOBSLEEPTIME; print "$jobnum "; }} while ( $jobnum > $JOBLIMIT );
      $jobid = `$SBATCH $jobfile`; chomp( $jobid ); $jobid =~ s/Submitted batch job //; $jobstr .= "$jobid ";
      print " $cntjob submit $jobid\t$jobfile\n"; sleep $SLEEPTIME;
    } else { print " $cntjob skipped\t$jobfile\n"; }
    print LG $exec;

    # create and submit job: map quality 0 and duplication 
    $exec = " $SCRIPT_HOME/scripts/dupNum_mq0_bam.pl -i $indir/$dir_list/$dir_list.bam -o $tmpdir/$dir_list.dup.mq0 -b $cvgTargetBed -s $SAMTOOLS/samtools\n";
    $jobname = "mq0.$dir_list"; $jobfile = "$logdir/$jobname.$T.job"; $cntjob++;
    $JOBHEAD = "\#!\/bin\/sh\n\n\#SBATCH -N 1\n\#SBATCH -J $jobname\n\#SBATCH -e $tmpdir/$jobname.$T.err\n\#SBATCH -o $tmpdir/$jobname.$T.out\n\n";
    open ( JB, "> $jobfile" ) or print "\nFail to open $jobfile\n" and die; print JB "$JOBHEAD\n$exec"; close ( JB );
    if ( $rundup eq "yes" && $runjob eq "yes" ){ 
      do { $jobnum = &checkUserJob(); if ( $jobnum > $JOBLIMIT ){ sleep $JOBSLEEPTIME; print "$jobnum "; }} while ( $jobnum > $JOBLIMIT );
      $jobid = `$SBATCH $jobfile`; chomp( $jobid ); $jobid =~ s/Submitted batch job //; $jobstr .= "$jobid ";
      print " $cntjob submit $jobid\t$jobfile\n"; sleep $SLEEPTIME;
    } else { print " $cntjob skipped\t$jobfile\n"; }
    print LG $exec;

    # create and submit job: zcat to count fastq.gz for raw read numbers 
    $exec = " $SCRIPT_HOME/scripts/zcatdirfastq.pl -id $indir/$dir_list -o $tmpdir/$dir_list.zcat\n";
    $jobname = "pccZ.$dir_list"; $jobfile = "$logdir/$jobname.$T.job"; $cntjob++;
    $JOBHEAD = "\#!\/bin\/sh\n\n\#SBATCH -N 1\n\#SBATCH -J $jobname\n\#SBATCH -e $tmpdir/$jobname.$T.err\n\#SBATCH -o $tmpdir/$jobname.$T.out\n\n";
    open ( JB, "> $jobfile" ) or print "\nFail to open $jobfile\n" and die; print JB "$JOBHEAD\n$exec"; close ( JB );
    if ( $runzcat eq "yes" && $runjob eq "yes" ){ 
      do { $jobnum = &checkUserJob(); if ( $jobnum > $JOBLIMIT ){ sleep $JOBSLEEPTIME; print "$jobnum "; }} while ( $jobnum > $JOBLIMIT );
      $jobid = `$SBATCH $jobfile`; chomp( $jobid ); $jobid =~ s/Submitted batch job //; $jobstr .= "$jobid ";
      print " $cntjob submit $jobid\t$jobfile\n"; sleep $SLEEPTIME;
    } else { print " $cntjob skipped\t$jobfile\n"; }
    print LG $exec;

    # create and submit job: call picard tool for duplication 
    $exec = " $SCRIPT_HOME/scripts/dup_via_picard.pl -pp $PICARD -id $indir/$dir_list -o $tmpdir/$dir_list.picard.matrics\n";
    $jobname = "pccP.$dir_list"; $jobfile = "$logdir/$jobname.$T.job"; $cntjob++;
    $JOBHEAD = "\#!\/bin\/sh\n\n\#SBATCH -N 1\n\#SBATCH -J $jobname\n\#SBATCH -e $tmpdir/$jobname.$T.err\n\#SBATCH -o $tmpdir/$jobname.$T.out\n\n";
    open ( JB, "> $jobfile" ) or print "\nFail to open $jobfile\n" and die; print JB "$JOBHEAD\n$exec"; close ( JB );
    if ( $rundup eq "yes" && $runjob eq "yes" ){ 
      do { $jobnum = &checkUserJob(); if ( $jobnum > $JOBLIMIT ){ sleep $JOBSLEEPTIME; print "$jobnum "; }} while ( $jobnum > $JOBLIMIT );
      $jobid = `$SBATCH $jobfile`; chomp( $jobid ); $jobid =~ s/Submitted batch job //; $jobstr .= "$jobid ";
      print " $cntjob submit $jobid\t$jobfile\n"; sleep 3;
      print "jobnum=$jobnum\n";
    } else { print " $cntjob skipped\t$jobfile\n"; }
    print LG $exec;
  }
}
closedir DH;

if ( $runjob eq "yes" || $runjob eq "Yes" ){
  do { $jobnum = &checkKnownJob( $jobstr ); if ( $jobnum > 0 ){ sleep $JOBSLEEPTIME; print "$jobnum "; }} while ( $jobnum > 0 );
}

# --------------------------------------------------------------------
# pull all data together to create output summary file                
# --------------------------------------------------------------------
$exec = " $SCRIPT_HOME/scripts/sum_cvgFbamAnnotCont.pl -id $indir -c $configfile -od $outdir -td $tmpdir -ot $output -l $label -t $technology";
if ( $samplesheet ne "" ){ $exec .= " -ss $samplesheet "; } $exec .= "\n";
system( $exec ); print LG $exec; print $exec; close ( LG );
print "\n\n *** $Prg finished\n\n";

# ----------------------------------------------------------------------------------------------------
# sub functions
# ----------------------------------------------------------------------------------------------------
sub lastName {(my $sub_input) = @_; my @sub_inline = split(/\//, $sub_input); my $sub_fname = $sub_inline[scalar(@sub_inline)-1]; return $sub_fname;}

sub gettimestamp {
  my $syear  = (localtime(time-86400))[5]+1900;
  my $smonth = (localtime(time-86400))[4]+1; if ( $smonth < 10 ) { $smonth = "0" . $smonth; }
  my $sday   = (localtime(time-86400))[3];   if ( $sday   < 10 ) { $sday   = "0" . $sday;   }
  my $shour  = (localtime(time-86400))[2];   if ( $shour  < 10 ) { $shour  = "0" . $shour;  }
  my $smin   = (localtime(time-86400))[1];   if ( $smin   < 10 ) { $smin   = "0" . $smin;   }
  my $ssec   = (localtime(time-86400))[0];   if ( $ssec   < 10 ) { $ssec   = "0" . $ssec;   }
  my $sstamp = $syear . $smonth . $sday . $shour . $smin . $ssec;
  return $sstamp;
}

sub checkUserJob{
  my $sub_w = `whoami`; chomp ( $sub_w );
  my $sub_jobline = `$SQUEUE -u $sub_w`;

  my @sub_tp = split ( /\n/, $sub_jobline );
  my $sub_job = scalar(@sub_tp);
  return $sub_job;
}

sub checkKnownJob {
  my $sub_jobnum = 0;
  my ( $sub_jobstr ) = @_;
  my @sub_jobstrin = split ( /\s+/, $sub_jobstr );
  my $sub_jobline = `$SQUEUE`;
  my @sub_line = split ( /\n/, $sub_jobline );
  for ( my $sub_n = 0; $sub_n < scalar(@sub_line); $sub_n++ ){
    my @sub_tp = split ( /\s+/, $sub_line[$sub_n] );
    for ( my $sub_c = 0; $sub_c < scalar(@sub_tp); $sub_c++ ){
      for ( my $sub_i = 0; $sub_i < scalar(@sub_jobstrin); $sub_i++ ){ 
        if ( $sub_tp[$sub_c] eq $sub_jobstrin[$sub_i] ){ $sub_jobnum++; }
      }
    }
  }
  return $sub_jobnum;
}

