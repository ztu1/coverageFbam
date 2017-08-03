USED FOR: this package can be used for generate coverage report for sample under specific strucure (see INPUT DATA STRUCTURE)

---------------------------------------------------------------------------------------------------------------
USGAGE:
---------------------------------------------------------------------------------------------------------------

   run_cvgFbamAnnotCont_on_sampledir.pl
                 -id indir          # REQUIRED input top directory
                  -c configfile     # REQUIRED panel specific configure file
                 -od outdir         # REQUIRED outdir name
                 -ot output         # REQUIRED output name suggest runbatch name as output name
                [-ss samplesheet    # OPTION Illumina SampleSheet.csv or Ion Torrent SampleSheet.txt file
                  -t technology     # CHOICE: illumina, iontorrent; DEFAULT illumina
                  -l label          # Torrent sample name (I006-12345678) has label CHOICE: yes no, DEFAULT yes
                 -td tmpdir         # tmpdir for job files default outdir
                 -rb runbedtools    # run bedtools DEFAULT yes
                 -rd runduplication # run deplication DEFAULT yes
                 -rz runzcat        # run zcat for read num DEFAULT yes
                  -q queue          # job queue DEFAULT sandbox.q
                 -rj runjob         # run job, DEFAULT no
                  -l logfile        # logfile
                 -dg debug          # for development use DEFAULT no]

---------------------------------------------------------------------------------------------------------------
INPUT DATA STRUCTURE:
---------------------------------------------------------------------------------------------------------------

   samples/<sampleA>
          /<sampleB>

  under each sampleX, it will looks for
    <sampleX>/<sampleX>_xxx.fastq.gz
             /<sampleX>.bam

---------------------------------------------------------------------------------------------------------------
MAIN SUB STEPS:
---------------------------------------------------------------------------------------------------------------

 It will them call a few sub-program under <scripts> to do 
   call scripts/dlmp_zcatdirfastq.pl to count fastq.gz files for total raw reads number
   call scripts/dlmp_bam2cvg.pl with input of bed and bam files and call bedtools to generate coverage files
   call scripts/dlmp_dup_via_picard.pl with picard tools to generate duplicate data
   call dupNum_mq0_bam.pl to get duplicate rate on different percentile
   
   each above task will be a job submitted to the queue.
   when all jobs finished, it will call dlmp_sum_cvgFbamAnnotCont.pl to generate summary output file

---------------------------------------------------------------------------------------------------------------
CONFIG FILE:
---------------------------------------------------------------------------------------------------------------

  A config file is needed. check config/template.cfg as example.
  Below fields inside config file are needed for this package

    cvgFbamPipeRoot=/dlmp/sandbox/cgslIS/codes/coverageFbam/scripts
    genelist=ASXL1 BCOR BRAF CALR CBL
    bedfile1=/dlmp/sandbox/cgslIS/proj/DavidViswanatha/ref_tracks/NGS25-NGSHMv0.4.bed
    bedfile2=/dlmp/sandbox/cgslIS/ref_tracks/hg19_93phix.bed
    BEDTOOLS=/usr/local/biotools/bedtools/2.20.1/bin/bedtools
    lowcoveragelimitIllumina=500
    lowcoveragelimitTorrent=100
  
NOTE: the setting (queue) is only tested on cgsl linux cluster

---------------------------------------------------------------------------------------------------------------
RUN EXAMPLE:
---------------------------------------------------------------------------------------------------------------

 run command example:

  /dlmp/sandbox/cgslIS/codes/coverageFbam/run_cvgFbamAnnotCont_on_sampledir.pl \
  -id /dlmp/prod/runs/NGSHM/NGSHM_20160321_SSXT133_MS258B_AM8HB/samples \
  -od /dlmp/sandbox/runs/prodQC/NGSHM/NGSHM_20160321_SSXT133_MS258B_AM8HB/NGSHM_20160321_SSXT133_MS258B_AM8HB.cvgdup \
  -td /dlmp/sandbox/cgslIS/proj/DavidViswanatha/scratch/NGSHM_20160321_SSXT133_MS258B_AM8HB.20160324124123.log \
  -c /dlmp/sandbox/cgslIS/proj/DavidViswanatha/bin/NGSHM.cfg \
  -t illumina -ot NGSHM_20160321_SSXT133_MS258B_AM8HB -q sandbox.q -rj yes

---------------------------------------------------------------------------------------------------------------
Output file example
---------------------------------------------------------------------------------------------------------------
when process finished, there will be a tab delimited output file generated such as:
 /dlmp/sandbox/runs/prodQC/NGSHM/NGSHM_20160321_SSXT133_MS258B_AM8HB/NGSHM_20160321_SSXT133_MS258B_AM8HB.cvgdup/NGSHM_20160321_SSXT133_MS258B_AM8HB.coverage_summary.xls

SAMPLE-BASED MAPPING, VARIANCE COUNT, COVERAGE, AND PHIX CONTROL SUMMARY TABLE
Samples ReadsNum        ReadsMap2Ref    Pct2ReadsNum    ReadsMap2Target Pct2ReadsNum    totSNP  SNP2Target      totDIP  DIP2Target      MaxCVG  AvgCVG  MinCVG  TargetBasesWOphix       Bases>=x500     Pct     Bases<x500      Pct     PhixNum Pct2ReadNum     DupRate DupNum25%perMreads      DupNum50%perMreads      DupNum75%perMreads      DupNum100%perMreads     MapQuality0Pct  TsNum   TvNum   TsTvRatio       TsNumOnTarget   TvNumOnTarget   TsTvRatioOnTarget
10067670251     6014032 5995172 99.7    1714533 28.5    19162   23      5773    2       5800    3467.8  870     47763   47763   100.0   0       0.0     48      0.00    0.47    1       1       2       11      2.094   11156   8006    1.4     19      4       4.8
10067639750     5088940 5073420 99.7    1415272 27.8    18013   13      5700    2       4261    2808.0  704     47763   47763   100.0   0       0.0     34      0.00    0.44    1       1       2       11      2.211   10233   7780    1.3     10      3       3.3
10067655185     5563820 5546643 99.7    1787303 32.1    17449   14      5435    2       5017    3605.5  821     47763   47763   100.0   0       0.0     32      0.00    0.47    1       1       2       14      2.116   10168   7281    1.4     9       5       1.8
10067447564     5194580 5175258 99.6    1537131 29.6    16850   17      5352    2       4693    3088.0  600     47763   47763   100.0   0       0.0     43      0.00    0.46    1       1       2       14      2.125   9977    6873    1.5     14      3       4.7
10067561880     6004814 5987416 99.7    1798338 29.9    19757   21      6358    3       5711    3597.1  956     47763   47763   100.0   0       0.0     24      0.00    0.48    1       1       2       12      2.085   11457   8300    1.4     17      4       4.2
RunBatchTotal   27866186                                                                                3313            238815  238815  100.0   0       0.0     181     0.0

LANE-BASED COVERAGE SUMMARY TABLE
LANE    lane1   lane2   lane3   lane4   lane5   lane6   lane7   lane8
CVG     3313    0       0       0       0       0       0       0


GENE-BASED COVERAGE TABLE
Gene    ASXL1   BCOR    BRAF    CALR    CBL     CBL-G   CEBPA   CSF3R   DNMT3A  ETV6    EZH2    EZH2-G  FLT3    GATA1   GATA2   IDH1    IDH2    JAK2    KIT     KRAS    MPL     MYD88   NOTCH1  NPM1    NRAS    PHF6    PTPN11  RUNX1   SETBP1  SF3B1   SRSF2   TERT    TET2    TP53    U2AF1   WT1     ZRSR2ASXL1      ZRSR2   PHIX
ExonNum 5       13      1       2       2       2       2       3       17      7       19      1       7       3       7       1       1       5       5       3       2       2       4       4       3       10      4       7       1       4       3       16      11      6       4       12              13
CVG     4076    2375    2793    2873    3877    3470    2567    3761    3556    3670    3144    3012    3430    2007    3410    3714    3397    2903    3725    3050    3484    3359    3445    2654    3785    1902    3521    3399    4248    3453    3793    3249    3848    3550    3180    3352    0       2148    0

SAMPLE AND GENEorAMPLICON-BASED SUMMARY TABLE
Sample  Gene    Exon    ExonLength      MaxCVG  AvgCVG  MinCVG  Note    ExonCvg2GeneCvg cntAvgLow=0
10067670251     EZH2-G  3u      6       2651    2641    2622            1.00
10067670251     KIT     8       134     5330    4828    3643            1.01
10067670251     KIT     9       213     5032    4305    2563            1.13
10067670251     KIT     10      126     4470    4157    3305            0.91
10067670251     KIT     11      146     4234    4082    3626            1.03
............
............
............

---------------------------------------------------------------------------------------------------------------
MORE INFO ABOUT OUTPUT FILE 
---------------------------------------------------------------------------------------------------------------
output file has multiple tables:

First Table: SAMPLE-BASED MAPPING, VARIANCE COUNT, COVERAGE, AND PHIX CONTROL SUMMARY TABLE
  it has multiple columns with information about each samples as well as batch total and phix 
    Samples		sample name
    ReadsNum		raw reads from fastq.gz
    ReadsMap2Ref	read number mapped to reference genome
    Pct2ReadsNum	percentage of ReadsMap2Ref divided by ReadsNum
    ReadsMap2Target	read number mapped to target (bedfile1)
    Pct2ReadsNum	percentage of ReadsMap2Target divided by ReadsNum
    totSNP		cmb_vcf SNP variant number
    SNP2Target		cmb_vcf SNP variant number on target
    totDIP		cmb_vcf INDEL variant number
    DIP2Target		cmb_vcf INDEL variant number on target
    MaxCVG		Max coverage
    AvgCVG		Average coverage
    MinCVG		Min coverage
    TargetBasesWOphix	target bed (bedfile1) bases without phix
    Bases>=x500		target region bases that has coverage equal or more than lowcoveragelimitIllumina or lowcoveragelimitTorrent
    Pct			percentage Bases>=x500 divided by target bed size 
    Bases<x500		target region bases that has coverage less lowcoveragelimitIllumina or lowcoveragelimitTorrent
    Pct			percentage Bases<x500 divided by target bed size
    PhixNum		Phix read number
    Pct2ReadNum		percentage PhixNum divided by ReadsNum
    DupRate		Duplication rate from picards tools
    DupNum25%perMreads	Duplication number when  25% total read counted normalized by 1 million reads
    DupNum50%perMreads	Duplication number when  50% total read counted normalized by 1 million reads
    DupNum75%perMreads	Duplication number when  75% total read counted normalized by 1 million reads
    DupNum100%perMreads	Duplication number when 100% total read counted normalized by 1 million reads
    MapQuality0Pct	Reads that have map quality=0 divided by ReadsNum

Second Table: LANE-BASED COVERAGE SUMMARY TABLE
   it has total coverage per lane (up to 8 lanes). It will not apply for ion torrent.
   it reflect lane variation (designed for GAIIx initially)

     row1	LANE   		laneX lane number    
     row2	CVG             average coverage for particular lane on the target (bedfile1)

Third Table: GENE-BASED COVERAGE TABLE
  It has compiled average coverage per gene. Used to see gene coverage efficiency
  
     row1	Gene		List of gene from genelist inside configure file
     row2	ExonNum         exon or amplicon number per gene (read from bedfile1) 
     row3	CVG     	compiled average coverage for particular gene   

Fourth Table: SAMPLE AND GENEorAMPLICON-BASED SUMMARY TABLE
  It has multiple columns and split into per sample, per gene, per exon/amplicon base
  Used for overall coverage performance 

    Sample  		Sample name
    Gene    		Gene name
    Exon    		Exon number 
    ExonLength      	Exon size
    MaxCVG  		Max Coverage
    AvgCVG  		Average Coverage
    MinCVG  		Min Coverage
    Note    		Not if any
    ExonCvg2GeneCvg 	Exon average coverage divided by Gene average coverage
    cntAvgLow=X		If any ExonCvg2Gene less than 0.01, it will mark as low on that exon. X is the total of marks
			Initial design to monitor exon drop.  It should be replaced by dup/del detection method.

---------------------------------------------------------------------------------------------------------------
NOTE
---------------------------------------------------------------------------------------------------------------
Most of functionality has been re-engined into RoQCM
