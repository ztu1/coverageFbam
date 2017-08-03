#!/bin/bash

read -r -d '' DOCS <<DOCS
This script will semi-automatically unit test the runCoverage wrapper script.
Both negative (failures) and postive successes (previous runs) will be unit tested.
Previously generated runs will be compared via check sums and diffs. 
This is because unit tests should be fast and not be held up by several minute cluster runs.
You will need to make sure this script knows about the following locations--
1. Where the test configuration data is located
2. Where the test input data is located  
3. Where the source files are located (due to ${QC_HOME}/run_cvgFbamAnnotCont_on_sampledir.pl where QC_HOME will normally point to /dlmp/dev/scripts/sources/coverage<xxx>)
DOCS

#comment out if using bashdb 
set -o pipefail 
set -o nounset 
unset module

#Unit Tests setUp defines instructions that will be executed before each test
UNITTEST_ANYVALIDRUNDIR=/dlmp/dev/runs/NGSHM/NGSHM_20160321_SSXT133_MS258B_AM8HB_keesh
UNITTEST2_RUNDIR=/dlmp/dev/runs_away  #fake
UNITTEST3_RUNDIR=/home/m008480/testing/coverageFbam3/runs  #fake
UNITTEST7_RUNDIR=/dlmp/dev/runs/NGSHM/NGSHM_20160321_SSXT133_MS258B_AM8HB_keesh_unittest7
UNITTEST8_RUNDIR=/dlmp/dev/runs/NGSHM/NGSHM_20160321_SSXT133_MS258B_AM8HB_keesh_unittest8
UNITTEST9_RUNDIR=/dlmp/dev/runs/HCCP/HCCP_CGSLV135_V2R2_AGWL9_keesh_unittest
UNITTEST10_RUNDIR=/dlmp/dev/runs/NGSHM/NGSHM_20160321_SSXT133_MS258B_AM8HB_keesh_unittest10 

#The following unit test CONFIG dirs are also used (and must contain a copy of the complete source code-- perhaps a symbolic link is needed or we should live under
#/dlmp/dev/scripts/sources/coverage<xxx>)
#/home/m008480/testing/coverageFbam3  
#The following dirs are referred to by /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/testing/unittest<x>/tool.info
#/home/m008480/testing/coverageFbam6  
#/home/m008480/testing/coverageFbam8  
#/home/m008480/testing/coverageFbam9  

MISSING_QC_CFG_DIR=${UNITTEST9_RUNDIR}/QC
if [ -d "${MISSING_QC_CFG_DIR}" ] ; then
     rmdir "${MISSING_QC_CFG_DIR}"
fi
SAMPLE_NOT_FOUND_IN_RUN_DIR=${UNITTEST10_RUNDIR}/QC
if [ -d "${SAMPLE_NOT_FOUND_IN_RUN_DIR}" ] ; then
     rm -rf "${SAMPLE_NOT_FOUND_IN_RUN_DIR}"
fi
#end setUp


#Unit tests (negative)
#"The run directory was not provied"	exit 1
#Unit Test 1
./runCoverage.sh -c illumina
if [ ! $? -eq 1 ]; then
    exit 255
fi

#"The provided run directory does not exist"  exit 2
#Unit Test 2
./runCoverage.sh -r ${UNITTEST2_RUNDIR} -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/illumina
if [ ! $? -eq 2 ]; then
    exit 255
fi

#"The provided run directory does not have a samples dir and might not be valid" exit 21
#Unit Test 3
./runCoverage.sh -r ${UNITTEST3_RUNDIR} -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/illumina
if [ ! $? -eq 21 ]; then
    exit 255
fi

#"The -c config directory was not provied" exit 3 
#Unit Test 4
./runCoverage.sh -r ${UNITTEST_ANYVALIDRUNDIR}
if [ ! $? -eq 3 ]; then
    exit 255
fi

#"The given configuration directory does not exist"  exit 4
#Unit Test 5
./runCoverage.sh -r ${UNITTEST_ANYVALIDRUNDIR} -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/illuminati
if [ ! $? -eq 4 ]; then
    exit 255
fi

#"The QC out dir does not exist", ${QC_OUT_DIR} exit 5 
#Unit Test 6
#ex.) QC_OUT_DIR="${STAGE_PATH}/runs/QC"
#needed to edit unittest/tool.info for this test only
#may need another unittest dir?
./runCoverage.sh -r ${UNITTEST_ANYVALIDRUNDIR} -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/testing/unittest6
if [ ! $? -eq 5 ]; then
    exit 255
fi

#"The QC directory already exists, indicating this has already been run" exit 6 
#Unit Test 7
#QC_OUT_DIR="${RUN_DIR}/QC"
./runCoverage.sh -r ${UNITTEST7_RUNDIR} -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/illumina
if [ ! $? -eq 6 ]; then
    exit 255
fi

#"The ordered service in the config dir does not exist"  exit 7 
#Unit Test 8
#STAGE_PATH="/home/m008480/testing/coverageFbam2"
#export QC_HOME="${STAGE_PATH}/scripts/sources/coverageFbam"
#QC_TD_DIR="${QC_HOME}/config/QC"
#/home/m008480/testing/coverageFbam2/scripts/sources/coverageFbam/config/QC
#f [ ! -d "${QC_TD_DIR}/${QC_OS}" ] ; then
./runCoverage.sh -r ${UNITTEST8_RUNDIR} -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/testing/unittest8
if [ ! $? -eq 7 ]; then
    exit 255
fi


#"The QC.cfg in the test def does not exist"  exit 8 
#Unit Test 9
#STAGE_PATH="/home/m008480/testing/coverageFbam4cd "  
#export QC_HOME="${STAGE_PATH}/scripts/sources/coverageFbam"
#QC_TD_DIR="${QC_HOME}/config/QC"
#"{QC_TD_DIR}/${QC_OS}/QC.cfg" 
./runCoverage.sh -r ${UNITTEST9_RUNDIR} -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/testing/unittest9
if [ ! $? -eq 8 ]; then
    exit 255
fi

#"Sample not found in"  exit 20
#Unit Test 10
#if [[ ! -n $(find ${RUN_DIR} -maxdepth 1 -name "*SampleSheet.csv" -o -name "*SampleSheet.txt") ]] ; then
./runCoverage.sh -r ${UNITTEST10_RUNDIR} -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/illumina
if [ ! $? -eq 20 ]; then
    exit 255
fi




#Positive tests
#./runCoverage.sh -r /dlmp/dev/runs/NGSHM/NGSHM_20160321_SSXT133_MS258B_AM8HB_keesh -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/illumina
#./runCoverage.sh -r /dlmp/dev/runs/HCCP/HCCP_CGSLV135_V2R2_AGWL3_keesh -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/illumina

#./runCoverage.sh -r /dlmp/dev/runs/GISTP/GISTP_032216-mgt56-3_TSN131_MS259_AM8NE_keesh -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/illumina
  #Not run from Tu's README Word document
  #./runCoverage.sh -r /dlmp/dev/runs/NGSHM/NGSHM_20150928_SSXT057R_MS133D_AJRGR_keesh -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/illumina

#./runCoverage.sh -r /dlmp/dev/runs/NONTP/NONTP_V08p_Chip1_LIBDIT002_08P_Pet-28_keesh -c /dlmp/dev/scripts/sources/coverageFbam/config/wrapper/iontorrent


#
#case 1.  NGSHM_20160321_SSXT133_MS258B_AM8HB
RECIPE="NGSHM"
FLAVOR="NGSHM_20160321_SSXT133_MS258B_AM8HB" 
FOLDER_GOLD="/dlmp/dev/runs/GOLDEN/QC/${FLAVOR}/QC/${FLAVOR}.cvgdup"
FOLDER="/dlmp/dev/runs/${RECIPE}/${FLAVOR}_keesh/QC/${FLAVOR}_keesh.cvgdup"
SPREADSHEET_GOLD="${FOLDER_GOLD}/${FLAVOR}.coverage_summary.xls"
SPREADSHEET="${FOLDER}/${FLAVOR}_keesh.coverage_summary.xls"
FOLDER_GOLD="${FOLDER_GOLD}/out"
FOLDER="${FOLDER}/out"
#
CHKSUM1=`sha1sum ${SPREADSHEET_GOLD} | cut -d' ' -f1`
CHKSUM2=`sha1sum ${SPREADSHEET} | cut -d' ' -f1`
echo "${CHKSUM1} : ${CHKSUM2}"
if [ ${CHKSUM1} != ${CHKSUM2} ]; then
    echo "${SPREADSHEET_GOLD} not equal to ${SPREADSHEET} !"
    exit 255
fi
DIFF1=`diff -q -x '*.mq?' -x '*.matrics' ${FOLDER_GOLD} ${FOLDER}`
echo ${DIFF1}
if [ -n "${DIFF1}" ]; then
    echo "${FOLDER_GOLD} and ${FOLDER} are different!"
    exit 255
fi 

#case 2. HCCP_CGSLV135_V2R2_AGWL3
RECIPE="HCCP"
FLAVOR="HCCP_CGSLV135_V2R2_AGWL3"
FOLDER_GOLD="/dlmp/dev/runs/GOLDEN/QC/${FLAVOR}/QC/${FLAVOR}.cvgdup"
FOLDER="/dlmp/dev/runs/${RECIPE}/${FLAVOR}_keesh/QC/${FLAVOR}_keesh.cvgdup"
SPREADSHEET_GOLD="${FOLDER_GOLD}/${FLAVOR}.coverage_summary.xls"
SPREADSHEET="${FOLDER}/${FLAVOR}_keesh.coverage_summary.xls"
FOLDER_GOLD="${FOLDER_GOLD}/out"
FOLDER="${FOLDER}/out"
#
CHKSUM1=`sha1sum ${SPREADSHEET_GOLD} | cut -d' ' -f1`
CHKSUM2=`sha1sum ${SPREADSHEET} | cut -d' ' -f1`
echo "${CHKSUM1} : ${CHKSUM2}"
if [ ${CHKSUM1} != ${CHKSUM2} ]; then
    echo "${SPREADSHEET_GOLD} not equal to ${SPREADSHEET} !"
    exit 255
fi
DIFF1=`diff -q -x '*.mq?' -x '*.matrics' ${FOLDER_GOLD} ${FOLDER}`
if [ -n "${DIFF1}" ]; then
    echo "${FOLDER_GOLD} and ${FOLDER} are different!"
    exit 255
fi 

#case #3 GISTP_032216-mgt56-3_TSN131_MS259_AM8NE was NGSHM_20150928_SSXT057R_MS133D_AJRGR
RECIPE="GISTP" # was "NGSHM"
FLAVOR="GISTP_032216-mgt56-3_TSN131_MS259_AM8NE"  #was NGSHM_20150928_SSXT057R_MS133D_AJRGR" which was not run from Tu's case #3 from his README Word doc
FOLDER_GOLD="/dlmp/dev/runs/GOLDEN/QC/${FLAVOR}/QC/${FLAVOR}.cvgdup"
FOLDER="/dlmp/dev/runs/${RECIPE}/${FLAVOR}_keesh/QC/${FLAVOR}_keesh.cvgdup"
SPREADSHEET_GOLD="${FOLDER_GOLD}/${FLAVOR}.coverage_summary.xls"
SPREADSHEET="${FOLDER}/${FLAVOR}_keesh.coverage_summary.xls"
FOLDER_GOLD="${FOLDER_GOLD}/out"
FOLDER="${FOLDER}/out"
#
CHKSUM1=`sha1sum ${SPREADSHEET_GOLD} | cut -d' ' -f1`
CHKSUM2=`sha1sum ${SPREADSHEET} | cut -d' ' -f1`
echo "${CHKSUM1} :  ${CHKSUM2}"
if [ ${CHKSUM1} != ${CHKSUM2} ]; then
    echo "${SPREADSHEET_GOLD} not equal to ${SPREADSHEET} !"
    exit 255
fi
DIFF1=`diff -q -x '*.mq?' -x '*.matrics' ${FOLDER_GOLD} ${FOLDER}`
if [ -n "${DIFF1}" ]; then
    echo "${FOLDER_GOLD} and ${FOLDER} are different!"
    exit 255
fi 

#case #4 NONTP_V08p_Chip1_LIBDIT002_08P_Pet-28
RECIPE="NONTP"
FLAVOR="NONTP_V08p_Chip1_LIBDIT002_08P_Pet-28"
FOLDER_GOLD="/dlmp/dev/runs/GOLDEN/QC/${FLAVOR}/QC/${FLAVOR}.cvgdup"
FOLDER="/dlmp/dev/runs/${RECIPE}/${FLAVOR}_keesh/QC/${FLAVOR}_keesh.cvgdup"
SPREADSHEET_GOLD="${FOLDER_GOLD}/${FLAVOR}.coverage_summary.xls"
SPREADSHEET="${FOLDER}/${FLAVOR}_keesh.coverage_summary.xls"
FOLDER_GOLD="${FOLDER_GOLD}/out"
FOLDER="${FOLDER}/out"
#
CHKSUM1=`sha1sum ${SPREADSHEET_GOLD} | cut -d' ' -f1`
CHKSUM2=`sha1sum ${SPREADSHEET} | cut -d' ' -f1`
echo "${CHKSUM1} :  ${CHKSUM2}"
if [ ${CHKSUM1} != ${CHKSUM2} ]; then
    echo "${SPREADSHEET_GOLD} not equal to ${SPREADSHEET} !"
    exit 255
fi
DIFF1=`diff -q -x '*.mq?' -x '*.matrics' ${FOLDER_GOLD} ${FOLDER}`
if [ -n "${DIFF1}" ]; then
    echo "${FOLDER_GOLD} and ${FOLDER} are different!"
    exit 255
fi 

#tearDown allows you to define instructions that will be executed after each test

echo "${0##*/} : ALL unit tests passed!"

#EOT