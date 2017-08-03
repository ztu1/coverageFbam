#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset
# set -o xtrace
unset module

usage() {
cat << EOF
usage: $0 options

This script performs QC analysis on a given run placing its results in a <run>/QC/ subdirectory.

OPTIONS:
  -h    help, Show this message
  -r    run, is the path to the run that we want to find
  -c    config, the configuration to use
  -d    debug, this will output more verbose messages for debugging purposes


Examples
  $0 -r /path/to/run -c /path/to/configdir

EOF
}

while getopts "hdr:c:" OPTION
do
  case $OPTION in
    h) usage ; exit 1 ;;
    d) DEBUG="-d" ; set -x ;;
    r) RUN_DIR=`readlink -f "$OPTARG"` ;;
    c) QC_CONFIG_DIR="${OPTARG}" ;;
    ?) usage ; exit 0;;
  esac
done


if [ -z "${RUN_DIR+x}" ] ; then
  echo "The run directory was not provied: RUN_DIR"
  exit 1
fi

if [ ! -d "${RUN_DIR}" ] ; then
  echo "The provided run directory does not exist: $RUN_DIR"
  exit 2
fi

if [ ! -d "${RUN_DIR}/samples" ] ; then
  echo "The provided run directory does not have a samples dir and might not be valid: $RUN_DIR"
  exit 21
fi

if [ -z "${QC_CONFIG_DIR+x}" ] ; then
  echo "The -c config directory was not provied: QC_CONFIG_DIR"
  exit 3
fi

if [ ! -d "${QC_CONFIG_DIR}" ] ; then
  echo "The given configuration directory does not exist: ${QC_CONFIG_DIR}"
  exit 4
fi

source ${QC_CONFIG_DIR}/tool.info


#make sure output area exists
if [ ! -d "${QC_OUT_DIR}" ] ; then
    echo "The QC out dir does not exist: ${QC_OUT_DIR}"
    exit 5
fi

QC_RUN_NAME=$(basename ${RUN_DIR})
#TODO better way to determine ordered service
QC_OUT_DIR="${RUN_DIR}/QC"

QC_OS=`echo ${QC_RUN_NAME} | cut -d_ -f1`
QC_OUTPUT="${QC_OUT_DIR}/${QC_RUN_NAME}.cvgdup"
QC_TEMP_DIR="${QC_OUT_DIR}/logs"
QC_LOG_FILE="${QC_TEMP_DIR}/runCoverage.log"


#data is valid, create a log area
if [ -d "${QC_OUT_DIR}" ] ; then
    echo "The QC directory already exists, indicating this has already been run: ${QC_OUT_DIR}"
    exit 6
fi

#make sure test defs exist
if [ ! -d "${QC_TD_DIR}/${QC_OS}" ] ; then
    echo "The ordered service in the config dir does not exist: ${QC_TD_DIR}/${QC_OS}"
    exit 7 
fi

#make sure test defs exist
if [ ! -f "${QC_TD_DIR}/${QC_OS}/QC.cfg" ] ; then
    echo "The QC.cfg in the test def does not exist: ${QC_TD_DIR}/${QC_OS}/QC.cfg"
    exit 8
fi


mkdir -p ${QC_TEMP_DIR}
echo "Beggining runCoverage execution." >> "${QC_LOG_FILE}"

#log the config used
cp "${QC_CONFIG_DIR}/tool.info" "${QC_TEMP_DIR}/tool.info"
cp "${QC_TD_DIR}/${QC_OS}/QC.cfg" "${QC_TEMP_DIR}/QC.cfg"

SAMPLE_SHEET=
if [[ ! -n $(find ${RUN_DIR} -maxdepth 1 -name "*SampleSheet.csv" -o -name "*SampleSheet.txt") ]] ; then
    echo "Sample not found in ${RUN_DIR}"
    exit 20
else
    SAMPLE_SHEET=$(find ${RUN_DIR} -maxdepth 1 -name "*SampleSheet.csv")
fi

echo "Sample Sheet used: ${SAMPLE_SHEET}" >> "${QC_LOG_FILE}"


QC_CMD="${QC_HOME}/run_cvgFbamAnnotCont_on_sampledir.pl -id ${RUN_DIR}/samples -ot ${QC_RUN_NAME} -od ${QC_OUTPUT} -td ${QC_TEMP_DIR} -t ${TECHNOLOGY} -q ${QC_QUEUE} -c ${QC_TD_DIR}/${QC_OS}/QC.cfg -rj yes"
echo "cmd: ${QC_CMD}" >> "${QC_LOG_FILE}"

#TODO run the command in a qsub
echo "QC coverage running, see logs for details: ${QC_TEMP_DIR}"
`${QC_PERL} ${QC_CMD} >> ${QC_LOG_FILE}`

