#!/bin/sh
#
#  alleleload.sh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper around the entire Allele load process.
#
#  Usage:
#
#      alleleload.sh
#
#  Env Vars:
#
#      See the configuration file (alleleload.config)
#
#  Inputs:  None
#
#  Outputs:
#
#      - Log file (${LOG_DIAG})
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Fatal error occurred
#
#  Assumes:  Nothing
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Source the configuration file to establish the environment.
#      2) Establish the log file.
#      3) Call makeNomenFile.sh to generate the marker file from the input file.
#      4) Call loadNomen.sh to load the marker file into the database.
#      5) Call makeAlleleFile.sh to generate the allele files from the input file.
#      6) Call loadAllele.sh to load the allele files into the database.
#
#  Notes:  None
#
###########################################################################

cd `dirname $0`

CONFIG=alleleload.config

#
# Make sure the configuration file exists and source it.
#
if [ -f ../${CONFIG} ]
then
    . ../${CONFIG}
else
    echo "Missing configuration file: ${CONFIG}"
    exit 1
fi

NOMENCONFIG=$1

#
# Make sure the nomen configuration file exists and source it.
#
if [ -f ../${NOMENCONFIG} ]
then
    . ../${NOMENCONFIG}
else
    echo "Missing configuration file: ${NOMENCONFIG}"
    exit 1
fi

#
#  Source the DLA library functions.
#

if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" | tee -a ${LOG}
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." | tee -a ${LOG}
    exit 1
fi

#
# createArchive
#
preload ${OUTPUTDIR}

#
# Establish the log file.
#
LOG=${LOG_DIAG}
rm -rf ${LOG}
touch ${LOG}

#
# Create the Marker input file
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Call makeNomenFile.csh (nomenload/bin/)" | tee -a ${LOG}
${NOMENLOAD}/bin/makeNomenFile.sh ${NOMENCONFIG} 2>&1 >> ${LOG}
STAT=$?
checkStatus ${STAT} "makeNomenFile.sh (nomenload/bin/)"

#
# Load the Marker input file
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Run nomenload.py (nomenload/bin/)" | tee -a ${LOG}
${NOMENLOAD}/bin/nomenload.py 2>&1 >> ${LOG}
STAT=$?
checkStatus ${STAT} "loadNomen.sh (nomenload/bin/)"

#
# Create the Allele files.
#
#echo "" >> ${LOG}
#date >> ${LOG}
#echo "Call makeAlleleFile.sh (alleleload.sh)" | tee -a ${LOG}
#./makeAlleleFile.sh 2>&1 >> ${LOG}
#STAT=$?
#checkStatus ${STAT} "makeAlleleFile.sh (alleleload.sh)"

#
# Load Allele files
#
#echo "" >> ${LOG}
##date >> ${LOG}
#echo "Call loadAllele.sh (alleleload.sh)" | tee -a ${LOG}
#./loadAllele.sh ${JOBKEY} 2>&1 >> ${LOG}
#STAT=$?
#checkStatus ${STAT} "loadAllele.sh (alleleload.sh)"

#
# run postload cleanup and email logs
#
shutDown
exit 0
