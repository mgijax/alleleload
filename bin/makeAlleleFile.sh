#!/bin/sh
#
#  makeAlleleFiles.sh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper around the process that generates
#      the Allele bcp files.
#
#  Usage:
#
#      makeAlleleFiles.sh
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
#      2) Verify that the input files exist.
#      3) Establish the log file.
#      4) Call makeAlleleFiles.py to generate the Allele bcp files.
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

#
# Establish the log file.
#
LOG=${LOG_DIAG}

#
# Call the Python script to generate the Allele bcp files.
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Make the Allele bcp files (makeAlleleFiles.sh)" | tee -a ${LOG}
./makeAlleleFiles.py 2>&1 >> ${LOG}
STAT=$?
if [ ${STAT} -ne 0 ]
then
    echo "Error: Make Allele files (makeAlleleFiles.sh)" | tee -a ${LOG}
    exit 1
fi

exit 0
