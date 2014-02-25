#!/bin/sh
#
#  makeTest.sh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper around the process that creates the Test Allele file
#
#  Usage:
#
#      makeTest.sh
#
#  Env Vars:
#
#      See the configuration file (sanger.config)
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
#      2) Verify that the input file exists.
#      3) Establish the log file.
#      4) Call makeTest.py to create the association file.
#
#  Notes:  None
#
###########################################################################

cd `dirname $0`

CONFIG=${1}

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

IKMC_INPUT_FILE=${ALLELELOAD}/test/DAL_pass_fail.txt
IKMC_COPY_INPUT_FILE=${INPUTDIR}/DAL_pass_fail.txt
INPUTFILE=${INPUTDIR}/DAL_pass_fail_mgi_allele_test.txt
LOG_PROC=${LOGDIR}/DAL_pass_fail.proc.log
LOG_DIAG=${LOGDIR}/DAL_pass_fail.diag.log
LOG_CUR=${LOGDIR}/DAL_pass_fail.cur.log
LOG_VAL=${LOGDIR}/DAL_pass_fail.val.log
LOG_IKMC=${LOGDIR}/DAL_pass_fail.log
SKIP_DIAG=${LOGDIR}/DAL_pass_fail.skip.log
EXISTS_DIAG=${LOGDIR}/DAL_pass_fail.exist.log

export IKMC_INPUT_FILE IKMC_COPY_INPUT_FILE INPUTFILE
export LOG_PROC LOG_DIAG LOG_CUR LOG_VAL LOG_IKMC SKIP_DIAG EXISTS_DIAG

#
# Establish the log file.
#
LOG=${LOG_DIAG}
rm -rf ${LOG}
touch ${LOG}

echo "" >> ${LOG}
date >> ${LOG}
echo "Create the Test/Allele input file (makeTest.sh)" | tee -a ${LOG}
${ALLELELOAD}/bin/makeIKMC.py
#${ALLELELOAD}/bin/makeIKMC.py 2>&1 >> ${LOG}
STAT=$?
if [ ${STAT} -ne 0 ]
then
    echo "Error: Create the Test/Allele input file (makeTest.sh)" | tee -a ${LOG}
    exit 1
fi

exit 0
