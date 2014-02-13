#!/bin/sh
#
#  makeIKMC.sh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper around the process that creates the IKMC Allele file
#
#  Usage:
#
#      makeIKMC.sh
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
#      4) Call makeIKMC.py to create the association file.
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

#
# Establish the log file.
#
LOG=${LOG_DIAG}
rm -rf ${LOG}
touch ${LOG}

#
# find & copy download file to input directory
#
if [ ! -r ${IKMC_INPUT_FILE} ]
then
    echo "Error: IKMC input file could not be found in download directory." | tee -a ${LOG}
    exit 1
fi

cp ${IKMC_INPUT_FILE} ${IKMC_COPY_INPUT_FILE}

if [ ! -r ${IKMC_COPY_INPUT_FILE} ]
then
    echo "Error: IKMC input file did not get created proporly." | tee -a ${LOG}
    exit 1
fi

#
# Create the IKMC/Allele input file
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Create the IKMC/Allele input file (makeIKMC.sh)" | tee -a ${LOG}
./makeIKMC.py 2>&1 >> ${LOG}
STAT=$?
if [ ${STAT} -ne 0 ]
then
    echo "Error: Create the IKMC/Allele input file (makeIKMC.sh)" | tee -a ${LOG}
    exit 1
fi
exit 0

#
# Create the Alleles
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Create the Alleles (makeAllele.sh)" | tee -a ${LOG}
./makeAllele.sh ${CONFIG} 2>&1 >> ${LOG}
STAT=$?
if [ ${STAT} -ne 0 ]
then
    echo "Error: Create the Alleles (makeAllele.sh)" | tee -a ${LOG}
    exit 1
fi

exit 0
