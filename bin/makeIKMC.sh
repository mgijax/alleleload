#!/bin/sh
#
#  makeIKMC.sh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper around the process that creates the IKMC Allele file
#	and loads the Alleles into MGI
#
#  Usage:
#
#      makeIKMC.sh
#
#  Env Vars:
#
#      See the configuration file (ikmc.config)
#
#  Inputs:  ikmc.config
#
#  Outputs:
#
#      - Log files (${LOG_DIAG}, ${LOG_)
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
#      4) Call makeIKMC.py to create the allele file.
#      5) Call makeAllele.csh to create/load the allele file.
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
# Establish the log file.
#
LOG=${LOG_DIAG}
rm -rf ${LOG}
touch ${LOG}

#
# createArchive
#
echo "archiving..." >> ${LOG}
date >> ${LOG}
preload ${OUTPUTDIR}
rm -rf ${OUTPUTDIR}/*.diagnostics
rm -rf ${OUTPUTDIR}/*.error
echo "archiving complete" >> ${LOG}
date >> ${LOG}

#
# find & copy download file to input directory
#
if [ ! -r ${IKMC_INPUT_FILE} ]
then
    echo "Error: IKMC input file could not be found in download directory." | tee -a ${LOG}
    exit 1
fi

cp ${IKMC_INPUT_FILE} ${IKMC_COPY_INPUT_FILE} | tee -a ${LOG}

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

#
# copy ${OUTPUTDIR}/mgi_allele_ikmc.txt.new to ${IKMC_FTP} directory
# add a 'ln -s' to a static ftp file
#
#useDate=`date '+%m%d%y'`
#echo ${OUTPUTDIR}/mgi_allele_ikmc.txt.new ${IKMC_FTP}/mgi_allele_ikmc.txt.${useDate} | tee -a ${LOG}
#cp ${OUTPUTDIR}/mgi_allele_ikmc.txt.new ${IKMC_FTP}/mgi_allele_ikmc.txt.${useDate} | tee -a ${LOG}
STAT=$?
if [ ${STAT} -ne 0 ]
then
    echo "Error: problem copying output file to ftp site" | tee -a ${LOG}
    exit 1
fi

#
# run postload cleanup and email logs
#
shutDown
exit 0

