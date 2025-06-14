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
#      See the configuration file
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
rm -rf ${LOG_CUR}
touch ${LOG}
touch ${LOG_CUR}

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
    echo "Error: IKMC input file did not get created properly." | tee -a ${LOG}
    exit 1
fi

#
# Create the IKMC/Allele input file
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Create the IKMC/Allele input file (makeIKMC.sh)" | tee -a ${LOG}
${PYTHON} -W "ignore" ./makeIKMC.py 2>&1 >> ${LOG}
STAT=$?
checkStatus ${STAT} 'Create the IKMC/Allele input file (makeIKMC.sh)'

#
# Create the Alleles
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Create the Alleles (makeAllele.sh)" | tee -a ${LOG}
./makeAllele.sh ${CONFIG} 2>&1 >> ${LOG}
STAT=$?
checkStatus ${STAT} 'Create the Alleles (makeAllele.sh)'

# wts2-1030;11/07/2022;per Cindy, no longer needed at this time
# copy ${OUTPUTDIR}/mgi_allele_ikmc.txt.new to ${IKMC_FTP} directory
#
#echo "" >> ${LOG}
#date >> ${LOG}
#echo "Copying IKMC output file to ftp directory..." | tee -a ${LOG}
#useDate=`date '+%m%d%Y'`
#cd ${IKMC_FTP}
#rm -rf mgi_allele_ikmc.txt.${useDate} | tee -a ${LOG}
#rm -rf mgi_allele_ikmc.txt.current | tee -a ${LOG}
#cp ${OUTPUTDIR}/mgi_allele_ikmc.txt.new mgi_allele_ikmc.txt.${useDate} | tee -a ${LOG}
#ln -s mgi_allele_ikmc.txt.${useDate} mgi_allele_ikmc.txt.current | tee -a ${LOG}
#STAT=0
#checkStatus ${STAT} 'Copying IKMC output file to ftp directory'

#
# curator log
#
#wc -l ${INPUTDIR}/mgi_modification_allele_report.tsv | tee -a ${LOG_CUR}
#wc -l ${INPUTDIR}/mgi_allele_ikmc.txt | tee -a ${LOG_CUR}
#wc -l ${LOGDIR}/ikmc.exist.log | tee -a ${LOG_CUR}
#wc -l ${LOGDIR}/ikmc.skip.log | tee -a ${LOG_CUR}
#wc -l ${OUTPUTDIR}/* | tee -a ${LOG_CUR}

#
# run postload cleanup and email logs
#
shutDown
exit 0

