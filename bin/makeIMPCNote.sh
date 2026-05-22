#!/bin/sh

#
# wts2-1857/Allele IMPC Molecular Note
#
# Wrapper script for loading
#

cd `dirname $0`

CONFIG=${ALLELELOAD}/impcnote.config

#
# Make sure the configuration file exists and source it.
#
if [ -f ${CONFIG} ]
then
    . ${CONFIG}
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
# There should be a "lastrun" file in the input directory that was created
# the last time the gene model load was run. If this file exists and is more
# recent than the gene model file, the load does not need to be run.
#
LASTRUN_FILE=${INPUTDIR}/lastrun
if [ -f ${LASTRUN_FILE} ]
then
    if test ${LASTRUN_FILE} -nt ${INPUTFILE}
    then
        echo "Input files have not been updated - skipping load" | tee -a ${LOG}
        exit 0
    fi
fi

#
# copy file
#
rm -rf ${INPUTFILE}
cp ${IMPC_INPUT_FILE} ${INPUTFILE}

#
# process impc molecular note
#
date | tee -a ${LOG}
${PYTHON} makeIMPCNote.py | tee -a ${LOG}
STAT=$?
if [ ${STAT} -ne 0 ]
then
    echo "Error: impcnote.py" | tee -a ${LOG}
    exit 1
fi

touch ${LASTRUN_FILE}

exit 0
