#
# Program: impcmolnote.py
#
# Purpose:
#
#	To load the Allele IMPC Notes into _noteype_key = 1053
#
# Inputs:
#
#	A tab-delimited file in the format:
#
#   1) Gene Symbol
#   2) Allele Symbol
#   3) MGI Gene Id <-------
#   4) MGI Allele Id <-------
#   5) TPN
#   6) PIN
#   7) MIN
#   8) Produc>on Center
#   9) Endonuclease Type
#   10) Guide
#   11) Build
#   12) Background Strain
#   13) Ins Loca>on
#   14) Genome Browser Link
#   15) Dele>on Info
#   16) Descrip>on <------- this is the note field
#
# Sanity checks
#   a) The MGI Gene ID is not a valid MGI marker ID
#   b) The MGI Allele ID is not a valid MGI allele ID
#   c) The MGI Allele ID is not associated with the MGI Gene ID in the same row
#   d) Duplicate Lines in input file
#
# Outputs:
#
#       BCP files:
#
#       MGI_Note_IMPCMolecular.bcp
#
#       Diagnostics file of all input parameters and SQL commands
#       Error file
#
# 05/21/2026    sc
#       - wts2-1857/Allele Molecular Note Load
#

import sys
import os
import db
import mgi_utils
import loadlib

inputFileName = os.environ['INPUTFILE']
outputDir = os.environ['OUTPUTDIR']
BCP_COMMAND = os.environ['PG_DBUTILS'] + '/bin/bcpin.csh'

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor
inputFile = ''		# file descriptor
noteFile = ''		# file descriptor

noteTable = 'MGI_Note'
noteFileName = outputDir + '/' + noteTable + '.bcp'

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name

mgiTypeKey = 11
noteTypeKey = 1053
createdByKey = 1000

loaddate = loadlib.loaddate

#
# Purpose: prints error message and exits
#
def exit(
    status,          # numeric exit status (integer)
    message = None   # exit message (str.
    ):

    if message is not None:
        sys.stderr.write('\n' + str(message) + '\n')
 
    try:
        diagFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        errorFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        diagFile.close()
        errorFile.close()
        inputFile.close()
    except:
        print('issued closing files from exit function')
        pass

    sys.exit(status)
 
#
# Purpose: process command line options
#
def initialize():
    global diagFile, errorFile, inputFile, errorFileName, diagFileName
    global noteFile
    
    head, tail = os.path.split(inputFileName)

    diagFileName = outputDir + '/' + tail + '.diagnostics'
    errorFileName = outputDir + '/' + tail + '.error'

    try:
        diagFile = open(diagFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % diagFileName)
                
    try:
        errorFile = open(errorFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % errorFileName)
                
    try:
        inputFile = open(inputFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inputFileName)

    try:
        noteFile = open(noteFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (db.get_sqlServer()))
    diagFile.write('Database: %s\n' % (db.get_sqlDatabase()))

    errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    # delete existing note
    db.sql('delete from MGI_Note n where n._mgitype_key = 11 and n._notetype_key = 1053', None)
    db.commit()

#
# Purpose: Close files.
#
def closeFiles():

    noteFile.close()

#
# Purpose:  sets global primary key variables
#
def setPrimaryKeys():

    global noteKey

    results = db.sql(''' select nextval('mgi_note_seq') as maxKey ''', 'auto')
    noteKey = results[0]['maxKey']

#
# Purpose:  BCPs the data into the database
#
def bcpFiles():

    closeFiles()

    bcpI = '%s %s %s' % (BCP_COMMAND, db.get_sqlServer(), db.get_sqlDatabase())
    bcpII = '"|" "\\n" mgd'
    bcpCmd = '%s %s "/" %s %s' % (bcpI, noteTable, noteFileName, bcpII)
    diagFile.write('%s\n' % bcpCmd)
    os.system(bcpCmd)

    # update mgi_note_seq auto-sequence
    db.sql(''' select setval('mgi_note_seq', (select max(_Note_key) from MGI_Note)) ''', None)
    db.commit()

#
# Purpose:  processes data
#
def processFile():

    global noteKey

    lineNum = 0
    # For each line in the input file

    totalRows = 0
    totalSkipped = 0
    totalProcessed = 0

    # save markerID:alleleID duplicate set
    allSet = []
    dupSet = []
    for line in inputFile.readlines():

        tokens = line[:-1].split('\t')

        if tokens[0] == '' or tokens[0] == 'Gene Symbol':
            continue

        markerID = tokens[2]
        alleleID = tokens[3]

        if alleleID == '':
            continue

        id = markerID + '|' + alleleID
        if id in allSet:
            dupSet.append(id)
        else:
            allSet.append(id)

    # reset pointer to the start
    inputFile.seek(0)

    for line in inputFile.readlines():

        error = 0
        lineNum  += 1

        # Split the line into tokens
        tokens = line[:-1].split('\t')

        if tokens[0] == '' or tokens[0] == 'Gene Symbol':
            continue

        totalRows += 1

        try:
            markerID = tokens[2]
            alleleID = tokens[3]
            note = tokens[15]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

        # skip if there is no allele
        if alleleID == '':
            totalSkipped += 1
            continue

        # skip if duplicate
        id = markerID + '|' + alleleID
        if id in dupSet:
            print('Duplicate Lines in input file: ', id)
            totalSkipped += 1
            continue

        markerKey = loadlib.verifyMarker(markerID, lineNum, errorFile)
        alleleKey = loadlib.verifyObject(alleleID, 11, None, lineNum, errorFile)
        
        if markerKey == 0:
            print('The MGI Gene ID is not a valid MGI marker ID: ', markerID)
            totalSkipped += 1
            continue

        if alleleKey == None:
            print('The MGI Allele ID is not a valid MGI allele ID: ', alleleID)
            totalSkipped += 1
            continue

        totalProcessed += 1

        results = db.sql(''' 
            select count(*) from all_allele where _allele_key = %s and _marker_key = %s 
            ''' % (markerKey, alleleKey), 'auto')
        if len(results) == 0:
            print('The MGI Allele ID is not associated with the MGI Gene ID in the same row: ', markerID, alleleID)

        # if no errors, process the allele

        note += ' (<A href=https://www.informatics.jax.org/reference/J:384794>J:384794</A>)'
        noteFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (noteKey, alleleKey, mgiTypeKey, noteTypeKey, \
               note, createdByKey, createdByKey, loaddate, loaddate))

        noteKey = noteKey + 1

    print('rows input: ', str(totalRows))
    print('rows skipped: ', str(totalSkipped))
    print('rows processed: ', str(totalProcessed))

    #	end of "for line in inputFile.readlines():"

#
# Main
#

if __name__ == '__main__':
        print('initialize')
        initialize()
        print('setPrimaryKeys')
        setPrimaryKeys()
        print('processFile')
        processFile()
        print('bcpFiles')
        bcpFiles()

        exit(0)
