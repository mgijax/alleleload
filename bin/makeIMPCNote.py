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
duplicateFileName = os.environ['IMPC_DUPLICATE']
outputDir = os.environ['OUTPUTDIR']
BCP_COMMAND = os.environ['PG_DBUTILS'] + '/bin/bcpin.csh'

inputFile = ''		# file descriptor
duplicateFile = ''	# file descriptor
noteFile = ''		# file descriptor

noteTable = 'MGI_Note'
noteFileName = outputDir + '/' + noteTable + '.bcp'

duplicateSet = []

mgiTypeKey = 11
noteTypeKey = 1053
createdByKey = 1000

loaddate = loadlib.loaddate

#
# Purpose: process command line options
#
def initialize():
    global inputFile, duplicateFile
    global noteFile
    global noteKey
    global duplicateSet
    
    head, tail = os.path.split(inputFileName)

    try:
        inputFile = open(inputFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inputFileName)

    try:
        duplicateFile = open(duplicateFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % duplicateFileName)

    try:
        noteFile = open(noteFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    # duplicates
    duplicateSet = []
    for line in duplicateFile.readlines():
        tokens = line[:-1].split('\t')
        if tokens[1] == '':
            continue
        duplicateSet.append(tokens[0] + '|' + tokens[1])
    duplicateFile.close()

    # delete existing note
    db.sql('delete from MGI_Note n where n._mgitype_key = 11 and n._notetype_key = 1053', None)
    db.commit()

    # set max key
    results = db.sql(''' select nextval('mgi_note_seq') as maxKey ''', 'auto')
    noteKey = results[0]['maxKey']

#
# Purpose:  BCPs the data into the database
#
def bcpFiles():

    bcpI = '%s %s %s' % (BCP_COMMAND, db.get_sqlServer(), db.get_sqlDatabase())
    bcpII = '"|" "\\n" mgd'
    bcpCmd = '%s %s "/" %s %s' % (bcpI, noteTable, noteFileName, bcpII)
    print('%s\n' % bcpCmd)
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
        if id in duplicateSet:
            print('Duplicate Lines in input file: ', id)
            totalSkipped += 1
            continue

        markerKey = loadlib.verifyMarker(markerID, lineNum, None)
        alleleKey = loadlib.verifyObject(alleleID, 11, None, lineNum, None)
        
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
            select * from all_allele where _allele_key = %s and _marker_key = %s 
            ''' % (alleleKey, markerKey), 'auto')
        if len(results) == 0:
            print('The MGI Allele ID is not associated with the MGI Gene ID in the same row: ', markerID, alleleID)

        # if no errors, process the allele

        note += ' (<A href="https://www.informatics.jax.org/reference/J:384794">J:384794</A>)'
        noteFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (noteKey, alleleKey, mgiTypeKey, noteTypeKey, \
               note, createdByKey, createdByKey, loaddate, loaddate))

        noteKey = noteKey + 1

    print('\n')
    print('rows input: ', str(totalRows))
    print('rows skipped: ', str(totalSkipped))
    print('rows processed: ', str(totalProcessed))

    #	end of "for line in inputFile.readlines():"
    inputFile.close()
    noteFile.close()

#
# Main
#

if __name__ == '__main__':
        print('\ninitialize')
        initialize()
        print('\nprocessFile')
        processFile()
        print('\nbcpFiles')
        bcpFiles()

