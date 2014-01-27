#!/usr/local/bin/python

#
# Program: alleleload.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To load new Alleles into MGI
#
# Requirements Satisfied by This Program:
#
# Usage:
#	alleleload.py
#
# Envvars:
#
# Inputs:
#
#	A tab-delimited file in the format:
#		field 1:  Allele Symbol
#		field 2:  Allele Name
#		field 3:  MGI Marker ID
#		field 4:  Allele Type
#		field 5:  Strain of Origin
#	   	field 6:  Reference (original, marker)
#		field 7:  Molecular Mutations(s)	xxxx|xxxx|...
#	   	field 8:  Other References (transmission, mixed, molecular, ...)?
#		field 9:  Allele Synonyms (optional)
#		field 10: General Notes (optional)
#		field 11: Molecular Notes (optional)
#		field 12: Created By
#		field 13: Creation Date
#
# Outputs:
#
#       BCP files:
#
#       ALL_Allele.bcp                  master Allele records
#	ALL_Allele_CellLine.bcp         allele/cellline associations
#	ALL_Allele_Mutation.bcp         allele/mutation associations
#	ALL_Marker_Assoc.bcp		allele/marker associations
#
#       MGI_Reference_Assoc             allele/reference associations (all types)
#       MGI_Note/MGI_NoteChunk          allele notes (all types)
#
#       ACC_Accession.bcp               Accession records
#       ACC_AccessionReference.bcp      Accession Reference records
#
#       Diagnostics file of all input parameters and SQL commands
#       Error file
#
# Exit Codes:
#
# Assumes:
#
#	That no one else is adding such records to the database.
#
# Bugs:
#
# Implementation:
#
# History
#
# 11/24/2010	lec
#	- TR10267/Gensat Transgene Alleles
#

import sys
import os
import accessionlib
import db
import mgi_utils
import loadlib

#globals

#
# from configuration file
#
user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']
mode = os.environ['ALLELELOADMODE']
inputFileName = os.environ['ALLELEDATAFILE']
outputDir = os.environ['ALLELELOADDATADIR']

DEBUG = 0		# if 0, not in debug mode
TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline

bcpon = 1		# can the bcp files be bcp-ed into the database?  default is yes.

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor
inputFile = ''		# file descriptor
alleleFile = ''         # file descriptor
markerFile = ''		# file descriptor
refFile = ''            # file descriptor
accFile = ''            # file descriptor
accRefFile = ''         # file descriptor
noteFile = ''		# file descriptor
noteChunkFile = ''	# file descriptor

alleleTable = 'ALL_Allele'
refTable = 'MGI_Reference_Assoc'
accTable = 'ACC_Accession'
accRefTable = 'ACC_AccessionReference'
noteTable = 'MGI_Note'
noteChunkTable = 'MGI_NoteChunk'
newAlleleFile = 'newAllele.txt'

alleleFileName = outputDir + '/' + alleleTable + '.bcp'
markerFileName = outputDir + '/' + markerTable + '.bcp'
refFileName = outputDir + '/' + refTable + '.bcp'
accFileName = outputDir + '/' + accTable + '.bcp'
accRefFileName = outputDir + '/' + accRefTable + '.bcp'
noteFileName = outputDir + '/' + noteTable + '.bcp'
noteChunkFileName = outputDir + '/' + noteChunkTable + '.bcp'
newAlleleFileName = outputDir + '/' + newAlleleFile

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name

alleleKey = 0            # PRB_Probe._Probe_key
refKey = 0		# PRB_Reference._Reference_key
accKey = 0              # ACC_Accession._Accession_key
mgiKey = 0              # ACC_AccessionMax.maxNumericPart

NA = -2			# for Not Applicable fields
mgiTypeKey = 3		# Molecular Segment
mgiPrefix = "MGI:"

loaddate = loadlib.loaddate

# Purpose: prints error message and exits
# Returns: nothing
# Assumes: nothing
# Effects: exits with exit status
# Throws: nothing

def exit(
    status,          # numeric exit status (integer)
    message = None   # exit message (string)
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
        pass

    db.useOneConnection(0)
    sys.exit(status)
 
# Purpose: process command line options
# Returns: nothing
# Assumes: nothing
# Effects: initializes global variables
#          exits if files cannot be opened
# Throws: nothing

def init():
    global diagFile, errorFile, inputFile, errorFileName, diagFileName
    global alleleFile, markerFile, refFile, accFile, accRefFile, noteFile, noteChunkFile
    global newAlleleFile
 
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)
 
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
        alleleFile = open(alleleFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % alleleFileName)

    try:
        markerFile = open(markerFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % markerFileName)

    try:
        refFile = open(refFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % refFileName)

    try:
        accFile = open(accFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % accFileName)

    try:
        accRefFile = open(accRefFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % accRefFileName)

    try:
        noteFile = open(noteFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteFileName)

    try:
        noteChunkFile = open(noteFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteChunkFileName)

    try:
        newAlleleFile = open(newAlleleFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % newAlleleFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    # Set Log File Descriptor
    db.set_sqlLogFD(diagFile)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (db.get_sqlServer()))
    diagFile.write('Database: %s\n' % (db.get_sqlDatabase()))

    errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    return

# Purpose: verify processing mode
# Returns: nothing
# Assumes: nothing
# Effects: if the processing mode is not valid, exits.
#	   else, sets global variables
# Throws:  nothing

def verifyMode():

    global DEBUG

    if mode == 'preview':
        DEBUG = 1
        bcpon = 0
    elif mode != 'load':
        exit(1, 'Invalid Processing Mode:  %s\n' % (mode))

# Purpose:  sets global primary key variables
# Returns:  nothing
# Assumes:  nothing
# Effects:  sets global primary key variables
# Throws:   nothing

def setPrimaryKeys():

    global alleleKey, refKey, accKey, mgiKey

    results = db.sql('select maxKey = max(_Probe_key) + 1 from PRB_Probe', 'auto')
    alleleKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Reference_key) + 1 from PRB_Reference', 'auto')
    refKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Accession_key) + 1 from ACC_Accession', 'auto')
    accKey = results[0]['maxKey']

    results = db.sql('select maxKey = maxNumericPart + 1 from ACC_AccessionMax ' + \
        'where prefixPart = "%s"' % (mgiPrefix), 'auto')
    mgiKey = results[0]['maxKey']

# Purpose:  BCPs the data into the database
# Returns:  nothing
# Assumes:  nothing
# Effects:  BCPs the data into the database
# Throws:   nothing

def bcpFiles():

    bcpdelim = "|"

    if DEBUG or not bcpon:
        return

    alleleFile.close()
    markerFile.close()
    refFile.close()
    accFile.close()
    accRefFile.close()
    noteFile.close()

    bcpI = 'cat %s | bcp %s..' % (passwordFileName, db.get_sqlDatabase())
    bcpII = '-c -t\"|" -S%s -U%s' % (db.get_sqlServer(), db.get_sqlUser())

    bcp1 = '%s%s in %s %s' % (bcpI, alleleTable, alleleFileName, bcpII)
    bcp2 = '%s%s in %s %s' % (bcpI, markerTable, markerFileName, bcpII)
    bcp3 = '%s%s in %s %s' % (bcpI, refTable, refFileName, bcpII)
    bcp5 = '%s%s in %s %s' % (bcpI, accTable, accFileName, bcpII)
    bcp6 = '%s%s in %s %s' % (bcpI, accRefTable, accRefFileName, bcpII)
    bcp7 = '%s%s in %s %s' % (bcpI, noteTable, noteFileName, bcpII)

    for bcpCmd in [bcp1, bcp2, bcp3, bcp4, bcp5, bcp6, bcp7]:
	diagFile.write('%s\n' % bcpCmd)
	os.system(bcpCmd)

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processFile():

    global alleleKey, refKey, accKey, mgiKey

    lineNum = 0
    # For each line in the input file

    for line in inputFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = line[:-1].split('\t')

        try:
	    name = tokens[0]
	    jnum = tokens[1]
	    parentID = tokens[2]
	    sourceName = tokens[3]
	    organism = tokens[4]
	    strain = tokens[5]
	    tissue = tokens[6]
	    gender = tokens[7]
	    cellLine = tokens[8]
	    age = tokens[9]
	    vectorType = tokens[10]
	    segmentType = tokens[11]
	    regionCovered = tokens[12]
	    insertSite = tokens[13]
	    insertSize = tokens[14]
	    markerIDs = tokens[15].split('|')
	    relationship = tokens[16]
	    sequenceIDs = tokens[17]
	    notes = tokens[19]
	    createdBy = tokens[20]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

        referenceKey = loadlib.verifyReference(jnum, lineNum, errorFile)
	createdByKey = loadlib.verifyUser(createdBy, lineNum, errorFile)

	if referenceKey == 0:
	    errorFile.write('Invalid Reference:  %s\n' % (jnum))
	    error = 1

	if createdByKey == 0:
	    errorFile.write('Invalid Creator:  %s\n\n' % (createdBy))
	    error = 1

	# marker IDs

	markerList = []
	for markerID in markerIDs:

	    markerKey = loadlib.verifyMarker(markerID, lineNum, errorFile)

	    if len(markerID) > 0 and markerKey == 0:
	        errorFile.write('Invalid Marker:  %s, %s\n' % (name, markerID))
	        error = 1
            elif len(markerID) > 0:
		markerList.append(markerKey)

	# sequence IDs
	seqAccDict = {}
	for seqID in sequenceIDs.split('|'):
	    if len(seqID) > 0:
	        [logicalDB, acc] = seqID.split(':')
	        logicalDBKey = loadlib.verifyLogicalDB(logicalDB, lineNum, errorFile)
	        if logicalDBKey > 0:
		    seqAccDict[acc] = logicalDBKey

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process the probe

        alleleFile.write('%d|%s|%s|%s|%s|%s|||%s|%s|%s||%s|%s|%s|%s\n' \
            % (alleleKey, name, parentProbeKey, sourceKey, vectorKey, segmentTypeKey, mgi_utils.prvalue(regionCovered), \
	    mgi_utils.prvalue(insertSite), mgi_utils.prvalue(insertSize), createdByKey, createdByKey, loaddate, loaddate))

	for markerKey in markerList:
	    if markerList.count(markerKey) == 1:
                markerFile.write('%s|%s|%d|%s|%s|%s|%s|%s\n' \
		    % (alleleKey, markerKey, referenceKey, relationship, createdByKey, createdByKey, loaddate, loaddate))
            else:
		errorFile.write('Invalid Marker Duplicate:  %s, %s\n' % (name, markerID))

        refFile.write('%s|%s|%s|0|0|%s|%s|%s|%s\n' \
		% (refKey, alleleKey, referenceKey, createdByKey, createdByKey, loaddate, loaddate))

        # MGI Accession ID for the marker

        accFile.write('%s|%s%d|%s|%s|1|%d|%d|0|1|%s|%s|%s|%s\n' \
            % (accKey, mgiPrefix, mgiKey, mgiPrefix, mgiKey, alleleKey, mgiTypeKey, createdByKey, createdByKey, loaddate, loaddate))

	# Print out a new text file and attach the new MGI Probe IDs as the last field

        newAlleleFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%d\n' \
	    % (name, jnum, \
	    mgi_utils.prvalue(sourceName), \
	    organism, \
	    mgi_utils.prvalue(strain), \
	    mgi_utils.prvalue(tissue), \
	    mgi_utils.prvalue(gender), \
	    mgi_utils.prvalue(cellLine), \
	    mgi_utils.prvalue(age), \
	    mgi_utils.prvalue(vectorType), \
	    mgi_utils.prvalue(segmentType), \
	    mgi_utils.prvalue(regionCovered) + \
	    mgi_utils.prvalue(insertSite), \
	    mgi_utils.prvalue(insertSize), \
	    markerIDs.join('|'), \
	    relationship, \
	    mgi_utils.prvalue(sequenceIDs), \
	    mgi_utils.prvalue(notes), \
	    createdBy, mgiPrefix, mgiKey))

	# Notes

        noteSeq = 1
		
        while len(notes) > 255:
	    noteFile.write('%s|%d|%s|%s|%s\n' % (alleleKey, noteSeq, notes[:255], loaddate, loaddate))
            newnote = notes[255:]
            notes = newnote
            noteSeq = noteSeq + 1

        if len(notes) > 0:
	    noteFile.write('%s|%d|%s|%s|%s\n' % (alleleKey, noteSeq, notes, loaddate, loaddate))

        accKey = accKey + 1
        mgiKey = mgiKey + 1

	# sequence accession ids
	for acc in seqAccDict.keys():
	    prefixPart, numericPart = accessionlib.split_accnum(acc)
            accFile.write('%s|%s|%s|%s|%s|%d|%d|0|1|%s|%s|%s|%s\n' \
                % (accKey, acc, prefixPart, numericPart, seqAccDict[acc], alleleKey, mgiTypeKey, createdByKey, createdByKey, loaddate, loaddate))
            accRefFile.write('%s|%s|%s|%s|%s|%s\n' \
                % (accKey, referenceKey, createdByKey, createdByKey, loaddate, loaddate))
	    accKey = accKey + 1

	refKey = refKey + 1
        alleleKey = alleleKey + 1

    #	end of "for line in inputFile.readlines():"

    #
    # Update the AccessionMax value
    #

    if not DEBUG:
        db.sql('exec ACC_setMax %d' % (lineNum), None)

#
# Main
#

init()
verifyMode()
setPrimaryKeys()
processFile()
bcpFiles()
exit(0)

