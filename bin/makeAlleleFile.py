#!/usr/local/bin/python

#
# Program: makeAlleleFile.py
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
#	makeAlleleFile.py
#
# Envvars:
#
# Inputs:
#
#	A tab-delimited file in the format:
#
#	field 1:  MGI Marker ID
#	field 2:  Allele Symbol
#	field 3:  Allele Name
#	field 4:  Allele Status
#	field 5:  Allele Type
#	field 6:  Allele Attribute/SubType
#	field 7:  Germ Line Transmission
#	field 8:  Reference Type/J#
#	field 9:  Strain of Origin
#	field 10: Mutant Cell Line ID
#	field 11: Molecular Notes
#	field 12: Driver Notes
#	field 13: Molecular Mutation
#	field 14: Inheritance
#	field 15: Mixed
#	field 16: Extinct
#	field 17: Creation Date
#
# Outputs:
#
#       BCP files:
#
#       ALL_Allele.bcp                  master Allele records
#	ALL_Marker_Assoc.bcp		allele/marker associations
#	ALL_Allele_Mutation.bcp
#       ALL_Allele_CellLine.bcp
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
# 01/27/2014	lec
#	- TR11515/Sanger
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
mode = os.environ['ALLELEMODE']
inputFileName = os.environ['ALLELEDATAFILE']
outputDir = os.environ['ALLELEDATADIR']
jnum = os.environ['NOMENREFERENCE']

DEBUG = 0		# if 0, not in debug mode
TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline

bcpon = 1		# can the bcp files be bcp-ed into the database?  default is yes.

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor
inputFile = ''		# file descriptor
alleleFile = ''         # file descriptor
markerFile = ''		# file descriptor
mutationFile = ''	# file descriptor
mutantFile = ''		# file descriptor
refFile = ''            # file descriptor
accFile = ''            # file descriptor
accRefFile = ''         # file descriptor
noteFile = ''		# file descriptor
noteChunkFile = ''	# file descriptor

alleleTable = 'ALL_Allele'
markerTable = 'ALL_Marker_Assoc'
mutationTable = 'ALL_Allele_Mutation'
mutantTable = 'ALL_Allele_CellLine'
refTable = 'MGI_Reference_Assoc'
accTable = 'ACC_Accession'
accRefTable = 'ACC_AccessionReference'
noteTable = 'MGI_Note'
noteChunkTable = 'MGI_NoteChunk'
newAlleleFile = 'newAllele.txt'

alleleFileName = outputDir + '/' + alleleTable + '.bcp'
markerFileName = outputDir + '/' + markerTable + '.bcp'
mutationFileName = outputDir + '/' + mutationTable + '.bcp'
mutantFileName =  outputDir + '/' + mutantTable + '.bcp'
refFileName = outputDir + '/' + refTable + '.bcp'
accFileName = outputDir + '/' + accTable + '.bcp'
accRefFileName = outputDir + '/' + accRefTable + '.bcp'
noteFileName = outputDir + '/' + noteTable + '.bcp'
noteChunkFileName = outputDir + '/' + noteChunkTable + '.bcp'
newAlleleFileName = outputDir + '/' + newAlleleFile

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name

alleleKey = 0           # ALL_Allele._Allele_key
assocKey =  0		# ALL_Marker_Assoc._Assoc_key
mutantKey = 0  		# ALL_Allele_CellLine._Assoc_key
refAssocKey = 0		# MGI_Reference_Assoc._Assoc_key
accKey = 0              # ACC_Accession._Accession_key
noteKey = 0		# MGI_Note._Note_key
mgiKey = 0              # ACC_AccessionMax.maxNumericPart
mgiNoteObjectKey = 11   # MGI_Note._MGIType_key
mgiNoteSeqNum = 1       # MGI_NoteChunk.sequenceNum
mgiMolecularNoteTypeKey = 1021   # MGI_Note._NoteType_key
mgiDriverNoteTypeKey = 1034   	 # MGI_Note._NoteType_key

NA = -2			# for Not Applicable fields
mgiTypeKey = 11		# Allele
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
    global alleleFile, markerFile, mutationFile, mutantFile, refFile, 
    global accFile, accRefFile, noteFile, noteChunkFile
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
        mutationFile = open(mutationFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % mutationFileName)

    try:
        mutantFile = open(mutantFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % mutantFileName)

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
        noteChunkFile = open(noteChunkFileName, 'w')
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

    global alleleKey, assocKey, refAssocKey, accKey, noteKey, mgiKey

    results = db.sql('select maxKey = max(_Allele_key) + 1 from ALL_Allele', 'auto')
    alleleKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Assoc_key) + 1 from ALL_Marker_Assoc', 'auto')
    assocKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Assoc_key) + 1 from MGI_Reference_Assoc', 'auto')
    refAssocKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Accession_key) + 1 from ACC_Accession', 'auto')
    accKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Note_key) + 1 from MGI_Note', 'auto')
    noteKey = results[0]['maxKey']

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
    mutationFile.close()
    mutantFile.close()
    refFile.close()
    accFile.close()
    accRefFile.close()
    noteFile.close()
    noteChunkFile.close()

    bcpI = 'cat %s | bcp %s..' % (passwordFileName, db.get_sqlDatabase())
    bcpII = '-c -t\"|" -S%s -U%s' % (db.get_sqlServer(), db.get_sqlUser())

    bcp1 = '%s%s in %s %s' % (bcpI, alleleTable, alleleFileName, bcpII)
    bcp2 = '%s%s in %s %s' % (bcpI, markerTable, markerFileName, bcpII)
    bcp3 = '%s%s in %s %s' % (bcpI, mutationTable, mutationFileName, bcpII)
    bcp4 = '%s%s in %s %s' % (bcpI, mutantTable, mutantFileName, bcpII)
    bcp5 = '%s%s in %s %s' % (bcpI, refTable, refFileName, bcpII)
    bcp6 = '%s%s in %s %s' % (bcpI, accTable, accFileName, bcpII)
    bcp7 = '%s%s in %s %s' % (bcpI, accRefTable, accRefFileName, bcpII)
    bcp8 = '%s%s in %s %s' % (bcpI, noteTable, noteFileName, bcpII)
    bcp9 = '%s%s in %s %s' % (bcpI, noteChunkTable, noteChunkFileName, bcpII)

    for bcpCmd in [bcp1, bcp2, bcp3, bcp4, bcp5, bcp6, bcp7, bcp8, bcp9]:
	diagFile.write('%s\n' % bcpCmd)
	os.system(bcpCmd)

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processFile():

    global alleleKey, assocKey, refAssocKey, accKey, noteKey, mgiKey

    lineNum = 0
    # For each line in the input file

    for line in inputFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = line[:-1].split('\t')

	#	field 1:  MGI Marker ID
	#	field 2:  Allele Symbol
	#	field 3:  Allele Name
	#	field 4:  Allele Status
	#	field 5:  Allele Type
	#	field 6:  Germ Line Transmission
	#	field 7:  Reference Type/J#
	#	field 8:  Strain of Origin
	#	field 9:  Mutant Cell Line ID
	#	field 10: Molecular Notes
	#	field 11: Driver Notes
	#	field 12: Molecular Mutation
	#	field 13: Creation Date

        try:
	    markerID = tokens[0]
	    symbol = tokens[1]
	    name = tokens[2]
	    alleleStatus = tokens[3]
	    alleleType = tokens[4]
	    germLine = tokens[5]
	    references = tokens[6]
	    strainOfOrigin = tokens[7]
	    mutantCellLine = tokens[8]
	    molecularNotes = tokens[9]
	    driverNotes = tokens[10]
	    mutation = tokens[11]
	    createdBy = tokens[12]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	# marker key
	markerKey = loadlib.verifyMarker(markerID, lineNum, errorFile)

	# hard-coded
	# _vocab_key = 70 (Marker-Allele Association Qualifier)
	# _term_key = 4268547 (Not Specified)
	qualifierKey = 4268547

	# hard-coded
	# _vocab_key = 73 (Marker-Allele Association Status)
	# _term_key = 4268545 (Curated)
	markerStatusKey = 4268545

	# hard-coded
	# _vocab_key = 35 (Allele Inheritance Mode)
	# _term_key = 847095 (Not Applicable)
	modeKey = 847095

	#
	# select * from MGI_RefAssocType 
	# where _MGIType_key = 11
	# and accType = 'Original'
	#

	# Original|J:xxxx
	refKey = loadlib.verifyReference(jnum, lineNum, errorFile)
	refAssocTypeKey = 1011

	# _vocab_key = 37 (Allele Status)
	alleleStatusKey loadlib.verifyTerm('', 37, alleleStatus, lineNum, errorFile)

	# _vocab_key = 38 (Allele Type)
	alleleTypeKey = loadlib.verifyTerm('', 38, alleleType, lineNum, errorFile)

	# _vocab_key = 61 (Allele Transmission)
	germLineKey = loadlib.verifyTerm('', 61, germLine, lineNum, errorFile)

	# _vocab_key = 36 (Allele Molecular Mutation)
	mutationKey = loadlib.verifyTerm('', 36, mutation, lineNum, errorFile)

	createdByKey = loadlib.verifyUser(createdBy, lineNum, errorFile)

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process the allele

        alleleFile.write('%d|%s|%s|%s|%s|%s|%s|%s|%s||0|0|0|%s|%s|%s|%s|%s|%s\n' \
            % (alleleKey, markerKey, strainOfOriginKey, modeKey, alleleTypeKey, \
	    alleleStatusKey, germLineKey, symbol, name, \
	    createdByKey, createdByKey, createdByKey, loaddate, loaddate, loaddate))

        markerFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
	    % (assocKey, alleleKey, markerKey, qualifierKey, refKey, markerStatusKey, \
	       createdByKey, createdByKey, loaddate, loaddate))

        mutationFile.write('%s|%s|%s|%s\n' \
	    % (alleleKey, mutationKey, loaddate, loaddate))

        refFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
	    % (refAssocKey, refKey, alleleKey, mgiTypeKey, refAssocTypeKey, \
	       createdByKey, createdByKey, loaddate, loaddate))

        #
        # mutant cell line
        #
        if len(mutantCellLine) > 0):
            addMutantCellLine(alleleKey, mutantCellLine, createdByKey)

        # MGI Accession ID for the allelearker

        accFile.write('%s|%s%d|%s|%s|1|%d|%d|0|1|%s|%s|%s|%s\n' \
            % (accKey, mgiPrefix, mgiKey, mgiPrefix, mgiKey, alleleKey, mgiTypeKey, \
	       createdByKey, createdByKey, loaddate, loaddate))

	# storing data in MGI_Note/MGI_NoteChunk
	# molecular notes

	mgiNoteSeqNum = 1
	if len(molecularNotes) > 0:

	    noteFile.write('%s|%s|%s|%s|%s|%s|%s|%s\n' \
		% (noteKey, alleleKey, mgiNoteObjectKey, mgiMolecularNoteTypeKey, \
		   createdByKey, createdByKey, loaddate, loaddate))

	    while len(molecularNotes) > 255:
	        noteChunkFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
		    % (noteKey, mgiNoteSeqNum, molecularNotes[:255], createdByKey, createdByKey, loaddate, loaddate))
		molecularNotes = molecularNotes[255:]
		mgiNoteSeqNum = mgiNoteSeqNum + 1

	    if len(molecularNotes) > 0:
	        noteChunkFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
		    % (noteKey, mgiNoteSeqNum, molecularNotes, createdByKey, createdByKey, loaddate, loaddate))

	    noteKey = noteKey + 1

	# driver notes
	mgiNoteSeqNum = 1
	if len(driverNotes) > 0:

	    noteFile.write('%s|%s|%s|%s|%s|%s|%s|%s\n' \
		% (noteKey, alleleKey, mgiNoteObjectKey, mgiDriverNoteTypeKey, \
		   createdByKey, createdByKey, loaddate, loaddate))

	    while len(driverNotes) > 255:
	        noteChunkFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
		    % (noteKey, mgiNoteSeqNum, driverNotes[:255], createdByKey, createdByKey, loaddate, loaddate))
		driverNotes = driverNotes[255:]
		mgiNoteSeqNum = mgiNoteSeqNum + 1

	    if len(driverNotes) > 0:
	        noteChunkFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
		    % (noteKey, mgiNoteSeqNum, driverNotes, createdByKey, createdByKey, loaddate, loaddate))

	    noteKey = noteKey + 1

	# Print out a new text file and attach the new MGI Allele IDs as the last field

        newAlleleFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s\t%s\t%s\n' \
	    % (markerID, \
	    symbol, \
	    name, \
	    alleleStatus, \
	    alleleType, \
	    mgi_utils.prvalue(germLine), \
	    references, \
	    mgi_utils.prvalue(strainOfOrigin), \
	    mgi_utils.prvalue(molecularNotes), \
	    mgi_utils.prvalue(driverNotes), \
	    createdBy, mgiPrefix, mgiKey))

        accKey = accKey + 1
        mgiKey = mgiKey + 1
	refAssocKey = refAssocKey + 1
	assocKey = assocKey + 1
        alleleKey = alleleKey + 1

    #	end of "for line in inputFile.readlines():"

    #
    # Update the AccessionMax value
    #

    if not DEBUG:
        db.sql('exec ACC_setMax %d' % (lineNum), None)

def addMutantCellLine(alleleKey, mutantCellLine, createdByKey):

    global mutantKey

    mutantCellLineKey = 0

    results = db.sql('''
        select c._CellLine_key
        from ALL_CellLine c, ALL_CellLine_Derivation d
        where c.isMutant = 1
        and c._Derivation_key = d._Derivation_key
	and c.cellLine = '%s'
        ''' % (mutantCellLine) , 'auto')

    for r in results:
        mutantCellLineKey = r['_CellLine_key']

    mutantFile.write('%d|%s|%s|%s|%s|%s|%s\n' \
            % (mutantKey, alleleKey, mutantCellLineKey, \
               createdByKey, createdByKey, loaddate, loaddate))

    mutantKey = mutantKey + 1

#
# Main
#

init()
verifyMode()
setPrimaryKeys()
processFile()
bcpFiles()
exit(0)

