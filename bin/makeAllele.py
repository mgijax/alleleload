#!/usr/local/bin/python

#
# Program: makeAllele.py
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
#	makeAllele.py
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
#	field 5:  Allele Generation (Type)
#	field 6:  Allele Subtype (currently not used)
#	field 7:  Allele Collection (currently not used)
#	field 8:  Germ Line Transmission
#	field 9:  Reference Type/J#
#	field 10: Strain of Origin
#	field 11: Mutant Cell Line ID
#	field 12: Molecular Notes (_NoteType_key = 1021)
#	field 13: Driver Notes (_NoteType_key = 1034)
#	field 14: IKMC Colony Name (_NoteType_key = 1041)
#	field 15: Molecular Mutation
#	field 16: Inheritance Mode
#	field 17: Mixed
#	field 18: Extinct
#	field 19: Creation Date
#	field 20: Add mutant cell line/IKMC to new Allele (not yet in database)
#	field 21: Add mutant cell line/IKMC to existing Allele (in database)
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
#	New file (input + new accession ids)
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
#	- TR11515/IKMC
#
# 11/24/2010	lec
#	- TR10267/Gensat Transgene Alleles
#

import sys
import os
import db
import mgi_utils
import loadlib
import sourceloadlib

#globals

#
# from configuration file
#
user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']
inputFileName = os.environ['INPUTFILE']
outputDir = os.environ['OUTPUTDIR']
jnum = os.environ['JNUMBER']

DEBUG = 0		# if 0, not in debug mode

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
newAlleleFile = ''      # file descriptor

alleleTable = 'ALL_Allele'
markerTable = 'ALL_Marker_Assoc'
mutationTable = 'ALL_Allele_Mutation'
mutantTable = 'ALL_Allele_CellLine'
refTable = 'MGI_Reference_Assoc'
accTable = 'ACC_Accession'
accRefTable = 'ACC_AccessionReference'
noteTable = 'MGI_Note'
noteChunkTable = 'MGI_NoteChunk'

alleleFileName = outputDir + '/' + alleleTable + '.bcp'
markerFileName = outputDir + '/' + markerTable + '.bcp'
mutationFileName = outputDir + '/' + mutationTable + '.bcp'
mutantFileName =  outputDir + '/' + mutantTable + '.bcp'
refFileName = outputDir + '/' + refTable + '.bcp'
accFileName = outputDir + '/' + accTable + '.bcp'
accRefFileName = outputDir + '/' + accRefTable + '.bcp'
noteFileName = outputDir + '/' + noteTable + '.bcp'
noteChunkFileName = outputDir + '/' + noteChunkTable + '.bcp'

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name
newAlleleFileName = ''	# output file with new accession ids

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
mgiIKMCNoteTypeKey = 1041   	 # MGI_Note._NoteType_key
ikmcNoteSQL = ''

mgiTypeKey = 11		# Allele
mgiPrefix = "MGI:"

loaddate = loadlib.loaddate

#
# Purpose: prints error message and exits
#
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
	newAlleleFile.close()
    except:
        pass

    db.useOneConnection(0)
    sys.exit(status)
 
#
# Purpose: process command line options
#
def initialize():
    global diagFile, errorFile, inputFile, errorFileName, diagFileName
    global alleleFile, markerFile, mutationFile, mutantFile, refFile
    global accFile, accRefFile, noteFile, noteChunkFile
    global newAlleleFile
 
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)
 
    head, tail = os.path.split(inputFileName) 

    diagFileName = outputDir + '/' + tail + '.diagnostics'
    errorFileName = outputDir + '/' + tail + '.error'
    newAlleleFileName = outputDir + '/' + tail + '.new'

    try:
        diagFile = open(diagFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % diagFileName)
		
    try:
        errorFile = open(errorFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % errorFileName)
		
    try:
        newAlleleFile = open(newAlleleFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % newAlleleFileName)

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

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    # Set Log File Descriptor
    db.set_sqlLogFD(diagFile)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (db.get_sqlServer()))
    diagFile.write('Database: %s\n' % (db.get_sqlDatabase()))

    errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    return 0

#
# Purpose: Close files.
#
def closeFiles():

    alleleFile.close()
    markerFile.close()
    mutationFile.close()
    mutantFile.close()
    refFile.close()
    accFile.close()
    accRefFile.close()
    noteFile.close()
    noteChunkFile.close()

    return 0

#
# Purpose:  sets global primary key variables
#
def setPrimaryKeys():

    global alleleKey, assocKey, refAssocKey, accKey, noteKey, mgiKey, mutantKey

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

    results = db.sql('select maxKey = max(_Assoc_key) + 1 from ALL_Allele_CellLine', 'auto')
    mutantKey = results[0]['maxKey']

    return 0

#
# Purpose:  BCPs the data into the database
#
def bcpFiles():

    bcpdelim = "|"

    if DEBUG or not bcpon:
        return

    closeFiles()

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

    #db.sql(ikmcNoteSQL, None)

    return 0

#
# Purpose:  processes data
#
def processFileIKMC(createMutantCellLine1, createMutantCellLine2, \
	mutantCellLine, ikmcNotes, createdByKey):

    global noteKey, ikmcNoteSQL

    print mutantCellLine

    if len(createMutantCellLine2) > 0:

	tokens = createMutantCellLine2.split(':')
	alleleKey = tokens[0]
	newNotes = tokens[1]

    	#
    	# mutant cell line
    	#
    	addMutantCellLine(alleleKey, mutantCellLine, createdByKey)

    	# existing ikmc notes
	if len(newNotes) > 0:
        	tokens = newNotes.split('|')
		nKey = tokens[0]
		note = tokens[1] + '|' + ikmcNotes

		ikmcNoteSQL = ikmcNoteSQL + \
			'''
			update MGI_NoteChunk
			set note = '%s'
			where _Note_key = %s
			''' % (note, nKey)

	else:
	    noteFile.write('%s|%s|%s|%s|%s|%s|%s|%s\n' \
		% (noteKey, alleleKey, mgiNoteObjectKey, mgiIKMCNoteTypeKey, \
		   createdByKey, createdByKey, loaddate, loaddate))

	    noteChunkFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
		% (noteKey, 1, ikmcNotes, createdByKey, createdByKey, loaddate, loaddate))

	    noteKey = noteKey + 1

    return 1

#
# Purpose:  processes data
#
def processFile():

    global alleleKey, assocKey, refAssocKey, accKey, noteKey, mgiKey

    lineNum = 0
    # For each line in the input file

    for line in inputFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = line[:-1].split('\t')
	#print line
        try:
	    markerID = tokens[0]
	    symbol = tokens[1]
	    name = tokens[2]
	    alleleStatus = tokens[3]
	    alleleType = tokens[4]
	    alleleSubtype = tokens[5]
	    collection = tokens[6]
	    germLine = tokens[7]
	    references = tokens[8]
	    strainOfOrigin = tokens[9]
	    mutantCellLine = tokens[10]
	    molecularNotes = tokens[11]
	    driverNotes = tokens[12]
	    ikmcNotes = tokens[13]
	    mutations = tokens[14]
	    inheritanceMode = tokens[15]
	    isMixed = tokens[16]
	    isExtinct = tokens[17]
	    createdBy = tokens[18]
	    createMutantCellLine1 = tokens[19]
	    createMutantCellLine2 = tokens[20]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	# creator
	createdByKey = loadlib.verifyUser(createdBy, lineNum, errorFile)
        if createdByKey == 0:
            continue

	#if createMutantCellLine1 == 1 or len(createMutantCellLine2) > 0:
	#	processFileIKMC(createMutantCellLine1, createMutantCellLine2, \
	#		mutantCellLine, ikmcNotes, createdByKey)
	#	continue

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

	# _vocab_key = 37 (Allele Status)
	alleleStatusKey = loadlib.verifyTerm('', 37, alleleStatus, lineNum, errorFile)

	# _vocab_key = 38 (Allele Type)
	alleleTypeKey = loadlib.verifyTerm('', 38, alleleType, lineNum, errorFile)

	# _vocab_key = 61 (Allele Transmission)
	germLineKey = loadlib.verifyTerm('', 61, germLine, lineNum, errorFile)

	# _vocab_key = 36 (Allele Molecular Mutation)
	allMutations = mutations.split('|')

	# _vocab_key = 35 (Allele Status)
	inheritanceModeKey = loadlib.verifyTerm('', 35, inheritanceMode, lineNum, errorFile)

	# strains
	strainOfOriginKey = sourceloadlib.verifyStrain(strainOfOrigin, lineNum, errorFile)

	# reference
	refKey = loadlib.verifyReference(jnum, lineNum, errorFile)

        # if errors, continue to next record
	# errors are stored (via loadlib) in the .error log

        if markerKey == 0 \
		or qualifierKey == 0 \
		or markerStatusKey == 0 \
		or alleleStatusKey == 0 \
		or alleleTypeKey == 0 \
		or germLineKey == 0 \
		or allMutations == 0 \
		or inheritanceModeKey == 0 \
		or strainOfOriginKey == 0 \
		or refKey == 0 \
		or createdByKey == 0:
            continue

        # if no errors, process the allele

	# allele (master)
        alleleFile.write('%d|%s|%s|%s|%s|%s|%s|%s|%s||0|%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (alleleKey, markerKey, strainOfOriginKey, inheritanceModeKey, alleleTypeKey, \
	    alleleStatusKey, germLineKey, symbol, name, \
	    isExtinct, isMixed, \
	    createdByKey, createdByKey, createdByKey, loaddate, loaddate, loaddate))

	# allele/marker
        markerFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
	    % (assocKey, alleleKey, markerKey, qualifierKey, refKey, markerStatusKey, \
	       createdByKey, createdByKey, loaddate, loaddate))

	# molecular mutation
	for mutation in allMutations:
		mutationKey = loadlib.verifyTerm('', 36, mutation, lineNum, errorFile)
        	mutationFile.write('%s|%s|%s|%s\n' \
	    	% (alleleKey, mutationKey, loaddate, loaddate))

	allReferences = references.split('||')
	for reference in allReferences:
		refType, refID = reference.split('|')
		refKey = loadlib.verifyReference(refID, lineNum, errorFile)

		if refType == 'Original':
			refAssocTypeKey = 1011
		elif refType == 'Transmission':
			refAssocTypeKey = 1023
		elif refType == 'Molecular':
			refAssocTypeKey = 1012

        	refFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
	    		% (refAssocKey, refKey, alleleKey, mgiTypeKey, refAssocTypeKey, \
	       		createdByKey, createdByKey, loaddate, loaddate))
		refAssocKey = refAssocKey + 1

        #
        # mutant cell line
        #
        if len(mutantCellLine) > 0:
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

	# ikmc notes
	if len(ikmcNotes) > 0:

	    noteFile.write('%s|%s|%s|%s|%s|%s|%s|%s\n' \
		% (noteKey, alleleKey, mgiNoteObjectKey, mgiIKMCNoteTypeKey, \
		   createdByKey, createdByKey, loaddate, loaddate))

	    noteChunkFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
		% (noteKey, 1, ikmcNotes, createdByKey, createdByKey, loaddate, loaddate))

	    noteKey = noteKey + 1

	# Print out a new text file and attach the new MGI Allele IDs as the last field

        newAlleleFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n' \
	    % (mgi_utils.prvalue(markerID), \
	       mgi_utils.prvalue(symbol), \
	       mgi_utils.prvalue(name), \
	       mgi_utils.prvalue(alleleStatus), \
	       mgi_utils.prvalue(alleleType), \
	       mgi_utils.prvalue(alleleSubtype), \
	       mgi_utils.prvalue(collection), \
	       mgi_utils.prvalue(germLine), \
	       mgi_utils.prvalue(references), \
	       mgi_utils.prvalue(strainOfOrigin), \
	       mgi_utils.prvalue(mutantCellLine), \
	       mgi_utils.prvalue(allMutations), \
	       mgi_utils.prvalue(inheritanceMode), \
	       mgi_utils.prvalue(isMixed), \
	       mgi_utils.prvalue(isExtinct), \
	       mgi_utils.prvalue(createdBy), \
	       mgi_utils.prvalue(mgiPrefix), mgi_utils.prvalue(mgiKey)))

        accKey = accKey + 1
        mgiKey = mgiKey + 1
	assocKey = assocKey + 1
        alleleKey = alleleKey + 1

    #	end of "for line in inputFile.readlines():"

    #
    # Update the AccessionMax value
    #

    if not DEBUG:
        db.sql('exec ACC_setMax %d' % (lineNum), None)


    return 0

#
# Purpose: write 1 or more mutation cell line associations to bcp file
#
def addMutantCellLine(alleleKey, mutantCellLine, createdByKey):

    global mutantKey

    mutantCellLineKey = 0

    results = db.sql('''
	select _CellLine_key from ALL_CellLine
	where isMutant = 1 and _Derivation_key is not null
	and cellLine = '%s'
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


if initialize() != 0:
    sys.exit(1)

if setPrimaryKeys() != 0:
    sys.exit(1)

if processFile() != 0:
    sys.exit(1)

if bcpFiles() != 0:
    closeFiles()
    sys.exit(1)

closeFiles()
sys.exit(0)

