#!/usr/local/bin/python
#
#  makeIKMC.py
###########################################################################
#
#  Purpose:
#
#      This script will use the records in the IKMC input file to:
#
#      1) create an Allele input file
#
#  Usage:
#
#      makeIKMC.py
#
#  Env Vars:
#
#      The following environment variables are set by the configuration
#      file that is sourced by the wrapper script:
#
#	   IKMC_COPY_INPUT_FILE
#    	   INPUTFILE
#
#  Inputs:
#
#      IKMC file ($IKMC_COPY_INPUT_FILE)
#
#       field 1: Marker Symbol
#       field 2: MGI Marker ID
#       field 3: Mi Attempt Colony Name  
#	field 4: Mi Attempt Colony Background Strain     
#	field 5: Mi Attempt Production Centre    
#	field 6: Mi Attempt Allele Symbol        
#	field 7: Mi Attempt Es Cell Allele Symbol        
#	field 8: Mi Attempt Es Cell MGI Allele Accession 
#	field 9: Mi Attempt Es Cell Name 
#	field 10: Mi Attempt Es Cell Line 
#	field 11: Colony Name     
#	field 12: Excision Type   
#	field 13: Tat Cre Phenotype Attempt Deleter Strain        
#	field 14: Phenotype Attempt Deleter Strain
#	field 15: Phenotype Attempt Colony Background StraiN
#	field 16: Phenotype Attempt Production Centre     
#	field 17: MGI Allele Accession    
#	field 18: MGI Allele Name
#
#  Outputs:
#
#	Allele file ($INPUTFILE):
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
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Initialize variables.
#      2) Open files.
#      3) Morph the IKMCme input file into a general-Allele input file
#      4) Close files.
#
#  01/27/2014	lec
#	- TR11515/IKMC/allele derivation load
#
###########################################################################

import sys 
import os
import db
from sets import Set

# LOG_DIAG
# LOG_CUR
# OUTPUTDIR
logDiagFile = None
logCurFile = None
skipDiagFile = None
existsDiagFile = None

# IKMC_COPY_INPUT_FILE
ikmcFile = None

# INPUTFILE
alleleFile = None

# file pointers
fpLogDiag = None
fpLogCur = None
fpSkipDiag = None
fpExistsDiag = None
fpIKMC = None
fpAllele = None

alleleLookup = {}
childAlleleLookup = {}
markerLookup = []
cellLineLookup = {}
allelesAdded = []

jnumber = ''
createdBy = ''

header = 'error\tfield 1\tfield 2\tfield 6\tfield 7\tfield 8\tfield 9\tfield 12\tfield 13\tfield 17\tfull allele symbol\tnew allele symbol\n'

#
# Purpose: Initialization
#
def initialize():
    global logDiagFile, logCurFile, skipDiagFile, existsDiagFile
    global ikmcFile, alleleFile
    global fpLogDiag, fpLogCur, fpSkipDiag, fpExistsDiag
    global fpIKMC, fpAllele
    global alleleLookup, childAlleleLookup, markerLookup, cellLineLookup
    global jnumber, createdBy

    logDiagFile = os.getenv('LOG_DIAG')
    logCurFile = os.getenv('LOG_CUR')
    skipDiagFile = os.getenv('SKIP_DIAG')
    existsDiagFile = os.getenv('EXISTS_DIAG')
    ikmcFile = os.getenv('IKMC_COPY_INPUT_FILE')
    alleleFile = os.getenv('INPUTFILE')
    jnumber = os.getenv('JNUMBER')
    createdBy = os.getenv('CREATEDBY')

    rc = 0

    #
    # Make sure the environment variables are set.
    #
    if not ikmcFile:
        print 'Environment variable not set: IKMC_COPY_INPUT_FILE'
        rc = 1

    # Make sure the environment variables are set.
    #
    if not alleleFile:
        print 'Environment variable not set: INPUTFILE'
        rc = 1

    #
    # Initialize file pointers.
    #
    fpLogDiag = None
    fpLogCur = None
    fpSkipDiag = None
    fpExistsDiag = None
    fpIKMC = None
    fpAllele = None

    #
    # Allele Accession ID/Key/Symbol
    #
    # KOMP, EUCOMM only : existing parents
    # excluding NCOM (for now)
    #		or a.symbol like "%<tm%[ae](NCOM%"
    #		or a.symbol like "%<tm[0-9](NCOM%"
    #
    print 'querying for parents'
    results = db.sql('''
	select aa.accID, a._Allele_key, a.symbol, a.name, s.strain,
		am.accID as markerID, m.symbol as markerSym
	from ALL_Allele a, ACC_Accession aa, ACC_Accession am, MRK_Marker m, PRB_Strain s
	where (a.symbol like "%<tm[ae](KOMP%"
	        or a.symbol like "%<tm[0-9](KOMP%"
		or a.symbol like "%<tm%[ae](KOMP%"
		or a.symbol like "%<tm[0-9](EUCOMM%"
		or a.symbol like "%<tm%[ae](EUCOMM%"
		)
	and a._Allele_Status_key in (847114, 3983021)
	and a._Allele_key = aa._Object_key
	and a._Strain_key = s._Strain_key
	and aa._MGIType_key = 11
	and aa._LogicalDB_key = 1
	and aa.prefixPart = "MGI:"
	and aa.preferred = 1
	and a._Marker_key = am._Object_key
	and am._MGIType_key = 2
	and am.prefixPart = "MGI:"
	and am._LogicalDB_key = 1
	and am.preferred = 1
	and a._Marker_key = m._Marker_key
	''', 'auto')
    for r in results:
	key = r['accID']
	alleleLookup[key] = []
	alleleLookup[key].append(r)
	key = r['markerID']
	markerLookup.append(key)

    #
    # Allele Accession ID/Key/Symbol
    #
    # KOMP, EUCOMM only : existing children
    #	excluding NCOM (for now)
    #		or a.symbol like "%<tm%.%(NCOM%"
    #
    print 'querying for children'
    results = db.sql('''
	select aa.accID, a._Allele_key, a.symbol
	from ALL_Allele a, ACC_Accession aa
	where (a.symbol like "%<tm%.%(KOMP%"
		or a.symbol like "%<tm%.%(EUCOMM%"
		)
	and a._Allele_Status_key in (847114, 3983021)
	and a._Allele_key = aa._Object_key
	and aa._MGIType_key = 11
	and aa._LogicalDB_key = 1
	and aa.prefixPart = "MGI:"
	and aa.preferred = 1
	''', 'auto')
    for r in results:
	key = r['symbol']
	childAlleleLookup[key] = []
	childAlleleLookup[key].append(r)

    #
    # Mutant Cell Lines and their Alleles
    #
    # KOMP, EUCOMM only
    #	excluding NCOM (for now)
    #		or a.symbol like "%<tm%(NCOM%"
    #
    print 'querying for cell lines'
    results = db.sql('''
	select a._Allele_key, c._CellLine_key, c.cellLine
	from ALL_Allele a, ALL_Allele_CellLine ac, ALL_CellLine c
	where (a.symbol like "%<tm%(KOMP%"
		or a.symbol like "%<tm%(EUCOMM%"
		)
	and a._Allele_Status_key in (847114, 3983021)
	and a._Allele_key = ac._Allele_key
	and ac._MutantCellLine_key = c._CellLine_key
	''', 'auto')
    for r in results:
	key = r['cellLine']
	if not cellLineLookup.has_key(key):
		cellLineLookup[key] = []
	cellLineLookup[key].append(r)

    return rc


#
# Purpose: Open files.
#
def openFiles():
    global fpLogDiag, fpLogCur, fpSkipDiag, fpExistsDiag
    global fpIKMC, fpAllele

    #
    # Open the Log Diag file; append to existing file
    #
    try:
        fpLogDiag = open(logDiagFile, 'a+')
    except:
        print 'Cannot open file: ' + logDiagFile
        return 1

    #
    # Open the Log Cur file; append to existing file
    #
    try:
        fpLogCur = open(logCurFile, 'a+')
    except:
        print 'Cannot open file: ' + logCurFile
        return 1

    #
    # Open the Skip Diag file
    #
    try:
        fpSkipDiag = open(skipDiagFile, 'w')
    except:
        print 'Cannot open file: ' + skipDiagFile
        return 1

    #
    # Open the Exists Diag file
    #
    try:
        fpExistsDiag = open(existsDiagFile, 'w')
    except:
        print 'Cannot open file: ' + existsDiagFile
        return 1

    #
    # Open the IKMC/Biomart file
    #
    try:
        fpIKMC = open(ikmcFile, 'r')
    except:
        print 'Cannot open file: ' + ikmcFile
        return 1

    #
    # Open the IKMC file with genotype sequence #
    #
    try:
        fpAllele = open(alleleFile, 'w')
    except:
        print 'Cannot open file: ' + alleleFile
        return 1


    fpSkipDiag.write(header)
    fpExistsDiag.write(header)

    return 0


#
# Purpose: Close files.
#
def closeFiles():

    if fpLogDiag:
        fpLogDiag.close()

    if fpLogCur:
        fpLogCur.close()

    if fpSkipDiag:
        fpSkipDiag.close()

    if fpExistsDiag:
        fpExistsDiag.close()

    if fpIKMC:
        fpIKMC.close()

    if fpAllele:
        fpAllele.close()

    return 0


#
# Purpose: Read the IKMC file and re-format it to create a general-Allele input file
#
def createAlleleFile():

    lineNum = 1

    print 'reading input file'
    for line in fpIKMC.readlines():

	if lineNum == 1:
		lineNum += 1
		continue

	error = 0
        tokens = line[:-1].split('\t')

	ikmc_marker_symbol_1 = tokens[0]
        ikmc_marker_id_2 = tokens[1]
	ikmc_allele_symbol_6 = tokens[5]
	ikmc_allele_escell_symbol_7 = tokens[6]
	ikmc_allele_id_8 = tokens[7]
	ikmc_escell_name_9 = tokens[8]
	ikmc_colony_11 = tokens[10]
	ikmc_iscre_12 = tokens[11]
	ikmc_tatcre_13 = tokens[12]
	mgi_allele_id_17 = tokens[16]

	symbol = 'create this'
	name = 'create this'
	status = 'Autoload'
	inheritance = ''
	mixed = ''
	extinct = ''
	creator = ''

	if len(mgi_allele_id_17) > 0:
		logit = 'field 17: we have already processed this row\n'
		error = 1

	if ikmc_marker_id_2 not in markerLookup:
		logit = 'field 2 : marker is not in MGI\n'
		error = 1

	if not alleleLookup.has_key(ikmc_allele_id_8):
		logit = 'field 8: allele is not in MGI and is not a tmX, tmXa, tmXe\n'
		error = 1

	else:
		#
		# if ikmc_escell_name_9 cell line is not associated with 
		#		ikmc_allele_id_8, then skip
		#

		if not cellLineLookup.has_key(ikmc_escell_name_9):
			logit = 'field 9: es cell line is not associated with *any* allele in MGI\n'
			error = 1
		else:
			skipIt = 1
			aKey = alleleLookup[ikmc_allele_id_8][0]['_Allele_key']
			cellLine = cellLineLookup[ikmc_escell_name_9]
			for c in cellLine:
				cKey = c['_Allele_key']
				if aKey == cKey:
					skipIt = 0

			if skipIt:
				logit = 'ES Cell Name (field 9) is not associated with allele ID (field 8)\n'
				error = 1

	if ikmc_iscre_12 not in ('cre', 'flp'):
		logit = 'Excision Type (field 12) is not "cre" or "flp"\n'
		error = 1

	if ikmc_tatcre_13 not in ('t', 'f'):
		logit = 'TAT-Cre (field 13) is not "t" or "f"\n'
		error = 1

	if error:
		fpSkipDiag.write(logit + '\t' + \
			ikmc_marker_symbol_1 + '\t' + \
			ikmc_marker_id_2 + '\t' + \
		 	ikmc_allele_symbol_6 + '\t' + \
		 	ikmc_allele_escell_symbol_7 + '\t' + \
		 	ikmc_allele_id_8 + '\t' + \
		 	ikmc_escell_name_9 + '\t' + \
		 	ikmc_iscre_12 + '\t' + \
		 	ikmc_tatcre_13 + '\t' + \
		 	mgi_allele_id_17 + '\n')
		continue

	# Allele Symbol # Allele will be 'tmXa', 'tmX', 'tmXe'
	isA = 0		# tm1a, tm2a, etc.
	isX = 0		# tm1, tm2, etc.
	isE = 0		# tm1e, tm2e, etc.
	isCre = 0
	isFlp = 0

	alleleSym = alleleLookup[ikmc_allele_id_8][0]['symbol']
	alleleName = alleleLookup[ikmc_allele_id_8][0]['name']
	strainOfOrigin = alleleLookup[ikmc_allele_id_8][0]['strain']
	markerSym = alleleLookup[ikmc_allele_id_8][0]['markerSym']

	tokens1 = alleleSym.split('<')
	tokens2 = tokens1[1].split('(')
	newAlleleSym1 = alleleSym.replace(tokens2[0], tokens2[0] + '.1')
	newAlleleSym2 = alleleSym.replace(tokens2[0], tokens2[0] + '.2')

	tokens3 = alleleName.split('targeted mutation')
	tokens4 = tokens3[1].split(',')
	newAlleleName1 = alleleName.replace(tokens4[0], tokens4[0] + '.1')
	newAlleleName2 = alleleName.replace(tokens4[0], tokens4[0] + '.2')

	if alleleSym.find('a(') != -1:
		isA = 1
	elif alleleSym.find('e(') != -1:
		isE = 1
	else:
		isX = 1

	if ikmc_iscre_12 == 'cre':
		isCre = 1
	else:
		isFlp = 1

	# skipping for now...
	if not isX:
		continue

	#
	# if child:1 already exists
	#

	if isX and isCre and childAlleleLookup.has_key(newAlleleSym1):
		logit = 'Child already exists in MGI as Cre/tmX.1'
		fpExistsDiag.write(logit + '\t' + \
			ikmc_marker_symbol_1 + '\t' + \
			ikmc_marker_id_2 + '\t' + \
		 	ikmc_allele_symbol_6 + '\t' + \
		 	ikmc_allele_escell_symbol_7 + '\t' + \
		 	ikmc_allele_id_8 + '\t' + \
		 	ikmc_escell_name_9 + '\t' + \
		 	ikmc_iscre_12 + '\t' + \
		 	ikmc_tatcre_13 + '\t' + \
		 	mgi_allele_id_17 + '\t' + \
			alleleSym + '\t' + newAlleleSym1 + '\n')
		continue

	#
	# if child:2 already exists
	#

	elif isX and isFlp and childAlleleLookup.has_key(newAlleleSym2):
		logit = 'Child already exists in MGI as Flp/tmX.2\n'
		fpExistsDiag.write(logit + '\t' + \
			ikmc_marker_symbol_1 + '\t' + \
			ikmc_marker_id_2 + '\t' + \
		 	ikmc_allele_symbol_6 + '\t' + \
		 	ikmc_allele_escell_symbol_7 + '\t' + \
		 	ikmc_allele_id_8 + '\t' + \
		 	ikmc_escell_name_9 + '\t' + \
		 	ikmc_iscre_12 + '\t' + \
		 	ikmc_tatcre_13 + '\t' + \
		 	mgi_allele_id_17 + '\t' + \
			alleleSym + '\t' + newAlleleSym2 + '\n')
		continue

	if isX and isCre:
		alleleType = 'Targeted (Reporter)'
		molecularMutation = 'Insertion|Intragenic deletion'
		newAlleleSym = newAlleleSym1
		newAlleleName = newAlleleName1
	elif isX and isFlp:
		alleleType = 'Targeted (knock-out)'
		molecularMutation = 'Insertion'
		newAlleleSym = newAlleleSym2
		newAlleleName = newAlleleName2

	if newAlleleSym in allelesAdded:
		logit = 'Duplicate: child already added by this load\n'
		fpExistsDiag.write(logit + '\t' + \
			ikmc_marker_symbol_1 + '\t' + \
			ikmc_marker_id_2 + '\t' + \
		 	ikmc_allele_symbol_6 + '\t' + \
		 	ikmc_allele_escell_symbol_7 + '\t' + \
		 	ikmc_allele_id_8 + '\t' + \
		 	ikmc_escell_name_9 + '\t' + \
		 	ikmc_iscre_12 + '\t' + \
		 	ikmc_tatcre_13 + '\t' + \
		 	mgi_allele_id_17 + '\t' + \
			alleleSym + '\t' + newAlleleSym2 + '\n')
		continue

	allelesAdded.append(newAlleleSym)

	#
	# if we made it this far, then we can create an Allele input row
	#

	# Marker ID
	fpAllele.write(ikmc_marker_id_2 + '\t')

	# Allele Symbol
	fpAllele.write(newAlleleSym + '\t')

	# Allele Name
	fpAllele.write(newAlleleName + '\t')

	# Allele Status
	fpAllele.write('Approved' + '\t')

	# Allele Type
	fpAllele.write(alleleType + '\t')

	# Allele Subtype
	fpAllele.write('\t')

	# Allele Collection
	fpAllele.write('\t')

	# Transmission
	fpAllele.write('Germline' + '\t')

	# Reference
	fpAllele.write('Original|' + jnumber + '\t')

	# Strain of Origin
	fpAllele.write(strainOfOrigin + '\t')

	# Mutant Cell Line
	fpAllele.write(ikmc_escell_name_9 + '\t')

	# Molecular Notes
	fpAllele.write('\t')

	# Drive Note
	fpAllele.write('\t')

	# IKMC Allele Colony Name Note (1041)
	fpAllele.write(ikmc_colony_11 + '\t')

	# Molecular Mutation
	fpAllele.write(molecularMutation + '\t')

	# Inheritance Mode
	# Mixed
	# Exitinct
	# Created By
	fpAllele.write('Not Applicable' + '\t')
	fpAllele.write('0' + '\t')
	fpAllele.write('0' + '\t')
	fpAllele.write(createdBy + '\n')

	lineNum += 1

    return 0

#
#  MAIN
#

if initialize() != 0:
    sys.exit(1)

if openFiles() != 0:
    sys.exit(1)

if createAlleleFile() != 0:
    closeFiles()
    sys.exit(1)

closeFiles()
sys.exit(0)

