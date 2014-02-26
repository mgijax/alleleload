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
#      see ikmc.config
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
#	field 20: Add mutant cell line/IKMC to new Allele (not yet in database)
#	field 21: Add mutant cell line/IKMC to existing Allele (in database)
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
#      3) Morph the IKMC input file into a general-Allele input file
#      4) Close files.
#
#  01/27/2014	lec
#	- TR11515/IKMC/allele derivation load
#
###########################################################################

import sys 
import os
import db

# LOG_DIAG
# LOG_CUR
# SKIP_DIAG
# EXISTS_DIAG
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

alleleByID = {}
childAlleleBySymbol = {}
markerByID = []
cellLineBySymbol = {}
cellLineByKey = {}
alleleAdded = {}
colonyAdded = {}
ikmcNotes = {}

jnumber = ''
createdBy = ''
mgiIKMCNoteTypeKey = 1041

header = 'error\tfield 1\tfield 2\tfield 6\tfield 7\tfield 8\tfield 9\tfield 12\tfield 13\tfield 17\tfull allele symbol\tnew allele symbol\n'

note_tmX1 = 'Cre-mediated excision of the parental %s allele resulted in the removal of the neomycin selection cassette and critical exon(s) leaving behind the inserted lacZ reporter sequence.  Further information on targeting strategies used for this and other KOMP alleles can be found at http://www.knockoutmouse.org/aboutkompstrategies.'

note_tmX2 = 'Flp-mediated excision of the parental %s allele resulted in the removal of the promoter-driven neomycin selection cassette, the inserted lacZ reporter sequence, and the loxP-flanked critical exon(s). Further information on targeting strategies used for this and other KOMP alleles can be found at http://www.knockoutmouse.org/aboutkompstrategies.'

note_tmXe = 'Cre-mediated excision of the parental %s allele resulted in the removal of the promoter-driven neomycin selection cassette leaving behind the inserted lacZ reporter sequence. Further information on targeting strategies used for this and other KOMP alleles can be found at http://www.knockoutmouse.org/aboutkompstrategies.'

note_tmXb = 'Cre-mediated excision of the parental %s allele resulted in the removal of the promoter-driven neomycin selection cassette and critical exon(s) leaving behind the inserted lacZ reporter sequence. Further information on targeting strategies used for this and other KOMP alleles can be found at http://www.knockoutmouse.org/aboutkompstrategies.'

note_tmXc = 'Flp-mediated excision of the parental %s allele resulted in the removal of the promoter-driven neomycin selection cassette and the inserted lacZ reporter sequence, leaving behind the loxP-flanked critical exon(s). Further information on targeting strategies used for this and other KOMP alleles can be found at http://www.knockoutmouse.org/aboutkompstrategies.'

note_tmX2 = 'Flp-mediated excision of the parental %s allele resulted in the removal of the promoter-driven neomycin selection cassette, the inserted lacZ reporter sequence, and the loxP-flanked critical exon(s). Further information on targeting strategies used for this and other KOMP alleles can be found at http://www.knockoutmouse.org/aboutkompstrategies.'

#
# Purpose: Initialization
#
def initialize():
    global logDiagFile, logCurFile, skipDiagFile, existsDiagFile
    global ikmcFile, alleleFile
    global fpLogDiag, fpLogCur, fpSkipDiag, fpExistsDiag
    global fpIKMC, fpAllele
    global alleleByID, childAlleleBySymbol, markerByID
    global cellLineBySymbol, cellLineByKey
    global ikmcNotes
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
    if not logDiagFile:
        print 'Environment variable not set: LOG_DIAG'
        rc = 1

    #
    # Make sure the environment variables are set.
    #
    if not logCurFile:
        print 'Environment variable not set: LOG_CUR'
        rc = 1

    #
    # Make sure the environment variables are set.
    #
    if not skipDiagFile:
        print 'Environment variable not set: SKIP_DIAG'
        rc = 1

    #
    # Make sure the environment variables are set.
    #
    if not existsDiagFile:
        print 'Environment variable not set: EXISTS_DIAG'
        rc = 1

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
    # Allele Accession ID/Key/Symbol/Name/Strain/Marker Acc ID/Marker Symbol
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
	alleleByID[key] = []
	alleleByID[key].append(r)
	key = r['markerID']
	markerByID.append(key)

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
		or a.symbol like "%<tm%b(KOMP%"
		or a.symbol like "%<tm%c(KOMP%"
		or a.symbol like "%<tm%.%(EUCOMM%"
		or a.symbol like "%<tm%b(EUCOMM%"
		or a.symbol like "%<tm%c(EUCOMM%"
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
	childAlleleBySymbol[key] = []
	childAlleleBySymbol[key].append(r)

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

	# by cell line symbol
	key = r['cellLine']
	if not cellLineBySymbol.has_key(key):
		cellLineBySymbol[key] = []
	cellLineBySymbol[key].append(r)

	# by allele key
	key = r['_Allele_key']
	if not cellLineByKey.has_key(key):
		cellLineByKey[key] = []
	cellLineByKey[key].append(r)

    #
    # IKMC Notes
    #
    print 'quering for ikmc notes'
    results = db.sql('''
	select n._Note_key, n._Object_key, rtrim(c.note) as note
	from MGI_Note n, MGI_NoteChunk c
	where n._NoteType_key = %s
	and n._Note_key = c._Note_key
	''' % (mgiIKMCNoteTypeKey), 'auto')
    for r in results:
	key = r['_Object_key']
	ikmcNotes[key] = []
	ikmcNotes[key].append(r)

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
	ikmc_iscre_12 = tokens[11].lower()
	ikmc_tatcre_13 = tokens[12]
	mgi_allele_id_17 = tokens[16]

	if len(mgi_allele_id_17) > 0:
		logit = 'field 17: we have already processed this row\n'
		error = 1

	if ikmc_marker_id_2 not in markerByID:
		logit = 'field 2 : marker is not in MGI\n'
		error = 1

	if not alleleByID.has_key(ikmc_allele_id_8):
		logit = 'field 8: allele is not in MGI or is not a tmX, tmXa, tmXe\n'
		error = 1

	else:
		#
		# if ikmc_escell_name_9 cell line is not associated with 
		#		ikmc_allele_id_8, then skip
		#

		if not cellLineBySymbol.has_key(ikmc_escell_name_9):
			logit = 'field 9: es cell line is not associated with *any* allele in MGI\n'
			error = 1
		else:
			skipIt = 1
			aKey = alleleByID[ikmc_allele_id_8][0]['_Allele_key']
			cellLine = cellLineBySymbol[ikmc_escell_name_9]
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

	# Allele will be 'tmXa', 'tmX', 'tmXe'
	isXa = 0		# tm1a, tm2a, etc.
	isX = 0			# tm1, tm2, etc.
	isXe = 0		# tm1e, tm2e, etc.
	isCre = 0
	isFlp = 0

	alleleSym = alleleByID[ikmc_allele_id_8][0]['symbol']
	alleleName = alleleByID[ikmc_allele_id_8][0]['name']
	alleleKey = alleleByID[ikmc_allele_id_8][0]['_Allele_key']
	strainOfOrigin = alleleByID[ikmc_allele_id_8][0]['strain']
	markerSym = alleleByID[ikmc_allele_id_8][0]['markerSym']
	alleleSym_6 = ikmc_marker_symbol_1 + '<' + ikmc_allele_symbol_6 + '>'

	tokens1 = alleleSym.split('<')
	tokens2 = tokens1[1].split('(')
	newAlleleSym1 = alleleSym.replace(tokens2[0], tokens2[0] + '.1')
	newAlleleSym2 = alleleSym.replace(tokens2[0], tokens2[0] + '.2')
	newAlleleSymB = alleleSym.replace('a(', 'b(')
	newAlleleSymC = alleleSym.replace('a(', 'c(')

	tokens3 = alleleName.split('targeted mutation')
	tokens4 = tokens3[1].split(',')
	newAlleleName1 = alleleName.replace(tokens4[0], tokens4[0] + '.1')
	newAlleleName2 = alleleName.replace(tokens4[0], tokens4[0] + '.2')
	newAlleleNameB = alleleName.replace('a,', 'b,')
	newAlleleNameC = alleleName.replace('a,', 'c,')

	molecularNote = ''

	childExists = 0
	cellLineExists = 0
	childKey = ''

	# determine isXa, isXe, isX

	if alleleSym.find('a(') != -1:
		isXa = 1
	elif alleleSym.find('e(') != -1:
		isXe = 1
	else:
		isX = 1

	# determine cre/flp

	if ikmc_iscre_12 == 'cre':
		isCre = 1
	else:
		isFlp = 1

	#
	# if tmXe is not Cre
	#
	if isXe and not isCre:
		logit = 'This tmXe allele is not Cre.\n'
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

	#
	# if isXa and field 8 != field 6, then this requires some special handling
	# 

	elif isXa and alleleSym != alleleSym_6:

		# special logging requested by Kim
		#if len(ikmc_allele_symbol_6) > 4 and ikmc_allele_symbol_6[3] != "e":
		if len(ikmc_allele_symbol_6) > 4 and ikmc_allele_symbol_6.find('e(') != -1:
			logit = "field 8 and field 6 symbols do not match"
		else:
			logit = 'Must handle special tmXa/tmXe case\n'

		fpSkipDiag.write(logit + '\t' + \
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

	#
	# if the child already exist:
	#	if mutant cell line is not attached to child:
	# 		add additional mutant cell line and IMKC colony name
	#

	else:
		if (isX or isXe) and isCre and childAlleleBySymbol.has_key(newAlleleSym1):
			childExists = 1
			newAlleleSym = newAlleleSym1

		elif isX and isFlp and childAlleleBySymbol.has_key(newAlleleSym2):
			childExists = 1
			newAlleleSym = newAlleleSym2

		elif isXa and isCre and childAlleleBySymbol.has_key(newAlleleSymB):
			childExists = 1
			newAlleleSym = newAlleleSymB

		elif isXa and isFlp and childAlleleBySymbol.has_key(newAlleleSymC):
			childExists = 1
			newAlleleSym = newAlleleSymC

		if childExists:
			cellLineExists = 0
			childKey = childAlleleBySymbol[newAlleleSym][0]['_Allele_key']

			if cellLineByKey.has_key(childKey):
				for c in cellLineByKey[childKey]:
					if c['cellLine'] == ikmc_escell_name_9:
						cellLineExists = 1

			if cellLineExists:
				logit = 'Child & Cell Line already exists in MGI'
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
					alleleSym + '\t' + newAlleleSym + '\n')
				continue
			#else:
			#	print alleleSym, ikmc_escell_name_9
			# OK to go!

	#
	# new Allele has passed the rules...ready to create the new allele
	#

	if isX and isCre:
		alleleType = 'Targeted (Reporter)'
		molecularMutation = 'Insertion|Intragenic deletion'
		newAlleleSym = newAlleleSym1
		newAlleleName = newAlleleName1
		n = alleleSym.replace('>', '</sup>')
		n = n.replace('<tm', '<sup>tm')
		molecularNote = note_tmX1 % (n)

	elif isXe and isCre:
		alleleType = 'Targeted (Reporter)'
		molecularMutation = 'Insertion'
		newAlleleSym = newAlleleSym1
		newAlleleName = newAlleleName1
		n = alleleSym.replace('>', '</sup>')
		n = n.replace('<tm', '<sup>tm')
		molecularNote = note_tmXe % (n)

	elif isX and isFlp:
		alleleType = 'Targeted (knock-out)'
		molecularMutation = 'Insertion'
		newAlleleSym = newAlleleSym2
		newAlleleName = newAlleleName2
		n = alleleSym.replace('>', '</sup>')
		n = n.replace('<tm', '<sup>tm')
		molecularNote = note_tmX2 % (n)

	elif isXa and isCre:
		alleleType = 'Targeted (Reporter)'
		molecularMutation = 'Insertion|Intragenic deletion'
		newAlleleSym = newAlleleSymB
		newAlleleName = newAlleleNameB
		n = alleleSym.replace('>', '</sup>')
		n = n.replace('<tm', '<sup>tm')
		molecularNote = note_tmXb % (n)

	elif isXa and isFlp:
		alleleType = 'Targeted (Floxed/Frt)'
		molecularMutation = 'Insertion'
		newAlleleSym = newAlleleSymC
		newAlleleName = newAlleleNameC
		n = alleleSym.replace('>', '</sup>')
		n = n.replace('<tm', '<sup>tm')
		molecularNote = note_tmXc % (n)
		
	#
	# if the new Allele has already been created (it's a duplicate)
	#

	attachCellLine = 0
	duplicateAllele = 0

	if alleleAdded.has_key(newAlleleSym):

		attachCellLine = 1
		duplicateAllele = 0

		for a in alleleAdded[newAlleleSym]:
			if a == ikmc_escell_name_9:
				attachCellLine = 0
				duplicateAllele = 1

		if duplicateAllele:
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
				alleleSym + '\n')
			continue
		else:
			alleleAdded[newAlleleSym].append(ikmc_escell_name_9)
			colonyAdded[newAlleleSym].append(ikmc_colony_11)

	# update the new allele list
	else:
		alleleAdded[newAlleleSym] = []
		alleleAdded[newAlleleSym].append(ikmc_escell_name_9)
		colonyAdded[newAlleleSym] = []
		colonyAdded[newAlleleSym].append(ikmc_colony_11)

	#
	# ready to create the Allele
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
	alleleType = alleleType.replace('Targeted (Reporter)', 'Targeted')
	fpAllele.write(alleleType + '\t')

	# Allele Subtype
	fpAllele.write('\t')

	# Allele Collection
	fpAllele.write('Not Specified\t')

	# Transmission
	fpAllele.write('Germline' + '\t')

	# Reference
	fpAllele.write('Original|' + jnumber + '||Transmission|' + jnumber + '||Molecular|' + jnumber + '\t')

	# Strain of Origin
	fpAllele.write(strainOfOrigin + '\t')

	# Mutant Cell Line
	fpAllele.write(ikmc_escell_name_9 + '\t')

	# Molecular Notes
	fpAllele.write(molecularNote + '\t')

	# Drive Note
	fpAllele.write('\t')

	# IKMC Allele Colony Name Note (1041)
	fpAllele.write(ikmc_colony_11 + '\t')

	# Molecular Mutation
	fpAllele.write(molecularMutation + '\t')

	# Inheritance Mode
	fpAllele.write('Not Applicable' + '\t')

	# Mixed
	fpAllele.write('0' + '\t')

	# Exitinct
	fpAllele.write('0' + '\t')

	# Created By
	fpAllele.write(createdBy + '\t')

	# Add additional mutant cell line to a new allele in the same input file
	if attachCellLine:
		if colonyAdded.has_key(newAlleleSym):
			fpAllele.write('|'.join(colonyAdded[newAlleleSym]) + '\t')
	else:
		fpAllele.write('\t')

	# Add additional mutant cell line to an existing allele that is already in the database
	if childExists and not cellLineExists:
		fpAllele.write(str(childKey) + '::')
		if ikmcNotes.has_key(childKey):
			ikmcNote = ikmcNotes[childKey]
			fpAllele.write(str(ikmcNote[0]['_Note_key']) + '||' + ikmcNote[0]['note'])
		fpAllele.write('\n')
	else:
		fpAllele.write('\n')

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

