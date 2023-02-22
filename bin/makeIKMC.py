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
#      IKMC file ($IKMC_COPY_INPUT_FILE) from GenTar
#
#       field 1: Marker Symbol
#       field 2: MGI Marker ID
#       field 3: Production Plan Colony Name  
#	field 4: Production Plan  Colony Background Strain     
#	field 5: Production Work Unit  (formerly IMITs Production Centre)
#	field 6: Production Mutation Allele Symbol (full symbol now)
#	field 7: Es Cell Allele Symbol (full symbol now)
#	field 8: ES Cell Name 
#	field 9: ES Cell Allele Accession ID 
#	field 10: ES Cell Parent 
#	field 11: Modification Plan Colony Name 
#	field 12: Excision Type   
#	field 13: Tat Cre Phenotype Attempt Deleter Strain        
#	field 14: Deleter Strain
#	field 15: Modification Plan Colony Background Strain
#	field 16: Modification Work Unit 
#	field 17: Modification Mutation Accession ID (MGI Allele Accession)
#	field 18: Modification Mutation Symbol 
#       field 19: Mutation Identifier
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
#       field 13: Driver Gene (MGI_Relationship._Category_key = 1006/_Object_key_2 = marker)
#	field 14: IKMC Colony Name (_NoteType_key = 1041)
#	field 15: Molecular Mutation
#	field 16: Inheritance Mode
#	field 17: Mixed
#	field 18: Extinct
#	field 19: Created By
#	field 20: Add mutant cell line
#	field 21: Add IKMC Colony Note
#	field 22: Set the child's Allele Status = Approved (847114)
#	field 23: Allele MGI ID (if child allele already exists)
#	field 24: Allele Symbol minus Marker Symbol (for IKMC format)
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
# 10/04/2022    sc
#       - https://mgi-jira.atlassian.net/browse/WTS2-558 iMITS/GenTar (TR13439)
#
# 01/25/2022    sc
#       - wts2-767 remove mgi_notechunk, note now in mgi_note
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

logitSkip = []
logitExists = []

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

note_tmX1 = 'Cre-mediated excision of the parental %s allele resulted in the removal of the neomycin selection cassette and critical exon(s) leaving behind the inserted lacZ reporter sequence.  Further information on targeting strategies used for this and other IKMC alleles can be found at http://www.informatics.jax.org/mgihome/nomen/IKMC_schematics.shtml.'

note_tmX2 = 'Flp-mediated excision of the parental %s allele resulted in the removal of the promoter-driven neomycin selection cassette, the inserted lacZ reporter sequence, and the loxP-flanked critical exon(s). Further information on targeting strategies used for this and other IKMC alleles can be found at http://www.informatics.jax.org/mgihome/nomen/IKMC_schematics.shtml.'

note_tmXe = 'Cre-mediated excision of the parental %s allele resulted in the removal of the promoter-driven neomycin selection cassette leaving behind the inserted lacZ reporter sequence. Further information on targeting strategies used for this and other IKMC alleles can be found at http://www.informatics.jax.org/mgihome/nomen/IKMC_schematics.shtml.'

note_tmXb = 'Cre-mediated excision of the parental %s allele resulted in the removal of the promoter-driven neomycin selection cassette and critical exon(s) leaving behind the inserted lacZ reporter sequence. Further information on targeting strategies used for this and other IKMC alleles can be found at http://www.informatics.jax.org/mgihome/nomen/IKMC_schematics.shtml.'

note_tmXc = 'Flp-mediated excision of the parental %s allele resulted in the removal of the promoter-driven neomycin selection cassette and the inserted lacZ reporter sequence, leaving behind the loxP-flanked critical exon(s). Further information on targeting strategies used for this and other IKMC alleles can be found at http://www.informatics.jax.org/mgihome/nomen/IKMC_schematics.shtml.'

note_tmX2 = 'Flp-mediated excision of the parental %s allele resulted in the removal of the promoter-driven neomycin selection cassette, the inserted lacZ reporter sequence, and the loxP-flanked critical exon(s). Further information on targeting strategies used for this and other IKMC alleles can be found at http://www.informatics.jax.org/mgihome/nomen/IKMC_schematics.shtml.'

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
        print('Environment variable not set: LOG_DIAG')
        rc = 1

    #
    # Make sure the environment variables are set.
    #
    if not logCurFile:
        print('Environment variable not set: LOG_CUR')
        rc = 1

    #
    # Make sure the environment variables are set.
    #
    if not skipDiagFile:
        print('Environment variable not set: SKIP_DIAG')
        rc = 1

    #
    # Make sure the environment variables are set.
    #
    if not existsDiagFile:
        print('Environment variable not set: EXISTS_DIAG')
        rc = 1

    #
    # Make sure the environment variables are set.
    #
    if not ikmcFile:
        print('Environment variable not set: IKMC_COPY_INPUT_FILE')
        rc = 1

    # Make sure the environment variables are set.
    #
    if not alleleFile:
        print('Environment variable not set: INPUTFILE')
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
    # Parent: Allele Accession ID/Key/Symbol/Name/Strain/Marker Acc ID/Marker Symbol
    #
    # KOMP, EUCOMM only : existing parents
    # excluding NCOM (for now)
    #		or a.symbol like "%<tm%[ae](NCOM%"
    #		or a.symbol like "%<tm[0-9](NCOM%"
    #
    # includes:  Approved
    #
    print('querying for parents')
    results = db.sql('''
        select aa.accID, a._Allele_key, a._Allele_Status_key, a.symbol, a.name, a._Collection_key,
                s.strain, am.accID as markerID, m.symbol as markerSym
        from ALL_Allele a, ACC_Accession aa, ACC_Accession am, MRK_Marker m, PRB_Strain s
        where (lower(a.symbol) ~ '.*<tm([0-9])\(komp.*'
                or lower(a.symbol) ~ '.*<tm.*([ae])\(komp.*'
                or lower(a.symbol) ~ '.*<tm([0-9])\(eucomm.*'
                or lower(a.symbol) ~ '.*<tm.*([ae])\(eucomm.*'
                )
        and a._Allele_Status_key in (847114)
        and a._Allele_key = aa._Object_key
        and a._Strain_key = s._Strain_key
        and aa._MGIType_key = 11
        and aa._LogicalDB_key = 1
        and aa.prefixPart = 'MGI:'
        and aa.preferred = 1
        and a._Marker_key = am._Object_key
        and am._MGIType_key = 2
        and am.prefixPart = 'MGI:'
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
    # Child: Allele Accession ID/Key/Symbol
    #
    # KOMP, EUCOMM only : existing children
    #	excluding NCOM (for now)
    #		or a.symbol like "%<tm%.%(NCOM%"
    #
    # includes:  Approved, Reserved
    #
    print('querying for children')
    results = db.sql('''
        select aa.accID, a._Allele_key, a._Allele_Status_key, a.symbol, a._Collection_key
        from ALL_Allele a, ACC_Accession aa
        where (lower(a.symbol) ~ '.*<tm.*\..*\(komp.*'
                or lower(a.symbol) ~ '.*<tm.*b\(komp.*'
                or lower(a.symbol) ~ '.*<tm.*c\(komp.*'
                or lower(a.symbol) ~ '.*<tm.*\..*\(eucomm.*'
                or lower(a.symbol) ~ '.*<tm.*b\(eucomm.*'
                or lower(a.symbol) ~ '.*<tm.*c\(eucomm.*'
                )
        and a._Allele_Status_key in (847114, 847113)
        and a._Allele_key = aa._Object_key
        and aa._MGIType_key = 11
        and aa._LogicalDB_key = 1
        and aa.prefixPart = 'MGI:'
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
    # includes:  Approved, Reserved
    #
    print('querying for cell lines')
    results = db.sql('''
        select a._Allele_key, c._CellLine_key, c.cellLine
        from ALL_Allele a, ALL_Allele_CellLine ac, ALL_CellLine c
        where (lower(a.symbol) ~ '.*<tm.*\(komp.*'
                or lower(a.symbol) ~ '.*<tm.*\(eucomm.*'
                )
        and a._Allele_Status_key in (847114, 847113)
        and a._Allele_key = ac._Allele_key
        and ac._MutantCellLine_key = c._CellLine_key
        ''', 'auto')
    for r in results:

        # by cell line symbol
        key = r['cellLine']
        if key not in cellLineBySymbol:
                cellLineBySymbol[key] = []
        cellLineBySymbol[key].append(r)

        # by allele key
        key = r['_Allele_key']
        if key not in cellLineByKey:
                cellLineByKey[key] = []
        cellLineByKey[key].append(r)

    #
    # IKMC Notes
    #
    print('quering for ikmc notes')
    results = db.sql('''
        select n._Note_key, n._Object_key, rtrim(n.note) as note
        from MGI_Note n
        where n._NoteType_key = %s
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
        print('Cannot open diag file: ' + logDiagFile)
        return 1

    #
    # Open the Log Cur file; append to existing file
    #
    try:
        fpLogCur = open(logCurFile, 'a+')
    except:
        print('Cannot open cur file: ' + logCurFile)
        return 1

    #
    # Open the Skip Diag file
    #
    try:
        fpSkipDiag = open(skipDiagFile, 'w')
    except:
        print('Cannot open skip file: ' + skipDiagFile)
        return 1

    #
    # Open the Exists Diag file
    #
    try:
        fpExistsDiag = open(existsDiagFile, 'w')
    except:
        print('Cannot open exists file: ' + existsDiagFile)
        return 1

    #
    # Open the IKMC file
    #
    try:
        fpIKMC = open(ikmcFile, encoding='utf-8', errors='replace')
    except:
        print('Cannot open ikmc file: ' + ikmcFile)
        return 1

    #
    # Open the IKMC file with genotype sequence #
    #
    try:
        fpAllele = open(alleleFile, 'w')
    except:
        print('Cannot open allele file: ' + alleleFile)
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

    lineNum = 0
    header = 1
    print('reading input file')
    for line in fpIKMC.readlines():
        lineNum += 1
        print('line: %s' % line)
        if header == 1:
                header += 1
                continue

        error = 0
        tokens = line[:-1].split('\t')

        ikmc_marker_symbol_1 = tokens[0]
        ikmc_marker_id_2 = tokens[1]
        ikmc_allele_symbol_6 = tokens[5]
        ikmc_allele_escell_symbol_7 = tokens[6]
        ikmc_allele_id_9 = tokens[8] # 8/9 swapped in in GenTar from IMITs
        ikmc_escell_name_8 = tokens[7]
        ikmc_colony_11 = tokens[10]
        ikmc_iscre_12 = tokens[11].lower()
        ikmc_tatcre_13 = tokens[12]
        mgi_allele_id_17 = tokens[16]

        if len(mgi_allele_id_17) > 0:
                logit = 'field 17 line %s: we have already processed this row: ' % lineNum
                error = 1

        if ikmc_marker_id_2 not in markerByID:
                logit = 'field 2 line %s : marker is not in MGI: ' % lineNum
                error = 1

        if ikmc_allele_id_9 not in alleleByID:
                logit = 'field 9 line %s: allele is not in MGI or is not a tmX, tmXa, tmXe: ' % lineNum
                error = 1

        else:
                #
                # if ikmc_escell_name_8 cell line is not associated with 
                #		ikmc_allele_id_9, then skip
                #

                if ikmc_escell_name_8 not in cellLineBySymbol:
                        logit = 'field 8 line %s: es cell line is not associated with *any* allele in MGI: '% lineNum
                        error = 1
                else:
                        skipIt = 1
                        aKey = alleleByID[ikmc_allele_id_9][0]['_Allele_key']
                        cellLine = cellLineBySymbol[ikmc_escell_name_8]
                        for c in cellLine:
                                cKey = c['_Allele_key']
                                if aKey == cKey:
                                        skipIt = 0

                        if skipIt:
                                logit = 'ES Cell Name (field 9) is not associated with allele ID (field 8) line %s: ' % lineNum
                                error = 1

        if ikmc_iscre_12 not in ('cre', 'flp'):
                logit = 'Excision Type (field 12) is not "cre" or "flp" line %s: ' % lineNum
                error = 1

        if ikmc_tatcre_13 not in ('true', 'false'):
                logit = 'TAT-Cre (field 13) is not "true" or "false" line %s: ' % lineNum
                error = 1

        if error:
                logitSkip.append(logit + ikmc_marker_symbol_1 + '\t' + \
                        ikmc_marker_id_2 + '\t' + \
                        ikmc_allele_symbol_6 + '\t' + \
                        ikmc_allele_escell_symbol_7 + '\t' + \
                        ikmc_allele_id_9 + '\t' + \
                        ikmc_escell_name_8 + '\t' + \
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

        alleleSym = alleleByID[ikmc_allele_id_9][0]['symbol']
        print('alleleSym: %s' % alleleSym)
        alleleName = alleleByID[ikmc_allele_id_9][0]['name']
        alleleKey = alleleByID[ikmc_allele_id_9][0]['_Allele_key']
        collectionKey = alleleByID[ikmc_allele_id_9][0]['_Collection_key']
        strainOfOrigin = alleleByID[ikmc_allele_id_9][0]['strain']
        markerSym = alleleByID[ikmc_allele_id_9][0]['markerSym']

        #alleleSym_6 = ikmc_marker_symbol_1 + '<' + ikmc_allele_symbol_6 + '>'
        alleleSym_6 = ikmc_allele_symbol_6

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
        childKey = 0
        isReserved = 0

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
                logit = 'This tmXe allele is not Cre line %s: ' % lineNum
                logitExists.append(logit + ikmc_marker_symbol_1 + '\t' + \
                        ikmc_marker_id_2 + '\t' + \
                        ikmc_allele_symbol_6 + '\t' + \
                        ikmc_allele_escell_symbol_7 + '\t' + \
                        ikmc_allele_id_9 + '\t' + \
                        ikmc_escell_name_8 + '\t' + \
                        ikmc_iscre_12 + '\t' + \
                        ikmc_tatcre_13 + '\t' + \
                        mgi_allele_id_17 + '\t' + \
                        alleleSym + '\t' + newAlleleSym2 + '\n')
                continue

        #
        # if isXa and field 9 != field 6, then this requires some special handling
        # 

        elif isXa and alleleSym != alleleSym_6:

                # special logging requested by Kim
                #if len(ikmc_allele_symbol_6) > 4 and ikmc_allele_symbol_6[3] != "e":
                if len(ikmc_allele_symbol_6) > 4 and ikmc_allele_symbol_6.find('e(') != -1:
                        logit = "field 9 and field 6 symbols do not match line %s: " % lineNum
                else:
                        logit = 'Must handle special tmXa/tmXe case line %s: ' % lineNum

                logitSkip.append(logit + ikmc_marker_symbol_1 + '\t' + \
                        ikmc_marker_id_2 + '\t' + \
                        ikmc_allele_symbol_6 + '\t' + \
                        ikmc_allele_escell_symbol_7 + '\t' + \
                        ikmc_allele_id_9 + '\t' + \
                        ikmc_escell_name_8 + '\t' + \
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
                if (isX or isXe) and isCre and newAlleleSym1 in childAlleleBySymbol:
                        childExists = 1
                        newAlleleSym = newAlleleSym1

                elif isX and isFlp and newAlleleSym2 in childAlleleBySymbol:
                        childExists = 1
                        newAlleleSym = newAlleleSym2

                elif isXa and isCre and newAlleleSymB in childAlleleBySymbol:
                        childExists = 1
                        newAlleleSym = newAlleleSymB

                elif isXa and isFlp and newAlleleSymC in childAlleleBySymbol:
                        childExists = 1
                        newAlleleSym = newAlleleSymC

                if childExists:

                        cellLineExists = 0
                        colonyExists = 0
                        childKey = childAlleleBySymbol[newAlleleSym][0]['_Allele_key']

                        if childAlleleBySymbol[newAlleleSym][0]['_Allele_Status_key'] == 847113:
                                isReserved = 1

                        if childKey in cellLineByKey:
                                for c in cellLineByKey[childKey]:
                                        if c['cellLine'] == ikmc_escell_name_8:
                                                cellLineExists = 1

                        if childKey in ikmcNotes:
                                for c in ikmcNotes[childKey]:
                                        if c['note'].find(ikmc_colony_11) != -1:
                                                colonyExists = 1

                        #
                        # if the child exists 
                        # and the cell line exists
                        # and the colony exists
                        # and the child's status is *not* reserved
                        #
                        if cellLineExists and colonyExists and not isReserved:
                                logit = 'Child/Cell Line/Colony already exists in MGI line %s: ' % lineNum
                                logitExists.append(logit + ikmc_marker_symbol_1 + '\t' + \
                                        ikmc_marker_id_2 + '\t' + \
                                        ikmc_allele_symbol_6 + '\t' + \
                                        ikmc_allele_escell_symbol_7 + '\t' + \
                                        ikmc_allele_id_9 + '\t' + \
                                        ikmc_escell_name_8 + '\t' + \
                                        ikmc_iscre_12 + '\t' + \
                                        ikmc_tatcre_13 + '\t' + \
                                        mgi_allele_id_17 + '\t' + \
                                        alleleSym + '\t' + newAlleleSym + '\n')
                                continue
                        #else:
                        #	print alleleSym, ikmc_escell_name_8
                        # OK to go!

        #print alleleSym, ikmc_colony_11

        #
        # new Allele has passed the rules...ready to create the new allele
        #

        if isX and isCre:
                alleleType = 'Targeted'
                alleleSubType = 'Null/knockout|Reporter'
                molecularMutation = 'Insertion|Intragenic deletion'
                newAlleleSym = newAlleleSym1
                newAlleleName = newAlleleName1
                n = alleleSym.replace('>', '</sup>')
                n = n.replace('<tm', '<sup>tm')
                molecularNote = note_tmX1 % (n)

        elif isXe and isCre:
                alleleType = 'Targeted'
                alleleSubType = 'Null/knockout|Reporter'
                molecularMutation = 'Insertion'
                newAlleleSym = newAlleleSym1
                newAlleleName = newAlleleName1
                n = alleleSym.replace('>', '</sup>')
                n = n.replace('<tm', '<sup>tm')
                molecularNote = note_tmXe % (n)

        elif isX and isFlp:
                alleleType = 'Targeted'
                alleleSubType = 'Null/knockout'
                molecularMutation = 'Insertion'
                newAlleleSym = newAlleleSym2
                newAlleleName = newAlleleName2
                n = alleleSym.replace('>', '</sup>')
                n = n.replace('<tm', '<sup>tm')
                molecularNote = note_tmX2 % (n)

        elif isXa and isCre:
                alleleType = 'Targeted'
                alleleSubType = 'Null/knockout|Reporter'
                molecularMutation = 'Insertion|Intragenic deletion'
                newAlleleSym = newAlleleSymB
                newAlleleName = newAlleleNameB
                n = alleleSym.replace('>', '</sup>')
                n = n.replace('<tm', '<sup>tm')
                molecularNote = note_tmXb % (n)

        elif isXa and isFlp:
                alleleType = 'Targeted'
                alleleSubType = 'Conditional ready'
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
        attachColony = 0
        print('newAlleleSym: %s' % newAlleleSym)
        if newAlleleSym in alleleAdded:

                attachCellLine = 1
                attachColony = 1

                for a in alleleAdded[newAlleleSym]:
                        if a == ikmc_escell_name_8:
                                attachCellLine = 0

                for a in colonyAdded[newAlleleSym]:
                        if a == ikmc_colony_11:
                                attachColony = 0

                if not attachCellLine and not attachColony:
                        logit = 'Duplicate: child already added by this load line %s: ' % lineNum
                        logitExists.append(logit + ikmc_marker_symbol_1 + '\t' + \
                                ikmc_marker_id_2 + '\t' + \
                                ikmc_allele_symbol_6 + '\t' + \
                                ikmc_allele_escell_symbol_7 + '\t' + \
                                ikmc_allele_id_9 + '\t' + \
                                ikmc_escell_name_8 + '\t' + \
                                ikmc_iscre_12 + '\t' + \
                                ikmc_tatcre_13 + '\t' + \
                                mgi_allele_id_17 + '\t' + \
                                alleleSym + '\n')
                        continue
                else:
                        print('alleleAdded[%s].append(%s)' % (newAlleleSym, ikmc_escell_name_8))
                        alleleAdded[newAlleleSym].append(ikmc_escell_name_8)
                        print('colonyAdded[%s].append(%s)' % (newAlleleSym, ikmc_colony_11))
                        colonyAdded[newAlleleSym].append(ikmc_colony_11)

        # update the new allele list
        elif int(childKey) == 0:
                alleleAdded[newAlleleSym] = []
                alleleAdded[newAlleleSym].append(ikmc_escell_name_8)
                colonyAdded[newAlleleSym] = []
                colonyAdded[newAlleleSym].append(ikmc_colony_11)

        #
        # ready to create the Allele
        #

        print('ready to create the Allele')
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
        fpAllele.write(alleleSubType + '\t')

        # Allele Collection
        fpAllele.write(str(collectionKey) + '\t')

        # Transmission
        fpAllele.write('Germline' + '\t')

        # Reference
        fpAllele.write('Original|' + jnumber + '||Transmission|' + jnumber + '||Molecular|' + jnumber + '\t')

        # Strain of Origin
        fpAllele.write(strainOfOrigin + '\t')

        # Mutant Cell Line
        fpAllele.write(ikmc_escell_name_8 + '\t')

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

        #
        # Add additional mutant cell line to a new or existing allele
        #
        # 0 => use new allele key created by makeAllele.py
        # > 0 => allele/child key of existing allele
        # blank => do nothing
        #
        if attachCellLine:
                fpAllele.write('0')
        elif childExists and not cellLineExists:
                fpAllele.write(str(childKey))
        fpAllele.write('\t')

        #
        # Add IKMC Colony/Note to a new or existing allele
        #
        # child exists/ikmc note exists : update existing note
        # 	|| => _Note_key||existing colony notes
        #
        # child exists/ikmc note does not exist : add note
        # 	:: => allele/child key
        #
        # new allele/child/non-duplicate IKMC Colony
        #	0::colony(s)
        #
        # blank => do nothing
        #

        if childExists and childKey in ikmcNotes:
                ikmcNote = ikmcNotes[childKey]
                note = ikmcNote[0]['note']
                note = note.replace('\n', '')
                fpAllele.write(str(ikmcNote[0]['_Note_key']) + '||' + note)

        elif childExists and childKey not in ikmcNotes:
                fpAllele.write(str(childKey) + '::')

        elif attachColony:
                fpAllele.write('0::')
                fpAllele.write('|'.join(colonyAdded[newAlleleSym]))
        fpAllele.write('\t')

        #
        # Set the child's Allele Status = Approved
        #

        if isReserved:
                fpAllele.write(str(childKey))
        fpAllele.write('\t')

        #
        # Child Allele MGI ID
        #
        if newAlleleSym in childAlleleBySymbol:
                fpAllele.write(childAlleleBySymbol[newAlleleSym][0]['accID'])
        fpAllele.write('\t')

        #
        # Allele Symbol nomenclature minus the Marker name
        # 	Syt17<tm1b(KOMP)Wtsi> => tm1b(KOMP)Wtsi
        #

        p1 = newAlleleSym.find('<')
        p2 = newAlleleSym.find('>')
        fpAllele.write(newAlleleSym[p1+1:p2] + '\n')

    return 0

def writeReports():
    logitSkip.sort()
    fpSkipDiag.write('\n'.join(logitSkip))
    fpExistsDiag.write('\n'.join(logitExists))
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

writeReports()
closeFiles()
sys.exit(0)
