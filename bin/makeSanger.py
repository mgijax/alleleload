#!/usr/local/bin/python
#
#  makeSanger.py
###########################################################################
#
#  Purpose:
#
#      This script will use the records in the Sanger input file to:
#
#      1) create an Allele input file
#
#  Usage:
#
#      makeSanger.py
#
#  Env Vars:
#
#      The following environment variables are set by the configuration
#      file that is sourced by the wrapper script:
#
#	   SANGER_COPY_INPUT_FILE
#    	   INPUTFILE
#
#  Inputs:
#
#      Sanger file ($SANGER_COPY_INPUT_FILE)
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
#       field 1:  MGI Marker ID
#       field 2:  Allele Symbol
#       field 3:  Allele Name
#       field 4:  Allele Status
#       field 5:  Allele Generate-Type
#       field 6:  Allele Attribute/SubType
#       field 7:  Germ Line Transmission
#       field 8:  Reference Type/J#
#       field 9:  Strain of Origin
#       field 10  Mutant Cell Line ID
#       field 11: Molecular Notes
#       field 12: Driver Notes
#       field 13: Molecular Mutation
#       field 14: Inheritance
#       field 15: Mixed
#       field 16: Extinct
#       field 17: Creation Date
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#
#  Assumes:  Nothing
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Initialize variables.
#      2) Open files.
#      3) Morph the Sangerme input file into a Hight-Throughput MP input file
#      4) Close files.
#
#  Notes:  None
#
#  01/27/2014	lec
#	- TR11515/Sanger/allele derivation load
#
###########################################################################

import sys 
import os

# SANGER_COPY_INPUT_FILE
sangerFile = None

# INPUTFILE
alleleFile = None

# file pointers
fpSanger = None
fpAllele = None

#
# Purpose: Initialization
# Returns: 1 if file does not exist or is not readable, else 0
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def initialize():
    global sangerFile, alleleFile
    global fpSanger, fpAllele

    sangerFile = os.getenv('SANGER_COPY_INPUT_FILE')
    alleleFile = os.getenv('INPUTFILE')

    rc = 0

    #
    # Make sure the environment variables are set.
    #
    if not sangerFile:
        print 'Environment variable not set: SANGER_COPY_INPUT_FILE'
        rc = 1

    # Make sure the environment variables are set.
    #
    if not alleleFile:
        print 'Environment variable not set: INPUTFILE'
        rc = 1

    #
    # Initialize file pointers.
    #
    fpSanger = None
    fpAllele = None

    return rc


#
# Purpose: Open files.
# Returns: 1 if file does not exist or is not readable, else 0
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def openFiles():
    global fpSanger, fpAllele

    #
    # Open the Sanger/Biomart file
    #
    try:
        fpSanger = open(sangerFile, 'r')
    except:
        print 'Cannot open file: ' + sangerFile
        return 1

    #
    # Open the Sanger file with genotype sequence #
    #
    try:
        fpAllele = open(alleleFile, 'w')
    except:
        print 'Cannot open file: ' + alleleFile
        return 1

    return 0


#
# Purpose: Close files.
# Returns: 1 if file does not exist or is not readable, else 0
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def closeFiles():

    if fpSanger:
        fpSanger.close()

    if fpAllele:
        fpAllele.close()

    return 0


#
# Purpose: Read the Sanger file and re-format it to create a High-Throughpug MP input file
# Returns: 1 if file does not exist or is not readable, else 0
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def createAlleleFile():

    phenotypingCenter = 'Sangerme'
    annotationCenter = 'Sangerme'
    evidenceCode = 'EXP'

    for line in fpSanger.readlines():

        tokens = line[:-1].split('\t')

#       field 1: Marker Symbol
#       field 2: MGI Marker ID
#	field 6: Mi Attempt Allele Symbol        
#	field 7: Mi Attempt Es Cell Allele Symbol        
#	field 8: Mi Attempt Es Cell MGI Allele Accession 
#	field 9: Mi Attempt Es Cell Name 
#	field 11: Colony Name     
#	field 12: Excision Type   
#	field 13: Tat Cre Phenotype Attempt Deleter Strain        
#	field 16: Phenotype Attempt Production Centre     
#	field 17: MGI Allele Accession    

	sgr_marker_symbol = tokens[0]
        sgr_marker_id = tokens[1]
	sgr_allele_symbol = tokens[5]
	sgr_allele_escell_symbol = tokens[6]
	sgr_allele_id = tokens[7]
	sgr_iscre = tokens[11]
	mgi_allele_id = tokens[16]

	mgi_allele_symbol = 'create this'
	mgi_allele_name = 'create this'
	mgi_allele_status = 'Autoload'

#       field 5:  Allele Generate-Type
#       field 6:  Allele Attribute/SubType
#       field 7:  Germ Line Transmission
#       field 8:  Reference Type/J#
#       field 9:  Strain of Origin
#       field 10  Mutant Cell Line ID
#       field 11: Molecular Notes
#       field 12: Driver Notes
#       field 13: Molecular Mutation
#       field 14: Inheritance
#       field 15: Mixed
#       field 16: Extinct
#       field 17: Creation Date

        fpAllele.write(sgr_marker_id + '\t' + \
                     mgi_allele_symbol + '\t' + \
                     mgi_allele_name + '\t' + \
		     test2 + '\n')

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

