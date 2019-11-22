#!/usr/bin/env python3

import re
import os
import subprocess
import string
from os.path import expanduser
import pandas as pd
import time

#################################################
#######  1.0 FUNCTION: RETRIEVE TAXON ID  #######
#################################################

def taxonid ():
####### Getting the user's home directory path, saving it to the variable 'home'
    home = expanduser("~")
    choice = input('\nPlease specify the taxonomic group\n').strip().lower()
    print('\n')
####### Creation of alphabetical and numerical lists, to check against the user input (to ensure he/she input just taxon name, or taxonID, not both)
    numreflist = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    abcreflist = list(string.ascii_lowercase)
    checkabc = any(ele in choice for ele in abcreflist)
    check123 = any(ele in choice for ele in numreflist)
####### Checks if user input was alphabetical, it will perform esearch using the Taxon name
    if checkabc == True and check123 == False:
        print('\n====================================================================================')
        print('esearch results are as follows:\n')
####### Displays esearch results in genebank format
        subprocess.call("esearch -db taxonomy -query %s | efetch -format Gb" % (choice), shell=True)
        print('====================================================================================\n')
        outp = subprocess.check_output("esearch -db taxonomy -query %s | efetch -format xml | xtract -pattern Taxon -element ScientificName -element TaxId" % (choice), shell=True)
        outputdec = outp.decode("utf-8").strip()
####### Creates a list from the esearch result, containing the taxon and taxonIDs
        taxonidlist = re.split(r"[\t\n]", outputdec)
####### Checks if user input was numerical, it will perform the esearch using the taxon id, with the [UID] parameter
    elif checkabc == False and check123 == True:
        print('\n====================================================================================')
        print('esearch results are as follows:\n')
####### Displays esearch results in genebank format
        subprocess.call("esearch -db taxonomy -query %s[UID] | efetch -format Gb" % (choice), shell=True)
        print('\n====================================================================================')
        outp = subprocess.check_output("esearch -db taxonomy -query %s[UID] | efetch -format xml | xtract -pattern Taxon -element ScientificName -element TaxId" % (choice), shell=True)
        outputdec = outp.decode("utf-8").strip()
####### Creates a list from the esearch result, containing the taxons and taxonIDs
        taxonidlist = re.split(r"[\t\n]", outputdec)
####### Checks if user input was alphanumerical, it won't perform the esearch, returns to initial step, as either the taxonID or Taxon name should be input, not both
    else:
        print('\nInput requires only taxonID or taxonomy, not a mix between both (eg. Either "8782" or "aves")')
        taxonidlist = taxonid()
####### Checking if the esearch had any valid output, if the taxon esearch had 0 results, it will loop back the start of the function
    if len(taxonidlist) == 1:
        print('\nThe esearch for %s taxonomy returned 0 results, the input was probably invalid, we shall return to the initial step' % (choice))
        taxonidlist = taxonid()
    return taxonidlist, home

#####################################################################################
#######  1.1 FUNCTION: CHECKING TO SEE IF TAXONID ESEARCH GAVE MULTIPLE HITS  #######
#####################################################################################

def taxidcheck (idlist):
####### If length of list is greater than 3, then it means that the Taxon esearch generated more than one output
	if len(idlist) >  3:
####### Asks the user to input one of the available options from the esearch output that was printed on the screen.
		choice = input('\nYou have to be more specific, there were more than 1 result! Of the listed results, type the name of the desired output\n').strip().lower().capitalize()
		if choice in idlist:
			for elements in idlist:
####### If his input matched an element from the list, it would delete the other elements from the list, keeping the specified Taxon name, and it's TaxonID
				if elements == choice:
					index = idlist.index(elements)
					del idlist[:index]
					del idlist[index + 2:]
####### If his input was not one of the available options on the list, then loop back to start of function, asking him for his input again
		else:
			print('\nYou did not input one of the available choices')
			idlist = taxidcheck(idlist)
####### On occassion, taxonIDs get updated, this creates a list, with 3 elements, Taxon name, TaxonID new, TaxonID old. This will delete the old taxon ID from the list, keeping the updated one
	elif len(idlist) == 3:
		del idlist[2:]
####### If list length is 2, this means Taxon esearch generated a single taxon, and thus script can continue without further user refinement and input
	elif len(idlist) == 2:
		return idlist
	return idlist


####################################################
#######  2.0 FUNCTION: GETTING PROTEIN NAME  #######
####################################################

def proteindataset (idlist):
####### Asking the user to input his protein of interest, while showing him is selected taxon in the printed statement
	protchoice = input('\nPlease specify the family of Protein which you would like to analyse for the %s taxon\n' % (idlist[0])).strip().lower()
	return protchoice


###############################################################################################################
#######  3.0 FUNCTION: CREATING SPECIES, TAXID, ACCESSION AND PROTEIN LENGTH LISTS FROM ESEARCH RESULTS #######
###############################################################################################################

def creatinglistfromesearch (home, moment):
    def creatinglistfromesearch (home, moment):
    ####### Creating the list of species names from output of Protein + Taxon esearch (Function 3.1)
    	protspecies = subprocess.check_output("cat %s/Assignment2_%s/docsum.txt | xtract -pattern Organism -element Organism"    % (home, moment), shell=True)
    	protspeciesdec = protspecies.decode("utf-8").strip()
    	protspecieslist = re.split(r"[\n]", protspeciesdec)
    ####### Creating the list of taxonIDS from the output of Protein + Taxon esearch (Function 3.1)
    	protspeciestaxid = subprocess.check_output("cat %s/Assignment2_%s/docsum.txt | xtract -pattern DocumentSummary -element TaxId" % (home, moment), shell=True)
    	protspeciestaxiddec = protspeciestaxid.decode("utf-8").strip()
    	protspeciestaxidlist = protspeciestaxiddec.split()
    ####### Creating the list of Protein Accession numbers from the output of Protein + Taxon esearch (Function 3.1)
    	protspeciesaccession = subprocess.check_output("cat %s/Assignment2_%s/docsum.txt | xtract -pattern DocumentSummary -element AccessionVersion" % (home, moment), shell=True)
    	protspeciesaccessiondec = protspeciesaccession.decode("utf-8").strip()
    	protspeciesaccessionlist = protspeciesaccessiondec.split()
    ####### Creating the list of Protein lengths from the output of Protein + Taxon esearch (Function 3.1)
    	protlength = subprocess.check_output("cat %s/Assignment2_%s/docsum.txt | xtract -pattern DocumentSummary -element Slen" % (home, moment), shell=True)
    	protlengthdec = protlength.decode("utf-8").strip()
    	protlengthlist = protlengthdec.split()
    ####### Converting the list of numbers from a string to integers, as stats will be performed on this collumn later, with panda dataframes
    	protlengthlistint = list(map(int, protlengthlist))
    	return protspecieslist, protspeciestaxidlist, protspeciesaccessionlist, protlengthlistint

####################################################################
####### 3.1 FUNCTION: CHECKING SEQUENCES WITH TXID AND PROT  #######
####################################################################

def checkforseq (idlist, protchoice, home):


#############################################################################
####### 3.2 FUNCTION: IS THIS THE DESIRED OUTPUT, IF NO, IF YES THEN  #######
#############################################################################

def changetaxprot (idlist, protchoice, choiceg, home, df, moment):


##################################################################################
#######  4.0 FUNCTION: X STANDARD DEVIATIONS ABOVE OR BELOW WARNING PROMPT #######
##################################################################################

def xstandarddeviationwarning (df):


##################################################################################################
#######  4.1 FUNCTION: SHOWING HOW MANY SEQUENCES ARE X STANDARD DEVIATIONS ABOVE THE MEAN #######
##################################################################################################

def xstandarddeviationabove (df, min, max, mean, std):


##################################################################################################
#######  4.2 FUNCTION: SHOWING HOW MANY SEQUENCES ARE X STANDARD DEVIATIONS BELOW THE MEAN #######
##################################################################################################

def xstandarddeviationbelow (df, min, max, mean, std):


###########################################################################################
#######  4.3 FUNCTION: DISPLAYING UPDATED DATAFRAME, FOLLOWING REMOVAL OF SEQUNCES  #######
###########################################################################################

def theupdateddataframe (df, idlist, protchoice, home, moment):


###########################################################
#######  5.0 FUNCTION: DOWNLOADING FASTA SEQUENCES  #######
###########################################################

def downloadingfasta (idlist, protchoice, choiceg, home, moment):


###################################################################################################################################################################################
#######  6.0 FUNCTION: PULLSEQ TO GENERATE FASTAFILES FOR EACH SPECIES, OUTPUT TO NEW FOLDER, AND EMBOSS SKIPREDUNDANT TO CREATE A NON-REDUNDANT (WITHIN SPECIES) FASTAFILE #######
###################################################################################################################################################################################

def creatingnonredundantfastafiles (df, home, moment):


###################################################################################################################################################################
#######  6.1 FUNCTION: NON-REDUNDANT OUTPUT, SHOWING TOTAL NUMBER OF SEQUENCES, AND TOTAL UNIQUE SPECIES, AS WELL AS A LIST OF EACH SPECIES AND NUMBER OF FASTA #######
###################################################################################################################################################################

def nonredundantcheckforseq (idlist, protchoice, choiceg, home, moment):



###############################################################################################
#######  6.2 FUNCTION: CREATION OF CONSENSUS SEQUENCE FROM THE NON-REDUNDANT FASTA FILE #######
###############################################################################################

def creatingconsensusseq (home, moment):


#######################################################################################################
#######  7.0 FUNCTION: SELECTION OF HOW MANY SEQUENCES THE CONSERVATION PLOT WAS TO BE PERFORMED ON #######
#######################################################################################################

def selectingmaxseqanalysis (totalseq, home, idlist, prot, moment, df):


#######################################################
#######  8.0 FUNCTION: IDENTIFICATION OF MOTIFS #######
#######################################################

def motifanalysis (home, moment, df):


##########################################
#######  RUNNING ALL THE FUNCTIONS #######
##########################################

def runningallfunctions():
