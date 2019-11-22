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


#####################################################################################
#######  1.1 FUNCTION: CHECKING TO SEE IF TAXONID ESEARCH GAVE MULTIPLE HITS  #######
#####################################################################################

def taxidcheck (idlist):


####################################################
#######  2.0 FUNCTION: GETTING PROTEIN NAME  #######
####################################################

def proteindataset (idlist):


###############################################################################################################
#######  3.0 FUNCTION: CREATING SPECIES, TAXID, ACCESSION AND PROTEIN LENGTH LISTS FROM ESEARCH RESULTS #######
###############################################################################################################

def creatinglistfromesearch (home, moment):


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
