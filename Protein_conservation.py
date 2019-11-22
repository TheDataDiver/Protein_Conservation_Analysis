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
	choiceg = ""
	print('\n====================================================================================')
	print('Would you like to:\n\n1. eSearch with:\n\t"TaxonID: %s"\n\t"Protein name: %s"\n\nWARNING: This will perform an extensive search, but conservation plot would not be as biologically significant\n\n\n2. eSearch with:\n\t"TaxonID: %s"\n\t"Protein name: %s"\n\t"Gene name: to be specified"\n\nWARNING: This would produce a biologically significant conservation plot, however less species would be covered due to poor maintenance of NCBI gene name database' % (idlist[1], protchoice, idlist[1], protchoice))
	print('====================================================================================\n')
####### While True loop here as an error trap, in case user does not input one of the available options (1/2)
	while True:
		choice = input('Please type 1 or 2\n').strip()
####### ESEARCH WITH TAXID + PROTEIN + GENE NAME
		if choice == '2':
####### Getting the current hour and minute, then creating a folder using the current hour and minute.
			moment = time.strftime("%H_%M",time.localtime())
			os.mkdir('%s/Assignment2_%s' % (home, moment))
			choiceg = input('\nPlease input the gene name\n').strip()
####### executing the esearch
			subprocess.call("esearch -db protein -query 'txid%s[Organism:exp] AND %s AND %s[Gene Name] NOT PARTIAL' | efetch -format docsum > %s/Assignment2_%s/docsum.txt" % (idlist[1], protchoice, choiceg, home, moment), shell=True)
####### Running function 3.0, which will generate all the lists required to create the panda dataframe from the esearch result
			protspecieslist, protspeciestaxidlist, protspeciesaccessionlist, protlengthlistint = creatinglistfromesearch(home, moment)
			break
####### ESEARCH WITH TAXID + PROT
		elif choice =='1':
####### Getting the current hour and minute, then creating a folder using the current hour and minute.
			moment = time.strftime("%H_%M",time.localtime())
			os.mkdir('%s/Assignment2_%s' % (home, moment))
####### executing the esearch
			subprocess.call("esearch -db protein -query 'txid%s[Organism:exp] AND %s NOT PARTIAL' | efetch -format docsum > %s/Assignment2_%s/docsum.txt" % (idlist[1], protchoice, home, moment), shell=True)
####### Running function 3.0, which will generate all the lists required to create the panda dataframe from the esearch result
			protspecieslist, protspeciestaxidlist, protspeciesaccessionlist, protlengthlistint = creatinglistfromesearch(home, moment)
			break
		else:
			print('Your input was invalid, please input one of the available choices')
####### Creating the panda Dataframe, using the above generated lists
	s1 = pd.Series(protspecieslist)
	s2 = pd.Series(protspeciestaxidlist)
	s3 = pd.Series(protspeciesaccessionlist)
	s4 = pd.Series(protlengthlistint)
####### Dataframe columns are as follows: Species Name, Species TaxID, Protein Accession, Protein length
	df = pd.DataFrame({'Species Name' : s1, 'Species TaxID' : s2, 'Prot Accession' : s3, 'Prot Length' : s4})
####### Getting the total number of sequences that were found, by using the first index of the df.shape
	totalseq = df.shape[0]
####### Getting the total number of unique species
	uniquespecies = len(df.drop_duplicates('Species Name'))
####### IF TOTAL NUMBER OF SEQUENCES FOUND IN ESEARCH IS < 10
	if totalseq < 10:
		print('The esearch found less than 10 total sequences')
		print('\n====================================================================================')
		print('WARNING: Esearch returned with less than 10 unique sequences:\n\nTaxon          :   %s\nProtein Family :   %s\nTotal sequences:   %s\nNo. of Species :   %s' % (idlist[0], protchoice, totalseq, uniquespecies))
		print('\n------------------------------------------------------------------------------------')
####### Creating a list of species names from the panda dataframe
		protspeciesl = df['Species Name'].tolist()
####### Creating a dictionary of Species (KEY) : Sequence counts (VALUE)
		numberoffasta = dict((x, protspeciesl.count(x)) for x in set(protspeciesl))
		print('The number of FASTA sequences per species are as follows:\n')
####### Prints all the unique species and their respective sequence counts
		for key, value in numberoffasta.items():
			print('Species: %-40s' 'Number of FASTA sequences: %s' %(key, numberoffasta[key]))
		print('\n------------------------------------------------------------------------------------\n')
		print('WARNING: The results above contains redundant sequences.\nWARNING: To improve the speed of the script, FASTA sequences will only be downloaded later, and redundant sequences removed.\nWARNING: We shall perform a further check later')
		print('====================================================================================\n')
		return protchoice, protspeciesl, choiceg, df, moment
####### IF TOTAL NUMBER OF UNIQUE SPECIES FOUND IN ESEARCH IS < 5
	if uniquespecies < 5:
		print('\n====================================================================================')
		print('WARNING: Esearch returned less than 5 unique species:\n\nTaxon          :   %s\nProtein Family :   %s\nTotal sequences:   %s\nNo. of Species :   %s' % (idlist[0], protchoice, totalseq, uniquespecies))
		print('\n------------------------------------------------------------------------------------')
####### Creating a list of species names from the panda dataframe
		protspeciesl = df['Species Name'].tolist()
####### Creating a dictionary of Species (KEY) : Sequence counts (VALUE)
		numberoffasta = dict((x, protspeciesl.count(x)) for x in set(protspeciesl))
		print('The number of FASTA sequences per species are as follows:\n')
####### Prints all the unique species and their respective sequence counts
		for key, value in numberoffasta.items():
		 	print('Species: %-40s' 'Number of FASTA sequences: %s' %(key, numberoffasta[key]))
		print('\n------------------------------------------------------------------------------------\n')
		print('WARNING: The results above contains redundant sequences.\nWARNING: To improve the speed of the script, FASTA sequences will only be downloaded later, and redundant sequences removed.\nWARNING: We shall perform a further check later')
		return protchoice, protspecieslist, choiceg, df, moment
####### Prints esearch output summary, showing total number of sequences found and number of unique species, as well as remind the user of his Taxon and Protein search input
	print('\n====================================================================================')
	print('Esearch Results are as follows:\n\nTaxon          :   %s\nProtein Family :   %s\nTotal sequences:   %s\nNo. of Species :   %s' % (idlist[0], protchoice, totalseq, uniquespecies))
	print('\n------------------------------------------------------------------------------------')
####### Creating a list of species names from the panda dataframe
	protspeciesl = df['Species Name'].tolist()
####### Creating a dictionary of Species (KEY) : Sequence counts (VALUE)
	numberoffasta = dict((x, protspecieslist.count(x)) for x in set(protspecieslist))
####### Printing only the top 3 most represented species for the user to see. This is performed through sorting the dictionary according to highest values first. Then appending the top 3 to a list
	print('The top 3 most represented species are as follows:\n')
	highestcountkey = sorted(numberoffasta, key=numberoffasta.get, reverse=True)[:3]
	highestcountvalue = []
	for elem in highestcountkey:
		va = numberoffasta.get(elem)
		highestcountvalue.append(va)
	print('\nSpecies: %-40s Sequences: %s\nSpecies: %-40s Sequences: %s\nSpecies: %-40s Sequences: %s\n' % (highestcountkey[0], highestcountvalue[0], highestcountkey[1], highestcountvalue[1], highestcountkey[2], highestcountvalue[2]))
	print('------------------------------------------------------------------------------------\n')
	print('WARNING: The results above contains redundant sequences.\nWARNING: To improve the speed of the script, FASTA sequences will only be downloaded later, and redundant sequences removed.\nWARNING: We shall perform a further check later')
	print('\n====================================================================================\n')
	viewful = input('\nWould you like to view the full list of species and their respective number of FASTA sequences? y/n\n')
####### If user would like to see the full list of species and their respective sequence counts, this will print the key of the dictionary (species), followed by the count (value)
	if viewful == 'y':
		for key, value in numberoffasta.items():
			print('Species: %-40s' 'Number of FASTA sequences: %s' %(key, numberoffasta[key]))
	return protchoice, protspeciesl, choiceg, df, moment


#############################################################################
####### 3.2 FUNCTION: IS THIS THE DESIRED OUTPUT, IF NO, IF YES THEN  #######
#############################################################################

def changetaxprot (idlist, protchoice, choiceg, home, df, moment):
####### While True loop acts as an error trap, if user does not input one of the available options (y or n)
	while True:
		finalanswer = input('\nWould you like to proceed with the analyis based on the above information? y/n\n').strip().lower()
####### Creating the list of available options for if the user wants to change his esearch input
		changelist = ['1', '2', '3']
		print('\n\n')
####### If finalanswer == 'y', then this ends the function, as this means the user is happy with the esearch output, and his input parameters
		if finalanswer == 'y':
			return protchoice, idlist, choiceg, df, moment
####### If final answer == 'n', then user will be asked which esearch input he would like to change
		elif finalanswer == 'n':
####### If user selected earlier that he wanted to perform an esearch with additional gene name parameter, this will print his current esearch input (Taxon + protein + gene name)
			if choiceg:
				print('\n====================================================================================')
				print('This is your current search parameters:\n\n1: [Taxon] %-25s [TaxonID] %-10s\n2: [Protein] %s\n3: [Gene name] %s' % (idlist[0], idlist[1], protchoice, choiceg))
				print('====================================================================================\n')
				tochange = input('\nWhich of the above would you like to change?\nType 1 for taxon,\nType 2 for Protein+ Gene(optional),\nType 3 for Taxon + Protein + Gene(optional).\n').strip()
####### If user did not select that he wanted to perform an esearch with additional gene name parameter, this will print his current esearch input (Taxon + protein)
			else:
				print('\n====================================================================================')
				print('This is your current search parameters:\n\n1: [Taxon] %-25s [TaxonID] %-10s\n2: [Protein] %s\n3: [Gene name] Unspecified by user' % (idlist[0], idlist[1], protchoice))
				print('====================================================================================\n')
####### While true loop acts as an error trap, in case user does not input one of the available options (1 or 2 or 3)
			while True:
				tochange = input('\nWhich of the above would you like to change?\nEnter the digit:\n\t1 : To change Taxon,\n\t2 : To change Protein + Gene(optional),\n\t3 : To change Taxon + Protein + Gene(optional).\n').strip()
####### Checks user input against the earlier created list ['1', '2', '3']. If input is 1, then run the following functions to change the TaxID
				if tochange == changelist[0]:
					subprocess.call('rm -fr %s/Assignment2_%s' % (home, moment), shell=True)
					idlist, home = taxonid()
					idlist = taxidcheck(idlist)
					protchoice, protspecieslist, choiceg, df, moment = checkforseq(idlist, protchoice, home)
					protchoice, idlist, choiceg, df, moment = changetaxprot(idlist, protchoice, choiceg, home, df, moment)
					return protchoice, idlist, choiceg, df, moment    ###################
####### Checks user input against the earlier created list ['1', '2', '3']. If input is 2, then run the following functions to change the Protein Name and Gene Name (Optional)
				elif tochange == changelist[1]:
					subprocess.call('rm -fr %s/Assignment2_%s' % (home, moment), shell=True)
					protchoice = proteindataset(idlist)
					protchoice, protspecieslist, choiceg, df, moment = checkforseq(idlist, protchoice, home)
					protchoice, idlist, choiceg, df, moment  = changetaxprot(idlist, protchoice, choiceg, home, df, moment)
					return protchoice, idlist, choiceg, df, moment  ################
####### Checks user input against the earlier created list ['1', '2', '3']. If input is 3, then run the following functions to change both the TaxID and Protein Name and Gene Name(optional)
				elif tochange == changelist[2]:
					subprocess.call('rm -fr %s/Assignment2_%s' % (home, moment), shell=True)
					idlist, home = taxonid()
					idlist = taxidcheck(idlist)
					protchoice  = proteindataset(idlist)
					protchoice, protspecieslist, choiceg, df, moment = checkforseq(idlist, protchoice, home)
					protchoice, idlist, choiceg, df, moment = changetaxprot(idlist, protchoice, choiceg, home, df, moment)
					return protchoice, idlist, choiceg, df, moment   ################
				else:
					print('you did not input one of the available option')
		else:
			print('you did not input one of the available options, please try again')
			return protchoice, idlist, choiceg, df, moment

##################################################################################
#######  4.0 FUNCTION: X STANDARD DEVIATIONS ABOVE OR BELOW WARNING PROMPT #######
##################################################################################

def xstandarddeviationwarning (df):
####### Gets the minimum protein length from dataframe and assigns to variable
	min = df['Prot Length'].min()
####### Gets the maximum protein length from dataframe and assigns to variable
	max = df['Prot Length'].max()
####### Calculates the mean protein length from dataframe and assigns to variable
	mean =df['Prot Length'].mean()
####### Calculates the standard deviation of the protein lengths from the dataframe, and assigns to variable
	std = df['Prot Length'].std()
	totala = 0
	totalb = 0
####### Prints the min, max, mean and standard deviation of the protein sequences from the esearch
	print('\n==============================================================================================================================================')
	print('Protein Length Statistics:\n\n\tMinimum Length:                  %s\n\tMaximum Length:                  %s\n\tMean Length:                     %.2f\n\tStandard Deviation:             %.2f' % (min, max, mean, std))
####### IF maximum sequence length is found to be >1 standard deviation above the mean, or >1 standard deviation below the mean, then it will print this warning message
	if max > (mean + std) or min < (mean - std):
		print('\n----------------------------------------------------------------------------------------------------------------------------------------------')
		print('PLEASE READ THIS SECTION:')
		print('\nWARNING: Maximum or minimum length sequence is >1 standarad deviation away from the mean')
		print('\nWARNING: What effect could this have on the output?')
		print('\n\t-If this sequence is significantly longer or shorter than other sequences, it could be an anomoly. Eg.human error while uploading to the NCBIs database')
		print('\n\t-Caution must be exercised by the user when deciding which sequences to remove. Please use the information printed below to make a wise decision')
		print('\n==============================================================================================================================================')
	return min, max, mean, std

##################################################################################################
#######  4.1 FUNCTION: SHOWING HOW MANY SEQUENCES ARE X STANDARD DEVIATIONS ABOVE THE MEAN #######
##################################################################################################

def xstandarddeviationabove (df, min, max, mean, std):
####### While true loop acts as an error trap, in case the user does not input one of the available options for removing sequences >1 standard deviations above/below the mean
	while True:
####### If maximum protein length is greater than 1 standard deviation over the mean, then...
		if max > (mean + std):
			print('WARNING: ||ABOVE THE MEAN|| Maximum length sequence is > 1 standard deviations above the mean\n')
####### Printing the options that the user has
			print('Your options are as follows for sequences |||ABOVE THE MEAN|||:')
			print('\n\t-Choice-\t-Action-')
			print('\t0\t:\t[DO NOT REMOVE]    any sequences')
			stdacheck = mean
			n = 1
####### This while loop searches the dataframe for sequences x standard deviations above the mean, increasing the number of standard deviations each time it loops. Till it reaches 5 standard deviations or finds 0 sequences
			while stdacheck < max and n <= 5:
				stdacheck = mean + std*n
				stda = df[df.apply( lambda x : x['Prot Length'] > (mean + std*n), axis=1 )].shape[0]
				if stda == 0:
					break
				print('\t%s\t:\t[REMOVE SEQUENCES] %s Sequences     that are     %s Standard Deviations above the Mean' % (n, stda, n))
				totala = n
				n += 1
####### Based on the above information, user is asked to input a digit for his choice
			choice = input('\nBased on the above information, please input a digit for your choice\n').strip()
####### Creates a list of options available to the user
			a = range(0, totala + 1)
			stdalist = list(map(str, a))
####### Checks if user input was one of the available options, by referencing the list that was created
			check = any(ele in choice for ele in stdalist)
			if check == True:
####### Choice 0 means that the user does not want to remove any sequences, thus it ends the function here, returning the new dataframe
				if choice == '0':
					print('\n==============================================================================================================================================')
					return df
####### Any other digit signifies x standard deviations above the mean
				else:
####### Turns the user input into int, so that calculations can be made with it
					choice = int(choice)
####### Finds the sequences to be removed and saves it as a new dataframe, and resets the index
					toberemoved = df[df.apply( lambda x : x['Prot Length'] > (mean + std*choice), axis=1 )]
					toberemoved1 = toberemoved.reset_index(drop=True, inplace=True)
####### Removes the selected sequences from the dataframe, saving it to a new dataframe
					ndf = df[~df.apply( lambda x : x['Prot Length'] > (mean + std*choice), axis=1 )]
					print('\n\n')
####### Prints the dataframe of sequences to be removed, and asks the user if he is sure that he wants to remove these sequences
					print(toberemoved)
					choices = input('\nThe above sequences will be removed. Would you like to continue? y/n\n').strip().lower()
####### If answer is y, then ends the function and saves the dataframe with sequences removed
					if choices == 'y':
						print('\n==============================================================================================================================================')
						return ndf
####### If answer is n, then returns back to the start of the function
					elif choices == 'n':
						print('\nLets go back to selecting which sequences are to be removed\n')
					else:
						print('Your input was invalid, please input one of the available options')
			else:
				print('\n----------------------------------------------------------------------------------------------------------------------------------------------')
				print('\nYour input was invalid, we shall try again\n')


##################################################################################################
#######  4.2 FUNCTION: SHOWING HOW MANY SEQUENCES ARE X STANDARD DEVIATIONS BELOW THE MEAN #######
##################################################################################################

def xstandarddeviationbelow (df, min, max, mean, std):
####### While true loop acts as an error trap, in case the user does not input one of the available options for removing sequences >1 standard deviations above/below the mean
	while True:
####### If minimum protein length is greater than 1 standard deviation below the mean, then...
		if min < (mean - std):
####### Printing the options that the user has
			print('WARNING: ||BELOW THE MEAN|| Maximum length sequence is > 1 standard deviations below the mean\n')
			print('Your options are as follows for sequences |||BELOW THE MEAN|||:')
			print('\n\t-Choice-\t-Action-')
			print('\t0\t:\t[DO NOT REMOVE]    any sequences')
			stdbcheck = mean
			n = 1
####### This while loop searches the dataframe for sequences x standard deviations below the mean, increasing the number of standard deviations each time it loops. Till it reaches 5 standard deviations or finds 0 sequences
			while stdbcheck > min and n <= 5:
				stdbcheck = mean - std*n
				stdb = df[df.apply( lambda x : x['Prot Length'] < (mean - std*n), axis=1 )].shape[0]
				if stdb == 0:
					break
				print('\t%s\t:\t[REMOVE SEQUENCES] %s Sequences     that are     %s Standard Deviations below the Mean' % (n, stdb, n))
				totalb = n
				n += 1
####### Based on the above information, user is asked to input a digit for his choice
			choice = input('\nBased on the above information, please input a digit for your choice. Sequences to be removed will then be displayed for confirmation\n').strip()
####### Creates a list of options available to the user
			b = range(0, totalb + 1)
			stdblist = list(map(str, b))
####### Checks if user input was one of the available options, by referencing the list that was created
			check = any(ele in choice for ele in stdblist)
			if check == True:
####### Choice 0 means that the user does not want to remove any sequences, thus it ends the function here, returning the new dataframe
				if choice == '0':
					print('\n==============================================================================================================================================')
					return df
####### Any other digit signifies x standard deviations below the mean
				else:
####### Turns the user input into int, so that calculations can be made with it
					choice = int(choice)
####### Finds the sequences to be removed and saves it as a new dataframe, and resets the index
					toberemoved = df[df.apply( lambda x : x['Prot Length'] < (mean - std*choice), axis=1 )]
					toberemoved1 = toberemoved.reset_index(drop=True, inplace=True)
####### Removes the selected sequences from the dataframe, saving it to a new dataframe
					ndf = df[~df.apply( lambda x : x['Prot Length'] < (mean - std*choice), axis=1 )]           #### new data frame without sequences that were to be removed
					print('\n\n')
####### Prints the dataframe of sequences to be removed, and asks the user if he is sure that he wants to remove these sequences
					print(toberemoved)
					choices = input('\nThe above sequences will be removed. Would you like to continue? y/n\n').strip().lower()
####### If answer is y, then ends the function and saves the dataframe with sequences removed
				if choices == 'y':
					print('\n==============================================================================================================================================')
					return ndf
####### If answer is n, then returns back to the start of the function
				elif choices == 'n':
					print('\n----------------------------------------------------------------------------------------------------------------------------------------------')
					print('\nLets go back to selecting which sequences are to be removed\n')
				else:
					print('Your input was invalid, please input one of the available options')
			else:
				print('\n----------------------------------------------------------------------------------------------------------------------------------------------')
				print('\nYour input was invalid, we shall try again\n')


###########################################################################################
#######  4.3 FUNCTION: DISPLAYING UPDATED DATAFRAME, FOLLOWING REMOVAL OF SEQUNCES  #######
###########################################################################################

def theupdateddataframe (df, idlist, protchoice, home, moment):
####### Gets number of total sequences and unique species from the updated dataframe
	totalseq = df.shape[0]
	uniquespecies = len(df.drop_duplicates('Species Name'))
####### Creates a list of all the species in the updated dataframe
	listofspecies = df['Species Name'].tolist()
####### Prints the summary information of the updated dataframe, being the taxon, protein family, total sequences and number of unique species
	print('\n\n\n====================================================================================')
	print('Filtered esearch results are as follows:\n\nTaxon          :   %s\nProtein Family :   %s\nTotal sequences:   %s\nNo. of Species :   %s' % (idlist[0], protchoice, totalseq, uniquespecies))
####### Creates a dictionary with Species (KEY) : Sequence counts (VALUE)
	numberoffasta = dict((x, listofspecies.count(x)) for x in set(listofspecies))
	print('\n------------------------------------------------------------------------------------')
####### Prints top 3 most represented species
	print('The top 3 most represented species are as follows:\n')
	highestcountkey = sorted(numberoffasta, key=numberoffasta.get, reverse=True)[:3]
	highestcountvalue = []
	for elem in highestcountkey:
		va = numberoffasta.get(elem)
		highestcountvalue.append(va)
	print('\nSpecies: %-40s Sequences: %s\nSpecies: %-40s Sequences: %s\nSpecies: %-40s Sequences: %s\n' % (highestcountkey[0], highestcountvalue[0], highestcountkey[1], highestcountvalue[1], highestcountkey[2], highestcountvalue[1]))
	print('------------------------------------------------------------------------------------\n')
	print('WARNING: The results above contains redundant sequences.\nWARNING: To improve the speed of the script, FASTA sequences will only be downloaded later, and redundant sequences removed.')
	print('\n====================================================================================\n')
####### If user wants to view the entire list of species and their respective sequence counts, then...
	viewful = input('\nWould you like to view the full list of species and their respective number of FASTA sequences? y/n\n')
	if viewful == 'y':
####### Prints the key and the value of the dictionary in a for loop, printing every line
		for key, value in numberoffasta.items():
			print('Species: %-40s' 'Number of FASTA sequences: %s' %(key, numberoffasta[key]))
			print('\n====================================================================================\n')
####### If user decides to go ahead with the analysis, then function returns choice2 and moment(the timestamp)
	choice2 = 'empty'
	while True:
		choice = input('\nWould you like to go ahead with the analysis with this data set? y/n\n').strip().lower()
		if choice == 'y':
			return choice2, moment
####### If user decides not to go ahead with the analysis, then...
		elif choice == 'n':
			choice2 = input('\n\n\nWhich dataset would you like to revert to?\n\n\t1 : Start again, change TaxID and Protein\n\t2 : Revert to dataset before removal of sequences x standarded deviations above/below mean\n\nPlease input one of the digits above\n').strip()
####### Choice 1, deletes the main output folder, as we are requiring a fresh shart, and function returns the choice he has made
			if choice2 == '1':
				subprocess.call("rm -fr %s/Assignment2_%s" % (home, moment), shell=True)
				print('Lets start fresh, here we go again')
				return choice2, moment
####### Choice 2, no need to delete any folder, as user is just reverting to dataframe which was prior to removing sequences above/below standard deviation
			elif choice2 == '2':
				print('Now reverting to dataset prior to removal of sequences x standard deviations above/below the mean')
				return choice2, moment
			else:
				print('Your input was invalid, lets try again')
		else:
			print('You input was invalid, lets try again')


###########################################################
#######  5.0 FUNCTION: DOWNLOADING FASTA SEQUENCES  #######
###########################################################

def downloadingfasta (idlist, protchoice, choiceg, home, moment):
####### Creates a subfolder within main output folder, this is where the main fastafile will be stored
	os.mkdir('%s/Assignment2_%s/fastafiles' % (home, moment))
	print('\n\nFasta sequences are now being downloaded, please be patient')
####### If the user had chose to use the gene name parameter in esearch, then fastafiles will be downloaded with the following parameters
	if choiceg:
		subprocess.call("esearch -db protein -query 'txid%s[Organism:exp] AND %s AND %s[Gene Name] NOT PARTIAL' | efetch -format fasta > %s/Assignment2_%s/fastafiles/unfiltered.fasta" % (idlist[1], protchoice, choiceg, home, moment), shell=True)
####### If not, fastafiles will be downloaded as follows
	else:
		subprocess.call("esearch -db protein -query 'txid%s[Organism:exp] AND %s NOT PARTIAL' | efetch -format fasta > %s/Assignment2_%s/fastafiles/unfiltered.fasta" % (idlist[1], protchoice, home, moment), shell=True)
####### Downloading fasta sequences now helps save time, as user will not have to redownload fasta sequences each time he changes his search parameters


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
