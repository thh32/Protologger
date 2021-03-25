#!/usr/bin/env python2.7

from __future__ import division
from subprocess import call
from datetime import datetime
import HTSeq
import glob
import math
import sys
import subprocess
import os
import operator
import HTSeq
import random
import argparse
import numpy as np
import time
import pickle
import random
import string







complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U':'A'} 

def reverse_complement(seq):    
    bases = list(seq.upper()) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def Similarity(seq1,seq2):
    #GAVLI, FYW, CM, ST, KRH, DENQ, P
    # Same groupings as the Sequence manipulation suite
    total = len(seq1)
    sim = 0
    for q,w in zip(list(seq1), list(seq2)):
        z = q.capitalize()
        x = w.capitalize()
        if z == '-' or x == '-':
            total = total - 1
        if z == 'T' or z == 'U':
            if x == 'U' or x == 'T':
                sim +=1
        if z == x:
            sim +=1
    return sim, total



def randomString(stringLength=8):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))




# Command line options




parser = argparse.ArgumentParser(description="Version 0.98. Protologger is an all-in-one genome description tool, aimed at simplifying the process of gathering the data required for writing protologues. For further information go to the github; https://github.com/thh32/Protologger") #simplifys the wording of using argparse as stated in the python tutorial

# Basic input output files
parser.add_argument("-s", type=str, action='store',  dest='input_16S', help="Input the 16S rRNA gene sequence in FASTA format") # allows input of the forward read
parser.add_argument("-g", type=str, action='store',  dest='input_Genome', help="Input the genome sequence file in FASTA format") # allows input of the forward read
parser.add_argument("-p", type=str, action='store',  dest='project_name', help="Name of your project which will be used for both the input and output folders") # allows input of the forward read


args = parser.parse_args()


# Place each of the input into a simple variable to call
File_16S = str(args.input_16S)
genome_file = str(args.input_Genome)
project_name = str(args.project_name)

if File_16S == 'None' or genome_file == 'None' or project_name == 'None':
    print 'ERROR: One or more of the required inputs are missing, please check you have provided your 16S rRNA gene sequence, genome and project name.'
    if File_16S == 'None':
        print 'ERROR: No 16S rRNA gene sequence file was provided using the -s argument.'
    if genome_file == 'None':
        print 'ERROR: No genome sequence file was provided using the -g argument.'
    if project_name == 'None':
        print 'ERROR: No project name was provided using the -p argument.'
    print 'Protologger will now end.'
    sys.exit()



print 'Input 16S file; ', File_16S
print 'Input genome file; ', genome_file
#project_name = randomString()
print 'Project name is; ', project_name



dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'

print 'Currently running in; ', dir_path



###########
#
#    Create standard input folder named based on the input --- if it already exists it will be deleted
#
######




bashCommand = 'mkdir -p '+ dir_path + 'Input/' + project_name
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


bashCommand = 'cp ' +  File_16S + ' '+ dir_path + 'Input/' + project_name +'/'
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = 'cp ' +  genome_file + ' '+ dir_path + 'Input/' + project_name +'/'
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


bashCommand = 'cp ' + genome_file + ' ' + dir_path + 'Input/' + project_name +'/' + genome_file.split('/')[-1:][0].replace('.dat','.fna').replace('.fa','.fna').replace('.fasta','.fna')
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

File16S =  dir_path + 'Input/' + project_name + '/' + File_16S.split('/')[-1:][0]
genome_file = dir_path + 'Input/' + project_name + '/' + genome_file.split('/')[-1:][0].replace('.dat','.fna').replace('.fa','.fna').replace('.fasta','.fna')



print 'Relocated 16S file; ', File_16S
print 'Relocated genome file; ', genome_file


# Output is named the same as the input folder
bashCommand = 'mkdir -p '+ dir_path + 'Output' 
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = 'mkdir -p '+ dir_path + 'Output/' + project_name
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


bashCommand = 'chmod 777 '+ dir_path + 'Output/' + project_name
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = 'mkdir ' + dir_path + 'Output/'+project_name+'/Genome_analysis' 
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


outputting_overview = open(dir_path + 'Output/'+project_name+'/Overview.txt','w')





















#############################################
# 
#       Annotate 16S
#
#######################






bashCommand = 'blastn -query ' + File_16S + ' -db ' + os.environ['PROTOLOGGER_DATA_DIR'] + '/16S-SILVA/LTP-DB/LTPDB -evalue 0.0000000000000000000000001 -perc_identity 60 -outfmt 6 -out ' + dir_path + 'Output/' + project_name +'/LTP-matches.m8' 
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()












































#############################################
# 
#       Idnetify top hitting 16S
#
#######################


# Read in bitscores for all matches

matches = {}
highest_identity = 0.0
for line in open(dir_path + 'Output/'+project_name+'/LTP-matches.m8'):
    timber = line.replace('\n','').split('\t')
    hit = timber[1].split('-')[0]
    bitscore = float(timber[11])
    identity = float(timber[2])
    if identity > highest_identity:
        highest_identity = identity
    matches[hit] = bitscore

#print '16S matches have been read in'



















# Get top 50 matches based on bitscore

sorted_bits = sorted(matches.items(), key=operator.itemgetter(1), reverse=True)

top_50 = []

num = 0

for i in sorted_bits:
    if num <50:
        top_50.append(i[0])
        num +=1
        
#print 'Identified the 50 top hits'


# Get full names for species

full_names = {}

for line in open(os.environ['PROTOLOGGER_DATA_DIR'] + '/16S-SILVA/LTP-DB/LTPs132_SSU_compressed.fasta','r'):
    if line.startswith('>'):
        timber = line.split('-')
        full_names[timber[0].replace('>','')] = timber[5].split(' subsp')[0] + '--' + timber[6].replace('\n','')


# Check the validity of the species names against the latest DSMZ database (every few months needs updating)
# from https://www.dsmz.de/services/online-tools/prokaryotic-nomenclature-up-to-date/downloads

Valid_names = []

num = 0
for line in open(os.environ['PROTOLOGGER_DATA_DIR'] + '/DSMZ-valid-list/DSMZ-August-2019.tab','r'):
    num +=1
    if num >1:
        timber = line.replace('\n','').split('\t')
        if len(timber) > 1:
            Valid_names.append(timber[0] + '_' + timber[1])
        
#print 'Number of valid species names in DSMZ list; ', len(Valid_names)




# Extract matching sequences

outputting = open(dir_path + 'Output/'+project_name+'/Matching-16S-sequences.fasta','w')

num = 0


for seq in HTSeq.FastaReader(os.environ['PROTOLOGGER_DATA_DIR'] + '/16S-SILVA/LTP-DB/LTPs132_SSU_compressed.fasta'):
    if seq.name.split('-')[0] in top_50:
        num +=1
        #print 'Extracted 16S; ' ,num
        name = full_names[seq.name.split('-')[0]].replace(' ','_')
        if name.split('--')[0] in Valid_names:
            outputting.write('>' + name + '--' + seq.name.split('-')[0] + '--Valid' + '\n')
        else:
            outputting.write('>' + name + '--' + seq.name.split('-')[0] + '--NOT_VALID' + '\n')
        outputting.write(seq.seq + '\n')
        
outputting.close()












































#############################################
# 
#       Calculate completness of 16S
#
#######################

outputting_overview.write('Quality check\n')
outputting_overview.write('-------------\n\n')



# Calculate completeness
top_16S_match = sorted_bits[0][0]


match_16S_length = 0
for read in HTSeq.FastaReader(dir_path + 'Output/'+project_name+'/Matching-16S-sequences.fasta'):
    if read.name.split('--')[2] == top_16S_match:
        match_16S_length = len(read.seq)
        
user_16S_length = 0
for read in HTSeq.FastaReader(File_16S):
    user_16S_length = len(read.seq)
    
completness_16S = (user_16S_length * 100.0)/float(match_16S_length)


if completness_16S > 85.0:
    #print 'The input 16S rRNA gene is acceptable at ' , completness_16S , '% complete compared to its closest relative'
    outputting_overview.write('The input 16S rRNA gene is acceptable at ' + str(completness_16S) + '% complete compared to its closest relative\n')
elif completness_16S < 85.0:
    #print 'WARNING: The input 16S rRNA gene may be incomplete, please check the 16S rRNA gene and be wary of the output results as they may be incorrect'
    #print 'The input 16S rRNA gene was determined to be ' + str(completness_16S)+ '% complete compared to its closest relative'
    outputting_overview.write('WARNING: The input 16S rRNA gene may be incomplete, please check the 16S rRNA gene and be wary of the output results as they may be incorrect\n')
    outputting_overview.write('The input 16S rRNA gene was determined to be ' + str(completness_16S)+ '% complete compared to its closest relative\n')











#############################################
# 
#       Chimera check
#
#######################

#Run UCHIME
bashCommand = 'usearch -uchime ' + File_16S + ' -db ' + os.environ['PROTOLOGGER_DATA_DIR'] + '/16S-SILVA/LTP-DB/LTPs132_SSU_compressed.fasta -uchimeout ' + dir_path + 'Output/' + project_name + '/Chimera_check.txt' 
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

#time.sleep(30)



for line in open(dir_path + 'Output/' +project_name + '/Chimera_check.txt','r'):
    if line.replace('\n','').split('\t')[-1:][0][0] == 'N':
        #print 'The input 16S rRNA gene was not identified as a chimera.\n'
        outputting_overview.write('The input 16S rRNA gene was not identified as a chimera.\n')
    else:
        #print 'WARNING: The input 16S rRNA gene was identified as a chimera, the following results may not be reliable due to this.\n'
        outputting_overview.write('WARNING: The input 16S rRNA gene was identified as a chimera, the following results may not be reliable due to this.\n')








#############################################
# 
#       Check genome quality -- CheckM
#
#######################

#Run GTDB-TK
bashCommand = 'checkm lineage_wf --reduced_tree -x fna ' + dir_path + 'Input/'+project_name+' '+ dir_path  + 'Output/'+project_name+'/CheckM_results' 
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

output, error = process.communicate()

print '1. ', output
print  '2. ',error
print  '3. ',output.split('\n')
print  '4. ',output.split('\n')[-4:][0]
print  '5. ',output.split('\n')[-4:][0].split(' ')

reduced_checkm = []

for i in output.split('\n')[-4:][0].split(' '):
    print i
    if i == '':
        qwep = 0
    else:
        reduced_checkm.append(i)

print  '6. ',reduced_checkm
complete = float(reduced_checkm[12])
contamination = float(reduced_checkm[13])



if complete > 95.0:
    if contamination < 3.0:
        #print 'The input genome is of high quality and suitable for analysis'
        outputting_overview.write('The input genome is of high quality and suitable for analysis'+'\n')
   
    
if contamination > 3.0:
    #print 'WARNING: The input genome is highly contaminated (' + str(contamination) + '%), please check the genome and be wary of the output results as they may be altered due to contamination'
    outputting_overview.write('WARNING: The input genome is highly contaminated (' + str(contamination) + '%), please check the genome and be wary of the output results as they may be altered due to contamination'+'\n')

if complete < 95.0:
    #print 'WARNING: The input genome is incomplete (' + str(complete) + '%), please check the genome and be wary of the output results as they may be altered due to missing genomic regions'
    outputting_overview.write('WARNING: The input genome is incomplete (' + str(complete) + '%), please check the genome and be wary of the output results as they may be altered due to missing genomic regions'+'\n')

bashCommand = 'rm -r ' + dir_path + 'Output/'+project_name+'/CheckM_results' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

outputting_overview.write('Genome completeness is ;' + str(complete) + '\n')
outputting_overview.write('Genome contamination is ;' + str(contamination) + '\n')










outputting_overview.flush()



























#############################################
# 
#       Identify identity between 16S
#
#######################

outputting_overview.write('\n\n16S rRNA gene analysis\n')
outputting_overview.write('-------------------------\n\n')


datas = {}
for prot1 in HTSeq.FastaReader(File_16S):
    #print prot1.name
    protdat = {}
    for prot2 in HTSeq.FastaReader(dir_path + 'Output/'+project_name+'/Matching-16S-sequences.fasta'):
        #print prot2.name
        if prot1.name == prot2.name:
            continue
        else:
            outfile = dir_path + 'Output/' + project_name + '/' + prot1.name + '--' + prot2.name 
            outputting = open(outfile + '.fasta','w')
            outputting.write('>' + prot1.name + '\n')
            outputting.write(prot1.seq + '\n')
            outputting.write('>' + prot2.name + '\n')
            outputting.write(prot2.seq + '\n')
            outputting.close()
            bashCommand = 'clustalw -INFILE=' + outfile + '.fasta -OUTPUT=FASTA -OUTFILE=' + outfile + '.afa'  
            #print bashCommand
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            for read in HTSeq.FastaReader(outfile + '.afa'):
                if read.name == prot1.name:
                    seq1 = read.seq
                elif read.name == prot2.name:
                    seq2 = read.seq
            ##print seq1
            ##print seq2
            sim,total = Similarity(seq1,seq2)
            #print sim
            #print total
            bashCommand = 'rm ' + outfile + '.fasta'  
            #print bashCommand
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            bashCommand = 'rm ' + outfile + '.afa'  
            #print bashCommand
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            bashCommand = 'rm ' + outfile + '.dnd'  
            #print bashCommand
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            percSim = (float(sim)/total)*100
            protdat[prot2.name] = percSim
    datas[prot1.name] = protdat
        


outputting = open(dir_path + 'Output/'+project_name+'/16S-identity-values.tab','w')

outputting.write('#Matching species' + '\t' + 'Similarity (%)' + '\n')

matches = datas.values()[0]

sorted_bits = sorted(matches.items(), key=operator.itemgetter(1), reverse=True)

valls = {}

for i in sorted_bits:
    k = i[0]
    v = i[1]
    valls[k.replace('_',' ').replace('--Valid', ' (Valid)')] = v
    outputting.write(k.replace('_',' ').replace('--Valid', ' (Valid)') + '\t' + str(v) + '\n')
    
outputting.close()






most_similar = max(valls.values())






##############
#
#   Sanity check
#

if most_similar < 70.0:
    output_rev = open(dir_path + 'Output/'+project_name+'/Reverse_complemenet_16S.fasta','w')
    for seq in HTSeq.FastaReader(File_16S):
        output_rev.write('>' + seq.name + '\n')
        output_rev.write(reverse_complement(seq.seq) + '\n')
    output_rev.close()

    rev_16S = dir_path + 'Output/'+project_name+'/Reverse_complemenet_16S.fasta'


    datas = {}
    for prot1 in HTSeq.FastaReader(rev_16S):
        #print prot1.name
        protdat = {}
        for prot2 in HTSeq.FastaReader(dir_path + 'Output/'+project_name+'/Matching-16S-sequences.fasta'):
            #print prot2.name
            if prot1.name == prot2.name:
                continue
            else:
                outfile = dir_path + 'Output/' + project_name + '/' + prot1.name + '--' + prot2.name 
                outputting = open(outfile + '.fasta','w')
                outputting.write('>' + prot1.name + '\n')
                outputting.write(prot1.seq + '\n')
                outputting.write('>' + prot2.name + '\n')
                outputting.write(prot2.seq + '\n')
                outputting.close()
                bashCommand = 'clustalw -INFILE=' + outfile + '.fasta -OUTPUT=FASTA -OUTFILE=' + outfile + '.afa'  
                #print bashCommand
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()
                for read in HTSeq.FastaReader(outfile + '.afa'):
                    if read.name == prot1.name:
                        seq1 = read.seq
                    elif read.name == prot2.name:
                        seq2 = read.seq
                ##print seq1
                ##print seq2
                sim,total = Similarity(seq1,seq2)
                #print sim
                #print total
                bashCommand = 'rm ' + outfile + '.fasta'  
                #print bashCommand
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()
                bashCommand = 'rm ' + outfile + '.afa'  
                #print bashCommand
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()
                bashCommand = 'rm ' + outfile + '.dnd'  
                #print bashCommand
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()
                percSim = (float(sim)/total)*100
                protdat[prot2.name] = percSim
        datas[prot1.name] = protdat
            


    outputting = open(dir_path + 'Output/'+project_name+'/16S-identity-values.tab','w')

    outputting.write('Matching species' + '\t' + 'Similarity (%)' + '\n')

    matches = datas.values()[0]

    sorted_bits = sorted(matches.items(), key=operator.itemgetter(1), reverse=True)

    valls = {}

    for i in sorted_bits:
        k = i[0]
        v = i[1]
        valls[k.replace('_',' ').replace('--Valid', ' (Valid)')] = v
        outputting.write(k.replace('_',' ').replace('--Valid', ' (Valid)') + '\t' + str(v) + '\n')
        
    outputting.close()






    most_similar_rev = max(valls.values())


        

    if most_similar_rev > most_similar:
        #print 'ERROR AVOIDED: Sequence was the reverse complement, corrected.'
        most_similar = most_similar_rev





























novel_level_16S = ''

if most_similar < 86.5:
    novel_level_16S = 'fam. nov.'
elif most_similar < 94.5:
    novel_level_16S = 'gen. nov.'
elif most_similar < 98.7:
    novel_level_16S = 'sp. nov.'
else:
    novel_level_16S = 'not novel'

#print 'The input 16S rRNA determined the isolate to be; ', novel_level_16S
outputting_overview.write('The input 16S rRNA determined the isolate to be; ' + str(novel_level_16S) +'\n')

best_match = ''
for k,v in valls.iteritems():
    if v == most_similar:
        best_match = k
outputting_overview.write('The best match was to; ' + str(best_match)+ ' (' + str(most_similar) +'%)'+'\n')












































#############################################
# 
#       Create 16S rRNA tree
#
#######################



# Create single file

outputting = open(dir_path + 'Output/'+project_name+'/Combined.fasta','w')
for seq in HTSeq.FastaReader(File_16S):
    outputting.write('>' + seq.name + '\n')
    outputting.write(seq.seq + '\n')
for seq in HTSeq.FastaReader(dir_path + 'Output/'+project_name+'/Matching-16S-sequences.fasta'):
    outputting.write('>' + seq.name + '\n')
    outputting.write(seq.seq + '\n')
outputting.close()





# Conduct alignment
bashCommand = 'muscle -in ' + dir_path + 'Output/'+project_name+'/Combined.fasta -out ' + dir_path + 'Output/'+project_name+'/Combined.afa -maxiters 10 -quiet'  
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()
#print '16S tree alignment completed'



outputting = open(dir_path + 'Output/'+project_name+'/16S-tree.nwk','w')

bashCommand = 'fasttree -quiet -gtr -nt ' + dir_path + 'Output/'+project_name+'/Combined.afa'
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

outputting.write(output)
outputting.close()

bashCommand = 'rm ' + dir_path + 'Output/'+project_name+'/Combined.fasta' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()
bashCommand = 'rm ' + dir_path + 'Output/'+project_name+'/Combined.afa' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()







































































outputting_overview.write('\n\nGenome analysis\n')
outputting_overview.write('---------------\n\n')






#############################################
# 
#       Run GTDBTK to identify similar genomes
#
#######################


#Run GTDB-TK

print 'Starting GTDB-TK'


bashCommand = 'gtdbtk classify_wf --genome_dir ' + dir_path + 'Input/'+project_name+' --out_dir ' + dir_path + 'Output/'+project_name+'/GTDB-TK-Results --cpus 20' 
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()
print output 
print error


# Read in meta-data

Representative_genomes = {}

ln = 0
for line in open(os.environ['PROTOLOGGER_DATA_DIR'] + '/GTDB-TK/bac120_metadata_r89.tsv','r'):
    ln +=1
    if ln >1:
        timber = line.replace('\n','').split('\t')
        genome = timber[0].replace('RS_','').replace('GB_','')
        representative = timber[15]
        try:
            taxonomy = timber[78].split(';s__')[1].replace('[','').replace(']','').replace(' ','_')
            if taxonomy == '':
                none = 1
            elif representative == 't':
                Representative_genomes[taxonomy] = genome
        except:
            continue
            
#print 'Representatives of species; ', len(Representative_genomes)




# Read in meta-data

Representative_species = {}

ln = 0
for line in open(os.environ['PROTOLOGGER_DATA_DIR'] + '/GTDB-TK/bac120_metadata_r89.tsv','r'):
    ln +=1
    if ln >1:
        timber = line.replace('\n','').split('\t')
        genome = timber[0].replace('RS_','').replace('GB_','')
        representative = timber[19]
        try:
            taxonomy = timber[78].split(';g__')[1].split(';s_')[0].replace('[','').replace(']','').replace(' ','_')
            if taxonomy == '':
                none = 1
            elif representative == 't':
                Representative_species[taxonomy] = genome
                ##print timber[78]
                ##print representative
        except:
            continue
            
#print 'Representatives of genera; ', len(Representative_species)




# Read in results
ln = 0


Genomes_for_comparison = []

matched_genera = []

for line in open(dir_path + 'Output/'+project_name+'/GTDB-TK-Results/classify/gtdbtk.bac120.summary.tsv','r'):
    ln +=1
    if ln > 1:
        timber = line.replace('\n','').split('\t')
        ##print timber
        #print 'Genome based assignment via GTDBK-TK placed the genome as; ', timber[1]
        outputting_overview.write('Genome based assignment via GTDBK-TK placed the genome as; ' + timber[1] +'\n')
        if timber[2] == 'N/A':
            #print 'GTDB-TK was unable to matched the input genome to that of a previously sequenced genome via FastANI'
            outputting_overview.write('GTDB-TK was unable to matched the input genome to that of a previously sequenced genome via FastANI'+'\n')
        else:
            #print 'GTDB-TK was able to assign the genome to; ' , timber[2]
            outputting_overview.write('GTDB-TK was able to assign the genome to; ' + timber[2] +'\n')
            
        if timber[14] != 'N/A':
            #print 'GTDB-TK close matching genomes identified'
            for i in timber[14].split(';'):
                gID = i.split(',')[0].replace(' ','')
                #print gID
                for k,v in Representative_genomes.iteritems():
                    if v == gID:
                        Genomes_for_comparison.append(gID)
                        matched_genera.append(k.split('_')[0])
        #else:
            #print 'GTDB-TK provided no close matches for analysis'





# Get 16S species matching genomes

for line in open(dir_path + 'Output/'+project_name+'/Matching-16S-sequences.fasta','r'):
    if line.startswith('>'):
        species = line.split('--')[0].replace('>','')
        try:
            Genomes_for_comparison.append(Representative_genomes[species])
            genus = species.split('_')[0]
            if genus in matched_genera:
                l = 1
            else:
                matched_genera.append(genus)
        except:
            continue
            
for i in matched_genera:
    try:
        genome = Representative_species[i]
        if genome in Genomes_for_comparison:
            l = 1
        else:
            Genomes_for_comparison.append(genome)
    except:
        continue
            

#print 'Genomes for comparison; ' + str(len(Genomes_for_comparison))




# Collect genomes

for genome in Genomes_for_comparison:
    taxa = ''
    for k,v in Representative_genomes.iteritems():
        if genome == v:
            taxa = k
            #print genome , taxa
            bashCommand = 'cp ' + os.environ['PROTOLOGGER_DATA_DIR'] + '/Genome-database/' + genome + '_genomic.fna.gz ' + dir_path + 'Output/'+project_name+'/Genome_analysis/' + taxa + '--' + genome + '_genomic.fna.gz' 
            #print bashCommand
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()

    
for i in glob.glob(dir_path + 'Output/'+project_name+'/Genome_analysis/*.gz'):    
    bashCommand = 'gunzip ' + i
    #print bashCommand
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()























outputting_overview.flush()





















#############################################
# 
#       Calculate genome statistics
#
#######################



gc = 0
tot = 0
for read in HTSeq.FastaReader(genome_file):
    for i in read.seq:
        tot +=1
        if i == 'G' or i == 'g' or i == 'C' or i == 'c':
            gc +=1
# Genome size

#print 'The genome size is; ', round(float(tot)/1000000,2) , 'Mbp base pairs'  
outputting_overview.write('The genome size is; ' + str(round(float(tot)/1000000,2)) + 'Mbp base pairs\n')

#print 'G+C percentage is; ' + str(round((float(gc)/tot)*100,2)) +  '%'
outputting_overview.write('G+C percentage is; ' + str(round((float(gc)/tot)*100,2)) + '%\n')


outputting_gc_content = open(dir_path + 'Output/'+project_name+'/GC_of_closest_relatives.tab','w')
outputting_gc_content.write('#Species\tGC%\n')
for genome in glob.glob(dir_path + 'Output/'+project_name+'/Genome_analysis/*.fna'):
    #print genome
    gc = 0
    tot = 0
    for read in HTSeq.FastaReader(genome):
        for i in read.seq:
            tot +=1
            if i == 'G' or i == 'g' or i == 'C' or i == 'c':
                gc +=1

    #print 'G+C percentage is; ' + str(round((float(gc)/tot)*100,2))
    outputting_gc_content.write(genome.split('/')[-1:][0].replace('.fna','').replace('_genomic','') + '\t' + str(round((float(gc)/tot)*100,2)) + '%\n')



outputting_gc_content.close()







































#############################################
# 
#       Calculate ANI statistics
#
#######################


outputting = open(dir_path + 'Output/'+project_name+'/comparison_list.txt','w')

for cfile in glob.glob(dir_path + 'Output/'+project_name+'/Genome_analysis/*.fna'):
    outputting.write(cfile + '\n')
    
outputting.close()


bashCommand = 'fastANI -q '+genome_file+' --rl ' + dir_path + 'Output/'+project_name+'/comparison_list.txt -o ' + dir_path + 'Output/'+project_name+'/ANI_values.tab'
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
time.sleep(120)





### #print highest ANi value

ani_best_match = ''
ani_best_match_perc = 0.0

for line in open(dir_path + 'Output/' + project_name + '/ANI_values.tab','r'):
    timber = line.replace('\n','').split('\t')
    if float(timber[2]) > ani_best_match_perc:
        ani_best_match_perc = float(timber[2])
        ani_best_match = timber[1].split('/')[-1:][0].split('--')[0]

#print 'Best matching species based on ANI was; ' + ani_best_match + '(' + str(ani_best_match_perc) + '%)'

if ani_best_match_perc > 95.0:
    outputting_overview.write('ANI identified the input genome as being a strain of ' + ani_best_match + ' due to an ANI value of ' + str(ani_best_match_perc) + '%\n')
elif ani_best_match_perc > 0.0:
    outputting_overview.write('ANI identified the input genome as novel, with the best match being to ' + ani_best_match + ' with a value of ' + str(ani_best_match_perc) + '%\n')
else:
    outputting_overview.write('Due to the low similarity of the input genome to its closest relatives, no ANI values were able to be calculated.\n')

# Read in meta-data

Representative_genomes = {}

ln = 0
for line in open(os.environ['PROTOLOGGER_DATA_DIR'] + '/GTDB-TK/bac120_metadata_r89.tsv','r'):
    ln +=1
    if ln >1:
        timber = line.replace('\n','').split('\t')
        genome = timber[0].replace('RS_','').replace('GB_','')
        representative = timber[15]
        try:
            taxonomy = timber[78].split(';s__')[1].replace('[','').replace(']','').replace(' ','_')
            if taxonomy == '':
                none = 1
            elif representative == 't':
                Representative_genomes[taxonomy] = genome
        except:
            continue



matches = 0
for line in open(dir_path + 'Output/'+project_name+'/ANI_values.tab','r'):
    timber = line.replace('\n','').split('\t')
    ani = float(timber[2])
    hit = timber[1].split('/')[1].split('_')[0]
    if ani > 95.0:
        matches +=1
        match_species = ''
        for k,v in Representative_genomes.iteritems():
            if v == hit:
                match_species = k
        #print 'ANI match >95% identified to; ' + match_species
        outputting_overview.write('ANI match >95% identified to; ' + match_species+'\n')
        
if matches == 0:
    #print 'No valid species with a sequenced genome within the GTDB-TK database was identified to have a ANI value >95% with the studied genome'
    outputting_overview.write('No valid species with a sequenced genome within the GTDB-TK database was identified to have a ANI value >95% with the studied genome'+'\n')



bashCommand = 'rm ' + dir_path + 'Output/'+project_name+'/comparison_list.txt '
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)





























#############################################
# 
#       Conduct PROKKA annotation and functional analysis
#
#######################

for cfile in glob.glob(dir_path + 'Output/'+project_name+'/Genome_analysis/*.fna'):
    bashCommand = 'mkdir ' + dir_path + 'Output/'+project_name+'/Genome_analysis/'+cfile.split('/')[-1:][0].replace('.fna','')
    #print bashCommand
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    bashCommand = 'prodigal -i ' + cfile + ' -c -m -a ' + dir_path + 'Output/'+project_name+'/Genome_analysis/'+cfile.split('/')[-1:][0].replace('.fna','') + '/' + cfile.split('/')[-1:][0].replace('.fna','') + '.faa'
    #print bashCommand
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()



bashCommand = 'prokka --centre C --locustag L '+genome_file+' -o ' + dir_path + 'Output/'+project_name+'/PROKKA-annotation --prefix Isolate'
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


# Output the number of CRISPR arrays


saved_lines = []
for line in open(dir_path + 'Output/'+project_name+'/PROKKA-annotation/Isolate.log'):
    if 'Found' in line:
        if 'CRISPR' in line:
            number_of_crisprs = line.split(']')[1].split(' ')[2]
            #print 'CRISPR arrays identified; ', number_of_crisprs
            saved_lines.append('CRISPR arrays identified; '+ str(number_of_crisprs)+'\n')

        elif 'CDS' in line:
            number_of_CDS = line.split(']')[1].split(' ')[2]
            #print 'Coding sequences identified; ', number_of_CDS            
            saved_lines.append('Coding sequences identified; '+  str(number_of_CDS) +'\n')














#############################################
# 
#       Create genome based tree
#
#######################
bashCommand = 'mkdir  ' + dir_path + 'Output/'+project_name+'/Genome_Tree' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

for genome in glob.glob(dir_path + 'Output/'+project_name+'/Genome_analysis/*/*.faa'):
    bashCommand = 'cp  '+genome+' ' + dir_path + 'Output/'+project_name+'/Genome_Tree/' 
    #print bashCommand
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

bashCommand = 'cp  ' + dir_path + 'Output/'+project_name+'/PROKKA-annotation/Isolate.faa ' + dir_path + 'Output/'+project_name+'/Genome_Tree/' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = 'phylophlan -i ' + dir_path + 'Output/'+project_name+'/Genome_Tree -d phylophlan --diversity medium  -f ' + os.environ['PROTOLOGGER_DATA_DIR'] + '/Phylophlan/supermatrix_aa.cfg --nproc 20 --output_folder ' + dir_path + 'Output/' + project_name
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = 'cp  ' + dir_path + 'Output/'+project_name+'/Genome_Tree_phylophlan/RAxML_bestTree.Genome_Tree_refined.tre ' + dir_path + 'Output/'+project_name+'/Genome_Tree.nwk' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()



#bashCommand = 'cp -r Output/'+project_name+'/Genome_Tree  bin/phylophlan/input/' + project_name 
##print bashCommand
#process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
#output, error = process.communicate()

#bashCommand = 'python2.7 bin/phylophlan/phylophlan.py -u '+project_name 
##print bashCommand
#process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
#output, error = process.communicate()

#bashCommand = 'cp bin/phylophlan/output/'+project_name+' Output/'+project_name+'/Genome_Tree/'
##print bashCommand
#process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
#utput, error = process.communicate()



#bashCommand = 'rm -r bin/phylophlan/input/'+project_name 
##print bashCommand
#process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
#output, error = process.communicate()
#bashCommand = 'rm -r bin/phylophlan/output/'+project_name 
##print bashCommand
#process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
#output, error = process.communicate()


































#############################################
# 
#       Calculate POCP values between genomes
#
#######################

#outputting_overview.write('\n\nPOCP analysis\n')
#outputting_overview.write('-----------------\n\n')



genomes = []
for genome in glob.glob(dir_path + 'Output/'+project_name+'/Genome_analysis/*/*.faa'):
    genomes.append(genome)

aimed_for = ''
for genome in glob.glob(dir_path + 'Output/'+project_name+'/PROKKA-annotation/*.faa'):
    genomes.append(genome)
    aimed_for = genome
#for i in genomes:
    #print i


    
    
#for i in genomes:
    #print i






    
raw_data = {}
for FILE1 in genomes:
    temp = {}
    for FILE2 in genomes:
        if aimed_for in [FILE1, FILE2]:
            bashCommand = 'blastp -evalue 0.00001 -qcov_hsp_perc 50.0 -subject ' + FILE2 + ' -query ' + FILE1 + ' -outfmt 6 > ' + dir_path + 'Output/' + project_name + '/temp.m8'  
            #print bashCommand
            notneeded = call(bashCommand, shell=True)
            #print notneeded
            total = 0
            for read in HTSeq.FastaReader(FILE1):
                total +=1

            matches = []
            for line in open(dir_path + 'Output/' + project_name + '/temp.m8'):
                timber = line.split('\t')
                if float(line.split('\t')[2]) > 40.0:
                    if line.split('\t')[0] in matches:
                        continue
                    else:
                        matches.append(line.split('\t')[0])
            temp[FILE2] = [total,len(matches)]
    raw_data[FILE1] = temp
        
        

bashCommand = 'rm ' + dir_path + 'Output/'+project_name+'/temp.m8' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()






number_tested = 0
same_genus = 0
different_genus = 0

outputting_pocp = open(dir_path + 'Output/'+project_name+'/POCP_results.tab','w')
outputting_pocp.write('#Matching species\tPOCP (%)\n')
assigned_genus = ''
assigned_pocp = 0.0

for genome1 in glob.glob(dir_path + 'Output/'+project_name+'/PROKKA-annotation/*.faa'):
    for genome2 in genomes:
        if genome1 != genome2:
            print genome1.split('/')[-1:][0].replace('.faa',''), genome2.split('/')[-1:][0].split('_genomic')[0]
            g1d = raw_data[genome1][genome2] #Get the information for genome1 V genome2
            g2d = raw_data[genome2][genome1] #Get the information for genome2 V genome1
            C1 = g1d[1]
            C2 = g2d[1]
            T1 = g1d[0]
            T2 = g2d[0]
            #print C1, C2
            #print T1,T2
            PCOP = ((C1+float(C2))/(T1+float(T2)))*100
            #print PCOP
            number_tested +=1
            if PCOP > 50.0:
                same_genus +=1
                print 'POCP suggests the input genome and ', genome2.split('/')[-1:][0].split('_genomic')[0] , ' belong to the same genus with a POCP of ',  str(PCOP) 
                outputting_pocp.write('POCP suggests the input genome and '+ genome2.split('/')[-1:][0].split('_genomic')[0] + ' belong to the same genus with a POCP was;\t' + str(PCOP) + '%\n')
                assigned_genus = genome2.split('/')[-1:][0].split('_genomic')[0].split('_')[0]
                assigned_pocp = PCOP
            else:
                print 'POCP suggests seperate genus between the input genome and ', genome2.split('/')[-1:][0].split('_genomic')[0]
                outputting_pocp.write('POCP suggests seperate genus between the input genome and '+ genome2.split('/')[-1:][0].split('_genomic')[0] + ' as POCP was;\t' + str(PCOP) + '%\n')

outputting_pocp.close()

if number_tested == same_genus:
    print 'The input genome belongs to the same genus as all database species, as determined by POCP values >50%'
    outputting_overview.write('The input genome belongs to the same genus as all database species, as determined by POCP values >50%'+'\n')
elif same_genus == 0:
    print 'No POCP values >50% were identified, indicating the input genome represents a novel genus'
    outputting_overview.write('No POCP values >50% were identified, indicating the input genome represents a novel genus'+'\n')
else:
    outputting_overview.write('The input genome was assigned to ' + assigned_genus +' with a POCP value of ' + str(assigned_pocp) +'%\n')














































































#############################################
# 
#       Identify if PHAGE are present
#
#######################
##print 'Phage search; started'
#bashCommand = 'wget --post-file="'+genome_file+'" "http://phaster.ca/phaster_api?contigs=1" -O Output/' + project_name + '/Phaster_status.txt'
##print bashCommand
#process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)


#completed = False



#while completed == False:
#    with open("Output_" + genome_file.split('/')[2].replace('.fna','')) as f:
#        for line in f:
#            timber = line.split(',')
#            #print timber[1]
#            if timber[1] == '"status":"Complete"':
#                completed = True
#            else:
#                jobid = timber[0].split('"')[-1:][0]
#                bashCommand = 'bash check_status-PHASTER.sh ' + jobid + ' ' + genome_file
#                #print bashCommand
#                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
#    time.sleep(600) # Check every 10 minutes
        
#phages = 0
#for i in line.split('\\n'):
#    if i.startswith('Totally'):
#        phages = int(i.split(' ')[1])
##print 'Number of intact phage detected; ', phages
#outputting_overview.write('Number of intact phage detected; ' + str(phages)+'\n')

##print 'Phage search; complete'



#bashCommand = 'rm Output_' + genome_file.split('/')[2].replace('.fna','')
##print bashCommand
#process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

outputting_overview.flush()




































for i in saved_lines:
    outputting_overview.write(i)





# Run PRokka2KEGG

#Incase an update is needed BUT DONT BECAUSE THE MAPPINGS MAY CHANGE!!!!
#  wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz

outputting_overview.write('\n\nFunctional analysis\n')
outputting_overview.write('-------------------\n\n')




bashCommand = 'python3.6 ' + os.environ['PROTOLOGGER_DATA_DIR'] + '/KEGG-analysis/prokka2kegg.py -i ' + dir_path + 'Output/'+project_name+'/PROKKA-annotation/Isolate.gbk -d ' + os.environ['PROTOLOGGER_DATA_DIR'] + '/KEGG-analysis/idmapping_KO.tab.gz -o ' + dir_path + 'Output/'+project_name+'/KEGG_profile.KOs' 
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()





KOs = []
for line in open(dir_path + 'Output/'+project_name+'/KEGG_profile.KOs','r'):
    timber = line.replace('\n','').split('\t')
    if len(timber) > 1:
        ko = timber[1]
        KOs.append(ko)


# Read in KEGG groupings

brite_groups = {}

for line in open(os.environ['PROTOLOGGER_DATA_DIR'] + '/KEGG-analysis/brite.map','r'):
    timber = line.replace('\n','').split('\t')
    brite = timber[0].replace('br:','')
    ko = timber[1].replace('ko:','')
    if brite in brite_groups:
        prev = brite_groups[brite]
        prev.append(ko)
        brite_groups[brite] = prev
    else:
        brite_groups[brite] = [ko]


# Transporters
transporters = 0
for ko in brite_groups['ko02000']:
    if ko in KOs:
        ##print ko
        transporters +=1
    
#print 'Number of transporters; ' + str(transporters)
outputting_overview.write('Number of transporters; ' + str(transporters)+'\n')





secrete = 0
for ko in brite_groups['ko02044']:
    if ko in KOs:
        ##print ko
        secrete +=1
    
#print 'Number of secretion genes; ' + str(secrete)
outputting_overview.write('Number of secretion genes; ' + str(secrete)+'\n')




enzymes = 0
for ko in brite_groups['ko01000']:
    if ko in KOs:
        ##print ko
        enzymes +=1
    
#print 'Number of unique enzymes; ' + str(enzymes)
outputting_overview.write('Number of unique enzymes; ' + str(enzymes)+'\n')

# Anoxygenic photosystem I

if 'K08940'in KOs:
	if 'K08941' in KOs:
		if 'K08942' in KOs:
			if 'K08943' in KOs:
				outputting_overview.write('Anoxygenic photosystem I was predicted from the genome (pscA, pscB, pscC, pscD)\n')

# Anoxygenic photosystem II

if 'K08928'in KOs:
	if 'K08929' in KOs:
		outputting_overview.write('Anoxygenic photosystem II was predicted from the genome (pufL, pufM)\n')

		
		
# Methanogenesis, CO2 -> methane
methanogensis = []
if 'K00200' in KOs:
	if 'K00201' in KOs:
		if 'K00202' in KOs:
			if 'K00203' in KOs:
				if 'K00205' in KOs or 'K11260' in KOs or 'K00204' in KOs:
					methanogensis.append('1.2.99.5')

if 'K00672' in KOs:
	methanogensis.append('2.3.1.101')

if 'K01499' in KOs:
	methanogensis.append('3.5.4.27')

if 'K00319' in KOs:
	methanogensis.append('1.5.99.9')
if 'K13942' in KOs:
	methanogensis.append('1.12.98.2')


if 'K00320' in KOs:
	methanogensis.append('1.5.99.11')
	
if 'K00577' in KOs:
	if 'K00578' in KOs:
		if 'K00579' in KOs:
			if 'K00580' in KOs:
				if 'K00581' in KOs:
					methanogensis.append('2.1.1.86')
					
if 'K00399' in KOs:
	if 'K00401' in KOs:
		if 'K00402' in KOs:
			methanogensis.append('2.8.4.1')

if len(methanogensis) == 7:
	outputting_overview.write('Methane production from CO2 was predicted from the genome (EC:1.2.99.5, 2.3.1.101, 3.5.4.27, 1.5.99.9, 1.12.98.2, 1.5.99.11, 2.1.1.86, 2.8.4.1)\n')


# Reductive pentose phosphate cycle

if 'K00855' in KOs: 
	if 'K01601' in KOs: 
		if 'K00927' in KOs: 
			if 'K05298' in KOs and 'K00150' in KOs and 'K00134' in KOs: 
				if 'K01623' in KOs and 'K01624' in KOs:
					if 'K03841' in KOs and 'K02446' in KOs and 'K11532' in KOs and 'K01086' in KOs:
						if 'K00615' in KOs: 
							if 'K01623' in KOs and 'K01624' in KOs:
								if 'K01100' in KOs and 'K11532' in KOs and 'K01086' in KOs:
									if 'K00615' in KOs: 
										if 'K01807' in KOs and 'K01808' in KOs:
											outputt.write('The reductive pentose phosphate cycle (Calvin cycle) was predicted in the genome (EC:2.7.1.19, 4.1.1.39, 2.7.2.3, 1.2.1.13, 1.2.1.59, 1.2.1.12, 4.1.2.13, 3.1.3.11, 2.2.1.1, 3.1.3.37, 5.3.1.6).')




# Siroheme biosynthesis

if 'K01885' in KOs or 'K14163' in KOs:
    if 'K02492' in KOs:
        if 'K01845' in KOs:
            if 'K01698' in KOs:
                if 'K01749' in KOs:
                    if 'K01719' in KOs or 'K13542' in KOs or 'K13543' in KOs:
                        if 'K02302' in KOs:
                            #print 'Siroheme biosynthesis from glutamate was identified from the genome (EC:6.1.1.17, 1.2.1.70, 5.4.3.8, 4.2.1.24, 2.5.1.61, 4.2.1.75, 2.1.1.107/1.3.1.76/4.99.1.4)'
                            outputting_overview.write('Siroheme biosynthesis from glutamate was predicted from the genome (EC:6.1.1.17, 1.2.1.70, 5.4.3.8, 4.2.1.24, 2.5.1.61, 4.2.1.75, 2.1.1.107/1.3.1.76/4.99.1.4)\n')

                        elif 'K02304' in KOs:
                            if 'K03794' in KOs:
                                if 'K00589' in KOs or 'K02303' in KOs or 'K02496' in KOs or 'K13542' in KOs or 'K13543' in KOs:
                                    #print 'Siroheme biosynthesis from glutamate was identified from the genome (EC:6.1.1.17, 1.2.1.70, 5.4.3.8, 4.2.1.24, 2.5.1.61, 4.2.1.75, 2.1.1.107, 1.3.1.76, 4.99.1.4)'
                                    outputting_overview.write('Siroheme biosynthesis from glutamate was predicted from the genome (EC:6.1.1.17, 1.2.1.70, 5.4.3.8, 4.2.1.24, 2.5.1.61, 4.2.1.75, 2.1.1.107, 1.3.1.76, 4.99.1.4)\n')




# Carbon source utilisation
carbon_sources = []

if 'K02777' in KOs or 'K20116' in KOs:
    if 'K02790' in KOs or 'K02778' in KOs or 'K20117' in KOs:
        if 'K02791' in KOs or 'K02779' in KOs or 'K20118' in KOs:
            carbon_sources.append('glucose')
            

if 'K02777' in KOs or 'K02752' in KOs or 'K02753' in KOs:
    if 'K01222' in KOs or 'K01223' in KOs:
        carbon_sources.append('arbutin')
        carbon_sources.append('salicin')
        
if 'K02759' and 'K02760' and 'K02761' in KOs:
    if 'K01222' and 'K01223' in KOs:
        carbon_sources.append('cellobiose')

if 'K02808' and 'K02809' and 'K02810' in KOs:
    carbon_sources.append('sucrose')

if 'K02777' in KOs or 'K02817' in KOs:
    if 'K02818' in KOs:
        if 'K02819' in KOs:
            carbon_sources.append('trehalose')
            
if 'K02777' in KOs:
    if 'K02790' in KOs or 'K20107' in KOs or 'K02749' in KOs:
        if 'K02791' in KOs or 'K20108' in KOs or 'K02750' in KOs:
            carbon_sources.append('maltose')
            
            
if 'K00688' in KOs or 'K16153' in KOs or 'K01196' in KOs or 'K00705' in KOs or 'K22451' in KOs or 'K02438' in KOs or 'K01200' in KOs:
    carbon_sources.append('starch')

if 'K05988' in KOs:
    carbon_sources.append('dextran')
    
    
if 'K01225' in KOs or 'K19668' in KOs:
    if 'K01188' in KOs or 'K05349' in KOs or 'K05350' in KOs or 'K00702':
        carbon_sources.append('cellulose')
elif 'K01179' in KOs or 'K19357' in KOs:
    if 'K01188' in KOs or 'K05349' in KOs or 'K05350':
        carbon_sources.append('cellulose')
    elif 'K00702' in KOs:
        carbon_sources.append('cellulose')
        
    
#print 'The following carbon sources were identified to be utilised; ', ', '.join(carbon_sources)
outputting_overview.write('The following carbon sources were predicted to be utilised; ' + ', '.join(carbon_sources) +'\n')



# VFA production

if 'K00625' in KOs or 'K13788' in KOs or 'K15024' in KOs:
    if 'K00925' in KOs:
        #print 'Acetate production identified from acetyl-CoA (EC:2.3.1.8, 2.7.2.1)'
        outputting_overview.write('Acetate production predicted from acetyl-CoA (EC:2.3.1.8, 2.7.2.1)\n')
      
        
if 'K01034' in KOs or 'K01035' in KOs or 'K19709' in KOs:
    #print 'Butanoate production from butanoyl-CoA (EC:2.8.3.8)'
    outputting_overview.write('Butyrate production predicted from butanoyl-CoA (EC:2.8.3.8)\n')
    
if 'K00634' in KOs:
    if 'K00929' in KOs:
        #print 'Butanoate production from butanoyl-CoA (EC:2.3.1.19, 2.7.2.7)'
        outputting_overview.write('Butyrate production from butanoyl-CoA (EC:2.3.1.19, 2.7.2.7)\n')

        
if 'K13788' in KOs or 'K00625' in KOs or 'K15024' in KOs:
    if 'K00925' in KOs:
        #print 'Propanoate production from propanoyl-CoA (EC:2.3.1.8, 2.7.2.1)'
        outputting_overview.write('Propionate production predicted from propanoyl-CoA (EC:2.3.1.8, 2.7.2.1)\n')
    elif 'K00932' in KOs or 'K19697' in KOs:
        #print 'Propanoate production from propanoyl-CoA (EC:2.3.1.8, 2.7.2.15)'
        outputting_overview.write('Propionate production predicted from propanoyl-CoA (EC:2.3.1.8, 2.7.2.15)\n')
elif 'K13923' in KOs:
    if 'K00925' in KOs:
        #print 'Propanoate production from propanoyl-CoA (EC:2.3.1.222, 2.7.2.1)'
        outputting_overview.write('Propionate production predicted from propanoyl-CoA (EC:2.3.1.222, 2.7.2.1)\n')
    elif 'K00932' in KOs or 'K19697' in KOs:
        #print 'Propanoate production from propanoyl-CoA (EC:2.3.1.222, 2.7.2.15)' 
        outputting_overview.write('Propionate production predicted from propanoyl-CoA (EC:2.3.1.222, 2.7.2.15)\n')
elif 'K01026' in KOs:
    #print 'Propanoate production from propanoyl-CoA (EC:2.8.3.1)' 
    outputting_overview.write('Propionate production predicted from propanoyl-CoA (EC:2.8.3.1)\n')
elif 'K01905' in KOs:
    #print 'Propanoate production from propanoyl-CoA (EC:6.2.1.13)' 
    outputting_overview.write('Propionate production predicted from propanoyl-CoA (EC:6.2.1.13)\n')




# Flagellar proteins
flagella_proteins = []

if 'K02400' in KOs:
    flagella_proteins.append('FlhA')
    
if 'K02401' in KOs or 'K13820' in KOs:
    flagella_proteins.append('FlhB')
    
    
    
    
if 'K02387' in KOs:
    flagella_proteins.append('FlgB')
    
if 'K02388' in KOs:
    flagella_proteins.append('FlgC')
    
if 'K02389' in KOs:
    flagella_proteins.append('FlgD')
    
if 'K02390' in KOs:
    flagella_proteins.append('FlgE')
    
if 'K02391' in KOs:
    flagella_proteins.append('FlgF')
    
if 'K02392' in KOs:
    flagella_proteins.append('FlgG')
    
if 'K02393' in KOs:
    flagella_proteins.append('FlgH')
    
if 'K02394' in KOs:
    flagella_proteins.append('FlgI')
    
if 'K02395' in KOs:
    flagella_proteins.append('FlgJ')
    
if 'K02396' in KOs:
    flagella_proteins.append('FlgK')
    
if 'K02397' in KOs:
    flagella_proteins.append('FlgL')
    
    
if 'K02406' in KOs:
    flagella_proteins.append('FliC')    
    
if 'K02407' in KOs:
    flagella_proteins.append('FliD')
    
if 'K02408' in KOs:
    flagella_proteins.append('FliE')
    
if 'K02409' in KOs:
    flagella_proteins.append('FliF')
    
if 'K02410' in KOs:
    flagella_proteins.append('FliG')
    
if 'K02414' in KOs:
    flagella_proteins.append('FliK')

if 'K02416' in KOs:
    flagella_proteins.append('FliM')
    
if 'K02417' in KOs:
    flagella_proteins.append('FliN')
    
    

if 'K02556' in KOs:
    flagella_proteins.append('MotA')
    
if 'K02557' in KOs:
    flagella_proteins.append('MotB')
    
if 'K21217' in KOs:
    flagella_proteins.append('MotX')
    
if 'K21218' in KOs:
    flagella_proteins.append('MotY')
    

if len(flagella_proteins) > 15:
    #print 'The following flagellar proteins were identified within the genome; ', ','.join(flagella_proteins)
    outputting_overview.write('The following flagellar proteins were identified within the genome; '+ ','.join(flagella_proteins)+'\n')


# Urease added
if 'K01428' in KOs and 'K01429' in KOs and 'K01430' in KOs:
    #print 'The urease cluster (alpha, beta and gamma subunits) were identified within the genome (EC:3.5.1.5)'
    outputting_overview.write('The urease cluster (alpha, beta and gamma subunits) were identified within the genome (EC:3.5.1.5)\n')



# Sulfide conversion to L-cysteine
if 'K00640' in KOs or 'K22304' in KOs:
    if 'K01738' in KOs or 'K13034' in KOs or 'K10150' in KOs or 'K17069' in KOs:
        #print 'Sulfide and L-serine utilised to produce L-cysteine and acetate (EC:2.3.1.30, 2.5.1.47)'
        outputting_overview.write('Sulfide and L-serine are predicted to be utilised to produce L-cysteine and acetate (EC:2.3.1.30, 2.5.1.47)\n')
    elif 'K10150' in KOs:
        #print 'Sulfide and L-serine utilised to produce L-cysteine and acetate (EC:2.3.1.30, 2.5.1.65)'
        outputting_overview.write('Sulfide and L-serine are predicted to be utilised to produce L-cysteine and acetate (EC:2.3.1.30, 2.5.1.65)\n')




# Assimilatory sulfate reduction

if 'K00955' in KOs or 'K00956' in KOs or 'K00957' in KOs or 'K00958' in KOs or 'K13811' in KOs:
    if 'K13811' in KOs or 'K00955' in KOs or 'K00860' in KOs:
        if 'K00390' in KOs:
            if 'K00380' in KOs or 'K00381' in KOs:
                #print 'Sulfate assimilatory reduction to sulfide identified via the genome (EC:2.7.7.4, 2.7.1.25, 1.8.4.8, 1.8.1.2)'
                outputting_overview.write('Sulfate assimilatory reduction to sulfide predicted via the genome (EC:2.7.7.4, 2.7.1.25, 1.8.4.8, 1.8.1.2)\n')
            elif 'K00392' in KOs:
                #print 'Sulfate assimilatory reduction to sulfide identified via the genome (EC:2.7.7.4, 2.7.1.25, 1.8.4.8, 1.8.7.1)'
                outputting_overview.write('Sulfate assimilatory reduction to sulfide predicted via the genome (EC:2.7.7.4, 2.7.1.25, 1.8.4.8, 1.8.7.1)\n')



# Ammonia utilisation

if 'K01915' in KOs:
    if 'K00265' in KOs or 'K00266' in KOs or 'K00264' in  KOs or 'K00284' in KOs:
        #print 'L-glutamate production from ammonia was identified via L-glutamine (EC:6.3.1.2, 1.4.1.-)'
        outputting_overview.write('L-glutamate production from ammonia was predicted via L-glutamine (EC:6.3.1.2, 1.4.1.-)\n')



# Nitrogen fixation

if 'K02586' in KOs or 'K02591' in KOs or 'K02588' in KOs:
    if 'K00531' in KOs:
        if 'K22896' in KOs or 'K22897' in KOs or 'K22899' in KOs:
            #print 'Nitrogen fixation via the production of ammonia identiied (EC:1.18.6.1)'
            outputting_overview.write('Nitrogen fixation via the production of ammonia was predicted (EC:1.18.6.1)\n')




# Nitrification

if 'K10944' in KOs and 'K10945' in KOs and 'K10946' in KOs:
    if 'K10535' in KOs:
        if 'K00370' in KOs and 'K00371' in KOs:
            #print 'Nitrification identified via conversion of ammonia into nitrate (EC:1.14.99.39, 1.7.2.6, 1.7.5.1)'
            outputting_overview.write('Nitrification predicted via conversion of ammonia into nitrate (EC:1.14.99.39, 1.7.2.6, 1.7.5.1)\n')




############################
#
# Vitamin biosynthesis
#
#####################


# Biotin biosynthesis

if 'K00652' in KOs:
    if 'K19562' in KOs:
        if 'K01012' in KOs:
            #print 'Biotin biosynthesis from pimeloyl-ACP/CoA identified (EC:2.3.1.47, 2.6.1.62, 6.3.3.3, 2.8.1.6)'
            outputting_overview.write('Biotin biosynthesis predicted from pimeloyl-ACP/CoA (EC:2.3.1.47, 2.6.1.62, 6.3.3.3, 2.8.1.6)\n')
    elif 'K00833' in KOs or 'K19563' in KOs:
        if 'K01935' in KOs:
            if 'K01012' in KOs:
                #print 'Biotin (vitamin B7) biosynthesis from pimeloyl-ACP/CoA identified (EC:2.3.1.47, 2.6.1.62, 6.3.3.3, 2.8.1.6)'
                outputting_overview.write('Biotin (vitamin B7) biosynthesis predicted from pimeloyl-ACP/CoA (EC:2.3.1.47, 2.6.1.62, 6.3.3.3, 2.8.1.6)\n')
               
# Cobalamin biosynthesis
if 'K00798' in KOs or 'K19221' in KOs:
    if 'K02232' in KOs:
        if 'K02225' in KOs or 'K02227' in KOs:
            if 'K02231' in KOs:
                #print 'Cobalamin (vitamin B12) biosynthesis from cobinamide identified (EC:2.5.1.17, 6.3.5.10, 6.2.1.10, 2.7.1.156)'
                outputting_overview.write('Cobalamin (vitamin B12) biosynthesis predicted from cobinamide (EC:2.5.1.17, 6.3.5.10, 6.2.1.10, 2.7.1.156)\n')
    
# Folate biosynthesis   
if 'K00287' in KOs or 'K18589' in KOs or 'K19643' in KOs or 'K18590' in KOs or 'K19644' in KOs or 'K18591' in KOs or 'K13998' in KOs:
    #print 'Folate (vitamin B9) biosynthesis from 7,8-dihydrofolate identified (EC:1.5.1.3)'
    outputting_overview.write('Folate (vitamin B9) biosynthesis predicted from 7,8-dihydrofolate (EC:1.5.1.3)\n')
    
# Riboflavin biosynthesis
if  'K01497' in KOs or 'K14652' in KOs:
    if 'K11752' in KOs:
        if 'K22912' in KOs or 'K20860' in KOs or 'K20861' in KOs or 'K20862' in KOs or 'K21063' in KOs or 'K21064' in KOs:
            if 'K14652' in KOs or 'K02858' in KOs:
                if 'K00794' in KOs:
                    if 'K00793' in KOs:
                        #print 'Riboflavin (vitamina B2) biosynthesis from GTP identified (EC:3.5.4.25, 3.5.4.26, 1.1.1.193, 3.1.3.104, 4.1.99.12, 2.5.1.78, 2.5.1.9, 2.7.1.26, 2.7.7.2)'
                        outputting_overview.write('Riboflavin (vitamin B2) biosynthesis predicted from GTP (EC:3.5.4.25, 3.5.4.26, 1.1.1.193, 3.1.3.104, 4.1.99.12, 2.5.1.78, 2.5.1.9, 2.7.1.26, 2.7.7.2)\n')
    elif 'K01498' in KOs:
        if 'K00082' in KOs:
            if 'K22912' in KOs or 'K20860' in KOs or 'K20861' in KOs or 'K20862' in KOs or 'K21063' in KOs or 'K21064' in KOs:
                if 'K14652' in KOs or 'K02858' in KOs:
                    if 'K00794' in KOs:
                        if 'K00793' in KOs:
                            #print 'Riboflavin (vitamina B2) biosynthesis from GTP identified (EC:3.5.4.25, 3.5.4.26, 1.1.1.193, 3.1.3.104, 4.1.99.12, 2.5.1.78, 2.5.1.9, 2.7.1.26, 2.7.7.2)'
                            outputting_overview.write('Riboflavin (vitamin B2) biosynthesis predicted from GTP (EC:3.5.4.25, 3.5.4.26, 1.1.1.193, 3.1.3.104, 4.1.99.12, 2.5.1.78, 2.5.1.9, 2.7.1.26, 2.7.7.2)\n')
elif 'K14652' in KOs or 'K02858' in KOs:
    if 'K00794' in KOs:
        if 'K00793' in KOs:
            #print 'Riboflavin (vitamina B2) biosynthesis from D-ribulose-5-phosphate identified (EC:4.1.99.12, 2.5.1.78, 2.5.1.9, 2.7.1.26, 2.7.7.2)'
            outputting_overview.write('Riboflavin (vitamin B2) biosynthesis predicted from D-ribulose-5-phosphate (EC:4.1.99.12, 2.5.1.78, 2.5.1.9, 2.7.1.26, 2.7.7.2)\n')



# Vitamin K2 (Menaquinone)
if 'K02548' in KOs:
    if 'K03183' in KOs:
        if 'K10106' in KOs:
            if 'K05357' in KOs:
                outputting_overview.write('Menaquinone (vitamin K2) biosynthesis predicted from polyprenyl-diphosphate (EC:2.5.1.74, 2.1.1.163, 4.1.1.90, 1.17.4.4)\n')
                
# Vitamin K1 (Phylloquinone)

if 'K23094' in KOs:
    if 'K17872' in KOs:
        if 'K23095' in KOs:
            if 'K10106' in KOs:
                if 'K05357' in KOs:
                    outputting_overview.write('Phylloquinone (vitamin K1) biosynthesis predicted from phytyl-diphosphate (EC:2.5.1.130, 1.6.5.12, 2.1.1.329, 4.1.1.90, 1.17.4.4)\n')
   


# Detection of EPS biosynthesis

if 'K19667' in KOs:
    if 'K01991' in KOs:
        if 'K16692' in KOs:
            if 'K01104' in KOs:
                if 'K01791' in KOs:
                    if 'K02472' in KOs:
                        if 'K16712' in KOs:
                            if 'K16713' in KOs:
                                    outputting_overview.write('The EPS gene cluster (EpsA, EpsB, EpsC, EpsD, EpsE, EpsF, EpsP) was predicted to be present along with the transcriptional activator of the cluster (XpsR)\n')
   
# Dection of cbb3 type cytochrome c oxidase
if 'K00404' in KOs:
    if 'K00405' in KOs:
        if 'K00407' in KOs:
            if 'K00406' in KOs:
                            outputting_overview.write('Cbb3-type cytochrome C oxidase was predicted based on the presence of sub-units I, II, III and IV\n')
elif 'K15862' in KOs:
    if 'K00407' in KOs:
        if 'K00406' in KOs:
                        outputting_overview.write('Cbb3-type cytochrome C oxidase was predicted based on the presence of sub-units I/II (ccoNO), III and IV\n')



# Detect sporulation genes

sporulation = []

spore_db = {'spo0B'   :   'K06375',
'spo0F'    :   'K02490',
'sigE'    :   'K03091',
'sigF'     :   'K03091',
'sigG'     :  'K03091',
'sigH'     :   'K03091',
'spoIIIC/spoIVCB'  :   'K03091',
'spoIVCA'   :  'K06400',
'spmA'      :  'K06373',
'spmB'      :  'K06374',
'obgE'      :  'K03979',
'spoIIAA'   :  'K06378',
'spoIID'    :  'K06381',
'spoIIE'    :  'K06382',
'spoIIGA'   :  'K06383',
'spoIIM'    :  'K06384',
'spoIIP'    :  'K06385',
'spoIIQ'    :  'K06386',
'spoIIR'    :  'K06387',
'spoIIIAA'   :     'K06390',
'spoIIIAB'    :    'K06391',
'spoIIIAC'     :   'K06392',
'spoIIIAD'      :  'K06393',
'spoIIIAE' :       'K06394',
'spoIIIAF'  :      'K06395',
'spoIIIAG'   :     'K06396',
'spoIIIAH'    :    'K06397',
'spoIIID' :    'K06283',
'spoIIIE'  :   'K03466',
'spoIIIJ'   :  'K03217',
'spoIVA'    :  'K06398',
'spoIVB'    :  'K06399',
'spoIVFA'   :  'K06401',
'spoIVFB'   :  'K06402',
'spoVAA'    :  'K06403',
'spoVAB'    :  'K06404',
'spoVAC'    :  'K06405',
'spoVAD'    :  'K06406',
'spoVAF'    :  'K06408',
'spoVB'     :  'K06409',
'spoVD'     :  'K08384',
'spoVE'     :  'K03588',
'spoVFA'    :  'K06410',
'spoVFB'    :  'K06411',
'spoVK'     :  'K06413',
'spoVM'     :  'K06414',
'spoVR'     :  'K06415',
'spoVT'     :  'K04769',
'bofA'      :  'K06317',
'dapA'      :  'K01714',
'etfA'      :  'K03522'}


if 'K07699' in KOs: # Master gene detected
    for gene, KO in spore_db.iteritems():
        if KO in KOs:
            sporulation.append(gene)
    if len(sporulation) > 40: # The majority of essential spore genes must be present 
        spore_genes = ','.join(sporulation)
        outputting_overview.write('Spore formation was predicted via the presence of Spo0A and the following essential sporulation proteins; ' + spore_genes + '\n')




outputting_overview.flush()



#############################################
# 
#       Annotate genome for CARD profile
#
#######################

outputting_overview.write('\n\nAntibiotic resistance analysis\n')
outputting_overview.write('------------------------------\n\n')




bashCommand = 'mkdir ' + dir_path + 'Output/'+project_name+'/CARD-annotation' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = 'diamond blastp --query-cover 80.0 --subject-cover 80.0 --id 80.0 -q ' + dir_path + 'Output/'+project_name+'/PROKKA-annotation/Isolate.faa -d ' + os.environ['PROTOLOGGER_DATA_DIR'] + '/CARD-annotation/CARDDB-DIAMOND -o ' + dir_path + 'Output/'+project_name+'/CARD-annotation/'+project_name+'.CARD.m8' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

Antibiotic_data = {}
ln = 0
for line in open(os.environ['PROTOLOGGER_DATA_DIR'] + '/CARD-annotation/aro_categories_index.tsv','r'):
    timber = line.replace('\n','').split('\t')
    ##print timber
    ln +=1
    if ln >1:
        Antibiotic_data[timber[0]] = timber[-1:][0] + ' resistance may be conferred by ' + timber[4].replace('\n','') + ' via detection of ' + timber[2]


Annotations = {}

for line in open(dir_path + 'Output/'+project_name+'/CARD-annotation/'+project_name+'.CARD.m8','r'):
    timber = line.replace('\n','').split('\t')
    query = timber[0]
    bitscore = float(timber[11])
    CARD = timber[1]
    if query in Annotations:
        if bitscore > Annotations[query][1]:
            Annotations[query] = [CARD,bitscore]
    else:
        Annotations[query] = [CARD,bitscore]

failure_ARG_not_in_metadata = 0

for k,v in Annotations.iteritems():
    try:
        outputting_overview.write(Antibiotic_data[v[0].split('|')[1]].capitalize() + '\n')
    except:
        failure_ARG_not_in_metadata +=1

if len(Annotations) == 0:
    outputting_overview.write('No antibiotic resistance genes were identified within the genome.\n')






#############################################
# 
#       Annotate genome for CAZyme profile
#
#######################

outputting_overview.write('\n\nCAZyme analysis\n')
outputting_overview.write('---------------\n\n')



bashCommand = 'mkdir ' + dir_path + 'Output/'+project_name+'/CAZy-annotation' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = 'diamond blastp --query-cover 40.0 --subject-cover 40.0 --min-score 100.0 -q ' + dir_path + 'Output/'+project_name+'/PROKKA-annotation/Isolate.faa -d ' + os.environ['PROTOLOGGER_DATA_DIR'] + '/CAZy-annotation/CAZyDB-DIAMOND -o ' + dir_path + 'Output/'+project_name+'/CAZy-annotation/'+project_name+'.CAZy.m8' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()









best_hits = {}

for line in open(dir_path + 'Output/'+project_name+'/CAZy-annotation/'+project_name+'.CAZy.m8','r'):
    timber = line.replace('\n','').split('\t')
    query = timber[0]
    hit = timber[1]
    bit = float(timber[11])
    pid = float(timber[2])
    if query in best_hits:
        prev = best_hits[query]
        if bit > prev[1]:
            best_hits[query] = [hit,bit]
    else:
        best_hits[query] = [hit,bit]





GTs = {}
GHs = {}
PLs = {}
CEs = {}
CBMs = {}

for k,v in best_hits.iteritems():
    name = v[0].split('|')
    for i in name[1:]:
        red = i.split('_')[0]
        if red.startswith('GT'):
            if red in GTs:
                GTs[red] +=1
            else:
                GTs[red] = 1
                
        elif red.startswith('GH'):
            if red in GHs:
                GHs[red] +=1
            else:
                GHs[red] = 1
                
        elif red.startswith('PL'):
            if red in PLs:
                PLs[red] +=1
            else:
                PLs[red] = 1
                
        elif red.startswith('CE'):
            if red in CEs:
                CEs[red] +=1
            else:
                CEs[red] = 1
                
        elif red.startswith('CBM'):
            if red in CBMs:
                CBMs[red] +=1
            else:
                CBMs[red] = 1




#print 'In total the following number of CAZymes were identified within the genome; ' + str(len(best_hits.keys()))
outputting_overview.write('In total the following number of CAZymes were identified within the genome; ' + str(len(best_hits.keys()))+'\n')

if len(GHs) > 0:
    #print 'The following occurances of glycoside hydrolase families were identified within the genome; '
    outputting_overview.write('The following occurances of glycoside hydrolase (GH) families were identified within the genome; '+'\n')
    for k,v in GHs.iteritems():
        #print k, ';', v
        outputting_overview.write(k+ ';'+ str(v)+'\n')
        
if len(GTs) > 0:
    #print 'The following occurances of glycoside transferase families were identified within the genome; '
    outputting_overview.write('The following occurances of glycoside transferase (GT) families were identified within the genome; '+'\n')
    for k,v in GTs.iteritems():
        #print k, ';', v
        outputting_overview.write(k+ ';'+ str(v)+'\n')
       
if len(PLs) > 0:
    #print 'The following occurances of polysaccharise lyase families were identified within the genome; '
    outputting_overview.write('The following occurances of polysaccharise lyase (PL) families were identified within the genome; '+'\n')
    for k,v in PLs.iteritems():
        #print k, ';', v
        outputting_overview.write(k+ ';'+ str(v)+'\n')
        
if len(CEs) > 0:
    #print 'The following occurances of carbohydrate esterase families were identified within the genome; '
    outputting_overview.write('The following occurances of carbohydrate esterase (CE) families were identified within the genome; '+'\n')
    for k,v in CEs.iteritems():
        #print k, ';', v
        outputting_overview.write(k+ ';'+ str(v)+'\n')
        
if len(CBMs) > 0:
    #print 'The following occurances of carbohydrate-binding module (CBM) families were identified within the genome; '
    outputting_overview.write('The following occurances of carbohydrate-binding module (CBM) families were identified within the genome; '+'\n')
    for k,v in CBMs.iteritems():
        #print k, ';', v
        outputting_overview.write(k+ ';'+ str(v)+'\n')













########################################################################################################################################################
#
#
#
#
#
#
#
#    Ecology
#
#
#
#
#
#
#
####################################################################
outputting_overview.write('\n\nEcological analysis\n')
outputting_overview.write('-------------------\n\n')




# MAG database
output_mash = open(dir_path + 'output-' + project_name + '.sh','w')
bashCommand = 'mash dist ' + os.environ['PROTOLOGGER_DATA_DIR'] + '/MAG_database/MAG_MASH_all.msh ' + genome_file + ' -t > ' + dir_path + 'Output/' + project_name + '/MAG_MASH.txt'
print bashCommand
output_mash.write(bashCommand )
output_mash.close()
new_command = 'bash '+ dir_path + 'output-' + project_name +'.sh '
process = subprocess.Popen(new_command.split(), stdout=subprocess.PIPE)
output, error = process.communicate()
bashCommand = 'rm '+ dir_path + 'output-' + project_name + '.sh '
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()



MASH_values = {}

matched_genomes = {}

ln = 0

for line in open(dir_path + 'Output/'+project_name + '/MAG_MASH.txt','r'):
    ln +=1
    timber = line.split('\t')
    if ln == 1:
        for num, splinter in enumerate(timber):
            MASH_values[num] = [splinter]
    if ln == 2:
        for num, splinter in enumerate(timber):
            if num > 0:
                if float(splinter) < 0.05:
                    #print MASH_values[num]
                    #print MASH_values[num][0].split('/')[-1:][0][0]
                    matched_genomes[MASH_values[num][0].split('/')[-1:][0]] = float(splinter)
    


MAG_output = open(dir_path + 'Output/'+project_name + '/MAG_comparison_overview.tab','w')

MAG_output.write('#MAG-ID\tHost-environment\tStudy\tSample meta-data\tMASH distance\n')

for line in open(os.environ['PROTOLOGGER_DATA_DIR'] + '/MAG_database/Meta_data_combined.txt','r'):
    timber = line.replace('\n','').split('\t')
    if timber[0] in matched_genomes.keys():
        MAG_output.write(line.replace('\n','') + '\t' + str(matched_genomes[timber[0]]) +'\n')

MAG_output.close()


if len(matched_genomes) > 0:
    outputting_overview.write('Metagenome assembled genomes (MAGs) clustering with your isolates genome were identified, in total '+ str(len(matched_genomes)) + ' MAGs were identified \n')
else:
    outputting_overview.write('No metagenome assembled genomes (MAGs) matching your genome were identified.\n')







#################################
#
# IMNGS
#
############


    
outputting_overview.write('\n\n')

bashCommand = 'mkdir  ' + dir_path + 'Output/'+project_name+'/IMNGS-analysis' 
#print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

    



for cenv in glob.glob(os.environ['PROTOLOGGER_DATA_DIR'] + '/IMNGS/FASTA-1000_files/FASTA_files/*.fasta'):
    bashCommand = 'blastn -subject ' + File_16S + '  -qcov_hsp_perc 80.0 -evalue 0.0000000000000000000000001 -perc_identity 97.0 -strand both -outfmt 6 -query '+cenv+' -out ' + dir_path + 'Output/'+project_name+'/IMNGS-analysis/' + cenv.split('/')[-1:][0].replace('.fasta','.m8') 
    print bashCommand
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()





matched_otus = {}

for cfile in glob.glob(dir_path + 'Output/'+project_name+'/IMNGS-analysis/*.m8'):
    env = cfile.split('/')[-1:][0].replace('.m8','')
    matches = []
    for line in open(cfile,'r'):
        matches.append(line.split('\t')[0].split(';')[0])
    matched_otus[env] = matches


with open(os.environ['PROTOLOGGER_DATA_DIR'] + '/IMNGS/FASTA-1000_files/Abundance_profiles/IMNGS_abundance_values.pickle', 'rb') as handle:
    abundances = pickle.load(handle)


# Prevalance

prev = {} #define which samples belong to which environment

for env,otus in matched_otus.iteritems():
    samples = []
    for i in otus:
        sample = i.split('.')[0]
        if sample not in samples:
            samples.append(sample)
    prev[env] = samples



# abundance calculation
remove_samples = []

abund = {}
for env,otus in matched_otus.iteritems():
    samples = {}
    for i in otus:
        sample = i.split('.')[0]
        try:
            abun = abundances[i]
            if sample in samples:
                samples[sample] += abun
            else:
                samples[sample] = abun
        except:
            lw = 0
            ##print sample, i
            if sample not in remove_samples:
                remove_samples.append(sample)
    abund[env] = samples



abundance = {}

for env,samples in abund.iteritems():
    total = []
    for i in samples.values():
        total.append(i)
    if len(samples.keys()) == 0:
    	abundance[env] = [0.0, 0.0]
    else:
	    abundance[env] = [np.mean(total),np.std(total)]


for k,v in prev.iteritems():
    try:
        abun = abundance[k][0]
        stdev = abundance[k][1]
        outputting_overview.write('The isolate was detected in ' + str((float(len(v))/1000)*100)  + '% of 1,000 amplicon samples from the ' + k.split('-')[0].replace('_',' ') + ' at a mean relative abundance of ' + str("{:.2f}".format(abun)) + '% with a standard deviation of ' + str("{:.2f}".format(stdev)) + '%.\n')
        print 'The isolate was detected in ' + str((float(len(v))/1000)*100)  + '% of ' + k.split('-')[0].replace('_',' ') + ' samples at a mean relative abundance of ' + str("{:.2f}".format(abun)) + '% with a standard deviation of ' + str("{:.2f}".format(stdev)) + '%.\n'
    except:
        outputting_overview.write('The isolate was detected in ' + str((float(len(v))/1000)*100)  + '% of ' + k.split('-')[0].replace('_',' ') + ' samples.\n')
        print 'The isolate was detected in ' + str((float(len(v))/1000)*100)  + '% of ' + k.split('-')[0].replace('_',' ') + ' samples.\n'




outputting_overview.flush()










outputting_overview.flush()














#print 'Analysis complete'



outputting_overview.close()
















output_overview = sys.argv[3] 
output_16tree = sys.argv[4] 
output_GenomeTree = sys.argv[5] 
output_16ID = sys.argv[6] 
output_POCP = sys.argv[7] 
output_ANI = sys.argv[8] 
output_GC = sys.argv[9] 
output_MAG = sys.argv[10]

bashCommand = 'cp '+dir_path + 'Output/' +project_name + '/Overview.txt ' + output_overview
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = 'cp '+dir_path + 'Output/' + project_name + '/16S-tree.nwk ' + output_16tree
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = 'cp '+dir_path + 'Output/' + project_name +  '/Genome_Tree.nwk ' + output_GenomeTree
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = 'cp '+dir_path + 'Output/' + project_name +  '/16S-identity-values.tab ' + output_16ID
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = 'cp '+dir_path + 'Output/' + project_name +  '/POCP_results.tab ' + output_POCP
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


outputting = open(output_ANI,'w')
for line in open(dir_path + 'Output/'  + project_name +  '/ANI_values.tab'):
    outputting.write(line.replace('/DATA/galaxy-run/galaxy/tools/protologger/Input/' +project_name + '/', '').replace('/DATA/galaxy-run/galaxy/tools/protologger/Output/' +project_name + '/', '').replace('_genomic.fna','').replace('Genome_analysis/','').replace('.fna','').replace('dataset_','Input-'))
outputting.close()


bashCommand = 'cp '+dir_path + 'Output/'  + project_name +  '/GC_of_closest_relatives.tab ' + output_GC
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

outputting = open(output_GC,'w')
for line in open(dir_path + 'Output/'  + project_name +  '/GC_of_closest_relatives.tab'):
    outputting.write(line.replace(',' , '\t'))
outputting.close()


bashCommand = 'cp '+dir_path + 'Output/'  + project_name +  '/MAG_comparison_overview.tab ' + output_MAG
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


################
#
# Not sure if this is useful? But if the files have been moved somewhere else the folder is no longer needed
#
#####
#bashCommand = 'rm -r '+ dir_path + 'Output/' + project_name
#print bashCommand
#process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
#output, error = process.communicate()

bashCommand = 'rm -r '+ dir_path + 'Input/' + project_name
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


bashCommand = 'rm -r '+ dir_path + 'Output/' + project_name
print bashCommand
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


print 'Old folder removed A-OK'




print 'Done'

exit
