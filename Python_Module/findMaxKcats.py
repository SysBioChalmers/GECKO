#!/usr/bin/python
################################################################################
# findMaxKcats
# Reads all EC files and finds the max Kcat for each substrate for the chosen
# microorganism. Writes a table with the following columns:
# 1)  EC number
# 2)  Substrate
# 3)  Max Kcat value for the organism [1/s]
# 4)  Kcat std for the organism
# 5)  Max Kcat value for the rest of the organisms [1/s]
# 6)  Kcat std for the rest of the organisms
#
# Benjamin Sanchez. Last edited: 2015-08-26
################################################################################

#INPUTS:
#1) Organism
organism = 'saccharomyces cerevisiae'
#2) Path in which the EC files are stored:
input_path = 'D:/Users/bensan.NET/Documents/BRENDA/EC_files'
#3) Path in which you wish to store the final table:
output_path = 'D:/Users/bensan.NET/Documents/Dropbox/PhD/Modelos/FBAwMC/GECKO'

################################################################################

#sub_max_std: Recieves a list of substrates///values, returns 3 lists:
#substrates - max - std
def sub_max_std(data):
    #Sort list, add a last empty line and initialize variables:
    data.sort()
    data.append('')
    previous   = ''
    substrates = []
    maxs       = []
    stdevs     = []
    
    for row in data:
        if row != '':
            substrate = row[0:row.find('///')]
            value     = float(row[row.find('///')+3:len(row)])
        
        else:
            substrate = ''
        
        #Change of substrate:
        if previous != substrate:
            if previous != '':
                maxv     = max(values)
                stdev    = numpy.std(values)
                
                substrates.append(previous)
                maxs.append(maxv)
                stdevs.append(stdev)
              
            values   = []
            previous = substrate
        
        if row != '':
            values.append(value)
        
    return (substrates,maxs,stdevs)

################################################################################

#find_in_list: Finds the position of a given substrate in a list and returns
#the average and standard deviation
def find_in_list(sub,subs,maxs,stdevs):
    try:
        pos   = subs.index(sub)
        maxv  = maxs[pos]
        stdev = stdevs[pos]
    
    except:
        maxv  = 0
        stdev = 0
    
    return (maxv,stdev)

################################################################################

#Main Script

#Read all EC file names:
import os
prev_path = os.getcwd()
os.chdir(input_path)
dir_files = os.listdir(input_path)

#Main loop:
output = ''
for ec in dir_files:
    ec_number  = ec[0:len(ec)-4]

    #Read all data in EC file and place it in two different lists:
    import csv
    fid       = open(ec,'r')
    csv_fid   = csv.reader(fid,delimiter='\t')
    kcat_org  = []
    kcat_rest = []
    
    for row in csv_fid:
        #Search in the commentary if entry is mutant:
        row[4] = row[4].lower()
        mutant = max(row[4].find('mutant'),row[4].find('mutated'))
        
        #Ignore invalid values:
        if row[2] != '-999' and mutant == -1:
            
            #Only allow Kcats <= 1e7 [Bar-Even et al. 2011]
            if row[0] == 'KCAT' and float(row[2]) <= 1e7:
                
                if row[1].lower() == organism:
                    kcat_org.append(row[3].lower() + '///' + row[2])
                
                else:
                    kcat_rest.append(row[3].lower() + '///' + row[2])
    
    fid.close()
    
    #Generate two lists: one for the organism and one for the rest.
    import numpy
    subs_kcat_org,maxs_kcat_org,stdevs_kcat_org    = sub_max_std(kcat_org)
    subs_kcat_rest,maxs_kcat_rest,stdevs_kcat_rest = sub_max_std(kcat_rest)
    
    #Form sorted merged list with all substrates (no repetitions):
    subs_all = subs_kcat_org + list(set(subs_kcat_rest)-set(subs_kcat_org))
    subs_all.sort()
    
    for sub in subs_all:            
        max_kcat_org,stdev_kcat_org   = find_in_list(sub,subs_kcat_org,maxs_kcat_org,stdevs_kcat_org)
        max_kcat_rest,stdev_kcat_rest = find_in_list(sub,subs_kcat_rest,maxs_kcat_rest,stdevs_kcat_rest)
        
        output = output + ec_number + '\t' + sub      + '\t' + str(max_kcat_org) 
        output = output + '\t' + str(stdev_kcat_org)  + '\t' + str(max_kcat_rest)
        output = output + '\t' + str(stdev_kcat_rest) + '\n'
        
    print 'Processed file ' + ec
    
#Write output:
os.chdir(output_path)
fid = open(organism + '_max_kCATs.txt','w')
fid.write(output)
fid.close()
os.chdir(prev_path)

################################################################################