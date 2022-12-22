#!/usr/bin/python
################################################################################
# retrieveKEGG
# Access the KEGG API and retrieves all data available for each protein-coding
# gene of the "n" organisms specified. Creates a file for each succesful query.

# Ivan Domenzain. Last edited: 2018-04-10
################################################################################

#INPUTS:
#1) Organism KEGG codes (as many as you want). Full list at:
#   http://rest.kegg.jp/list/organism
organism_codes = ['sce',...,...,...]
#2) Path for storing all generated files:
output_path = '.../GECKO/databases/KEGG'
#3) Last organism processed (if the program was interrupted)
#   Starting form scratch?, leave empty:
last_organism = ''
#4) Last gene entry processed (if the program was interrupted),
#   Starting form scratch?, leave empty:
last_entry = ''

################################################################################

#retrieve_org_genesData: Function that extracts all data available
#in KEGG database for the organism in turn.
def retrieve_org_genesData(organism, last_entry):
    
    #URL that returns the entire genes list for the organism
    url        = 'http://rest.kegg.jp/list/' + organism
    genes_list = []
    
    #Try/except for avoiding execution abortions ending in case of
    #querying timeout exceeded
    try:
        #Stores the queried genes list as a string
        data_str = urllib2.urlopen(url, timeout=20).read()
        
        #String division into substrings for each gene. Just the entry names are
        #saved on a list. Previously queried genes, if any, are removed from the list.
        separator  = organism + ':'
        substrings = data_str.split(separator)

        for i in substrings:
            if i[0:i.find('\t')]!=(' ' and '\0'and ''):
                genes_list.append(i[0:i.find('\t')])

        if last_entry!='':
            genes_list=genes_list[genes_list.index(last_entry):]
        
        #Retrieves gene data, if sucessfuly queried and a UniProt code is found
        #then a file .txt is created, otherwise, a warning is displayed
        for gene in genes_list:
            gene_query, gene_string = extract_geneEntry_data(organism, gene)
            
            if gene_query.find('UniProt:')!=-1:
                if gene_query!='':
                    fid = open(gene + '.txt','w')
                    fid.write(gene_query.decode('ascii','ignore'))
                    fid.close()
                    print 'Succesfully constructed ' + gene_string + '.txt'
                else:
                    print 'Unsuccesful query for gene ' + gene_string
            else:
                print 'No UniProt code for ' + gene_string
    except:
        print organism + ' not found or timeout exceeded'

################################################################################

#extract_geneEntry_data: Function that retrieves specific
#gene entries from KEGG
def extract_geneEntry_data(organism, gene):
    #URL that returns available data of the gene entry on KEGG
    gene_string = organism+ ':' + gene
    url         = 'http://rest.kegg.jp/get/' + gene_string
    
    #Try/except for avoiding timeout exceedings
    try:
        gene_query = urllib2.urlopen(url, timeout=20).read()
    except:
        gene_query=''

    return(gene_query, gene_string)

################################################################################

#Main script

#Get current path:
import os
prev_path = os.getcwd()

#Remove organisms already queried from the list
if last_organism!='':
    organism_codes=organism_codes[organism_codes.index(last_organism):]

#extensible library for opening URLs
import urllib2
#Main loop: retrieves all genes found for every organism
for organism in organism_codes:
    
    #Creates (if not present) a subfolder for the organism inside the
    #specified output path
    org_path = output_path + '/' + organism
    if not os.path.exists(org_path):
        os.makedirs(org_path)
    #access to the created organism subfolder
    os.chdir(org_path)

    #gets and creates files for all the gene entries found for the organism
    organism_genes=retrieve_org_genesData(organism, last_entry)

os.chdir(prev_path)

################################################################################


