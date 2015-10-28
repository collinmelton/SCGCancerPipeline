#!/usr/bin/env python

# Python code for annotating the SNVs and small indels in vcf format using Annovar
# The default gene annotation is based on refgene
# It returns a list of variant objects and output a tsv file with the variant annotations
# Only variants (SNV or indels) in or near (1kb upstream or downstream) the cancer genes are included in the .tsv file
# 
# usage: vcf_annotation.py -i <input fileName.vcf> -s <variant source>

# Created on July 14, 2014

import os
import sys
import inspect
import VariantFormatConverters
import argparse
from time import strftime, localtime

# SOME IMPORTANT GLOBAL VARIABLES
# ANNOVARPATH = '/srv/gs1/software/annovar/annovar-11122013/'
ANNOVARPATH = '/home/cmelton/cancerPatients/SharedSoftware/Annotation/annovar'
DATAPATH = '/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/CODES/Data/ForDB/'
# DATAPATH = '/Users/cmelton/Documents/AptanaStudio3Workspace/CODES/Data/ForDB/'

def loadSymbols(names):
    result = []
    for name in names:
        f = open(os.path.join(DATAPATH, "GeneAnnotations", name))
        result = result + map(lambda x: x.upper(), f.read().split("\n"))
        f.close
    return set(result)

def loadGeneSymbolLookup():
    gene_dict = dict()    
    gene_table = open(os.path.join(DATAPATH, "Genes and Gene Synonyms", 'hugo_gene_symbol_data.tsv'), 'r')
    
    for counter, line in enumerate(gene_table):
        if counter == 0: continue
        line = line.strip().split('\t')       
        hugo_symbol = line[1]
        gene_dict[hugo_symbol.upper()] = hugo_symbol
        if len(line) > 4:
            old_symbol = line[3].split(',')
            old_symbol = map(lambda x: x.strip().upper(), old_symbol)
            for x in old_symbol: 
                if x != '': gene_dict[x] = hugo_symbol
        if len(line) > 5:
            synonyms = line[4].split(',')
            synonyms = map(lambda x: x.strip().upper(), synonyms)
            for x in synonyms: 
                if x != '': gene_dict[x] = hugo_symbol 
                            
    gene_table.close() 
    return gene_dict

## MORE IMPORTANT GLOBAL VARIABLES
TS = loadSymbols(["CancerGeneCensusRec.txt", "TSGene.txt"])
ONCOGENES = loadSymbols(["CancerGeneCensusDom.txt", "PanCancerDriver.txt"])-TS
PHARMGKBGENES = loadSymbols(["PharmGKBVIP.txt"])
HUGOLOOKUP = loadGeneSymbolLookup()
        
def geneAnnotation(inputName, annovar_path):
    '''The function performs gene-based annotation. It takes a variant file in Annovar input format and 
    returns the name of the annotated file'''

    print '\nPerform gene-based annotations\n'
    
    outputName=inputName+".out"
    try:
        command_line = 'perl $path/annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene -outfile $outfile.refgene -exonsort $queryfile $dbloc'
        command = command_line.replace('$path', annovar_path).replace('$outfile', outputName).replace('$queryfile', inputName+'.avinput').replace('$dbloc', annovar_path+'/humandb')
        status = 'Return code for command '+command+':'+str(os.system(command))
        print status
        return outputName
    except IOError:
        print >> sys.stderr, 'Error with gene level annotation'    
        exit(1)

def variantAnnotation(inputName, annovar_path):  
    '''The input file contains variants in known tumor suppressor genes or oncongenes
    (including 1kb upstream and 1kb downstream). The input format follows: annovarinput format 
    with cancer gene info added in the format and sample fields (the last two columns)
    The inputName is the same as the output of the filterByCancerGenes function '''

    print '\nAnnotate the variants in cancer genes\n'
    
    try:
        command_line = 'perl $path/table_annovar.pl -buildver hg19 -protocol refGene,ljb23_all,cosmic68wgs,clinvar_20140211,segdup,popfreq_all,caddgt20 -operation g,f,f,f,r,f,f $queryfile $dbloc -nastring NA'
        command = command_line.replace('$path', annovar_path).replace('$queryfile', inputName).replace('$dbloc', annovar_path+ '/humandb')
        status = 'Return code for command '+command+':'+str(os.system(command))
        print status
    except IOError:
        print >> sys.stderr, 'Error with gene level annotation'    
        exit(1)

def get_HUGOsymbol(gene):
    if '(' in gene: gene = gene.split('(')[0] #for entries like CDKN2a(p14)
    HUGOsymbol = HUGOLOOKUP.get(gene.upper())
    if gene.upper() == 'ALO17': HUGOsymbol = 'RNF213'
    return HUGOsymbol

def is_suppressor(symbol):
    return symbol in TS

def is_oncogene(symbol):
    return symbol in ONCOGENES
    
def is_PharmGKBGene(symbol):
    return symbol in PHARMGKBGENES

def filterByCancerGenes(inputName):   
    '''The function performs cancer-gene-based filtering. The input file is the output from gene-based annotation.
    The output file contains the exonic, intronic and splicing variants, and variants in the upstream and downstream of 
    known tumor suppressor genes or oncongenes. It returns the name of the filted file'''

#     gene_dict, oncogenes, suppressors=getCancerGenes(hugoFile, cancerGeneFile)
    print '\nFilter by cancer genes and pharm_genes\n'
    
    out_file = open(inputName + '.selectedGenes.txt', 'w')
    annotated_vcf = inputName + '.refgene.variant_function'
    
    in_file = open(annotated_vcf, 'r')
    for line in in_file:
        line = line.strip().split('\t')
        if line[0] != 'intergenic':
            gene_name = ''
            #add columns for TS: tumor suppressor; OG: oncogene; PG: pharm gene; HG: gene symbol
            TS, OG, PG, HG='F', 'F', 'F', '.' 
             
            if line[0] == 'splicing': gene_name = line[1].split('(')[0]
            else: gene_name = line[1]
            
            if is_suppressor(gene_name):  TS= 'T'
            if is_oncogene(gene_name): OG= 'T'
            if is_PharmGKBGene(gene_name): PG= 'T'
            
            if TS == 'T' or OG == 'T' or PG == 'T':
                HUGO = get_HUGOsymbol(gene_name)
                if HUGO: HG = HUGO
                
                # only export variants in or near tumor suppressor genes, oncogenes, or pharm gene
                out_line = line[2:] + [TS, OG, PG, HG]
                out_line = '\t'.join(out_line) + '\n'
                out_file.write(out_line)
        
    out_file.close()
    in_file.close()
    return inputName + '.selectedGenes.txt' 
      
def vcf_annotation(inputName, source):
    annovar_path = ANNOVARPATH
    
    # reformat vcf file to be acceptable input for annovar
    VariantFormatConverters.vcfToAnnovarInput(inputName, annovar_path)
    
    # annotate with gene symbols
    outputName=geneAnnotation(inputName, annovar_path)
    
    # filter variants
    filtered_output=filterByCancerGenes(outputName)
    
    # annotate filtered variants
    variantAnnotation(filtered_output, annovar_path)     

if __name__=="__main__":
 
    #get args from command line
    parser = argparse.ArgumentParser(description = 'Wrapper for annotation')
    parser.add_argument('-i', '--input', required = True, help = 'input vcf file')
    parser.add_argument('-s', '--source', required = True, help = 'variant source: germline or somatic')
    parser.add_argument('-o', '--output1', required = True, help = 'output file 1')
    parser.add_argument('-p', '--output2', required = True, help = 'output file 2')
    args = parser.parse_args()
    source=args.source
        
    if args.input[-4:]==".vcf": 
        inputName=args.input[:len(args.input)-4]          
        vcf_annotation(inputName, source)
        # rename output files
        infile_multianno = inputName + '.out.selectedGenes.txt' + '.hg19_multianno.txt'
        infile_cancergenes = inputName + '.out.selectedGenes.txt' + '.refGene.variant_function'
        os.rename(infile_multianno, args.output1)
        os.rename(infile_cancergenes, args.output2)
        
    else:
        print >> sys.stderr, 'Error in annotate_vars.py: Invalid vcf input file format'
        exit(1)
