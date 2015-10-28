'''
Created on Jun 25, 2014
@author: cmelton

Edited on July 12, 2014 by Wei
'''

import math
import os
import sys

def ExomeCNVToAnnovarInput(inputName): #removed the outputName Argument, July 15, 2014
    f=open(inputName, 'r')
    outputName=inputName+".avinput" #added the outputName , July 15, 2014
    g=open(outputName, 'w')
    header=f.readline().strip().split("\t")
    line=f.readline()
    firstLine=True
    while line!="":
        vals=line.strip().split("\t")
        chrom=vals[0].replace("chr", "")
        start, end=vals[1], vals[2]
        additionalInfo=", ".join(map(lambda i: header[i]+":"+vals[i], range(0, len(vals))))
        if not firstLine: g.write("\n")
        g.write("\t".join([chrom, start, end, "0", "0", "comments: VariantType:CNV, "+additionalInfo]))
        firstLine=False
        line=f.readline()
    f.close()
    g.close()

# modified the code for the right margin of deletions
def CrestToAnnovarInput(inputName): #removed the outputName Argument, July 15, 2014
    f=open(inputName, 'r')
    outputName=inputName+".avinput" #added the outputName , July 15, 2014
    g=open(outputName, 'w')
    header=["left_chr","left_pos","left_strand","#_of_left_soft-clipped_reads","right_chr","right_pos","right_strand","#_right_soft-clipped_reads",
            "SV_type","coverage_at_left_pos","coverage_at_right_pos","assembled_length_at_left_pos","assembled_length_at_right_pos",
            "average_percent_identity_at_left_pos","percent_of_non-unique_mapping_reads_at_left_pos","average_percent_identity_at_right_pos",
            "percent_of_non-unique_mappingreads_at_right_pos","start_position_of_consensus_mapping_to_genome","starting_chromosome_of_consensus_mapping",
            "position_of_the_genomic_mapping_ofconsensus_starting_position","end_position_of_consensus_mapping_to_genome",
            "ending_chromsome_of_consnesus_mapping","position_of_genomic_mapping_ofconsensus_ending_posiiton","consensus_sequences"]
    for counter, line in enumerate(f):
        vals=line.strip().split("\t")
        varDict={}
        for i in range(0,len(vals)):
            varDict[header[i]]=vals[i]
        additionalInfo=", ".join(map(lambda i: header[i]+":"+vals[i], range(0, len(vals))))
        variant_type=varDict["SV_type"]
        chrom, start, end=varDict["left_chr"], varDict["left_pos"], varDict["left_pos"]
        
        if varDict["left_chr"] == "hs37d5" or varDict["right_chr"] == "hs37d5": continue
        if variant_type=="DEL": end=varDict["right_pos"] #for deletions, use "right_pos" as the right margin
        
        if chrom in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]:    
            g.write("\t".join([chrom, start, end, "0", "0", variant_type, "comments: Row:"+str(counter)+", "+additionalInfo])+"\n")
        
        if variant_type=="DEL": continue
        
        chrom, start, end=varDict["right_chr"], varDict["right_pos"], varDict["right_pos"]
        if chrom in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]:
            g.write("\t".join([chrom, start, end, "0", "0", variant_type, "comments: Row:"+str(counter)+", "+additionalInfo])+"\n")
    f.close()
    g.close()

def CrestInsertions(inputName): 
    '''Sort the insertions and try to find two closesest breakpoints, 
    then check and see if they are the two ends of the insertion (or tandem duplication in CREST)'''
    from operator import itemgetter
    f=open(inputName, 'r')
    outputName=inputName+".ins" 
    g=open(outputName, 'w')
    header=["left_chr","left_pos","left_strand","#_of_left_soft-clipped_reads","right_chr","right_pos","right_strand","#_right_soft-clipped_reads",
            "SV_type","coverage_at_left_pos","coverage_at_right_pos","assembled_length_at_left_pos","assembled_length_at_right_pos","average_percent_identity_at_left_pos",
            "percent_of_non-unique_mapping_reads_at_left_pos","average_percent_identity_at_right_pos","percent_of_non-unique_mappingreads_at_right_pos",
            "start_position_of_consensus_mapping_to_genome","starting_chromosome_of_consensus_mapping","position_of_the_genomic_mapping_ofconsensus_starting_position",
            "end_position_of_consensus_mapping_to_genome","ending_chromsome_of_consnesus_mapping","position_of_genomic_mapping_ofconsensus_ending_posiiton","consensus_sequences"]
    varList = []
    varDict={}
    for line in f:
        vals=line.strip().split("\t")       
        if vals[header.index('SV_type')]=='INS':         
            chrom, start=vals[header.index("left_chr")], int(vals[header.index("left_pos")])           
            additionalInfo=", ".join(map(lambda i: header[i]+":"+vals[i], range(0, len(vals))))
            var_ind = chrom + '_' + str(start)
            varDict[var_ind] = additionalInfo
            varList.append((chrom, start))

    sorted_vars = sorted(varList, key=itemgetter(0,1))    
    for item in sorted_vars:
        var_ind= item[0] + '_' + str(item[1])
        print var_ind
        additional_info = varDict.get(var_ind) 
        g.write(item[0] + '\t' + str(item[1]) + '\t' + additional_info + '\n') 
              
    f.close()
    g.close()
    
def VarscanSNPToAnnovarInput(inputName):
    ''' This function converts an SNP file generated by Varscan to two Annovar input files, one for 
    the germline variants and one for somatic variants. '''
    f=open(inputName, 'r')
    outputGermline=inputName+'.avinput.germline'
    outputSomatic=inputName+'.avinput.somatic'
    g_germline=open(outputGermline, 'w')
    g_somatic=open(outputSomatic, 'w')
    header=["CHROM", "POS", "REF", "ALT", "RD_normal", "AD_normal", "FREQ_normal", "GT_normal", "RD_tumor", "AD_tumor", "FREQ_tumor", "GT_tumor", "SS", "GPV", "SPV", "RD_p_tumor", "RD_n_tumor", "AD_p_tumor", "AD_n_tumor", "RD_p_normal", "RD_n_normal", "AD_p_normal", "AD_n_normal"]
    for counter, line in enumerate(f):
        if counter== 0: continue #skip the header line in the input file  
        vals=line.strip().split("\t")
        chrom=vals[0]
        start, end, ref, alt=vals[1], vals[1], vals[2], vals[3] 
        if vals[12]=="Germline":
            GT=vals[7]
            zygocity="homo" if GT==vals[3] else "het"
            Pvalue=float(vals[13]) if vals[13] != '0.0' else 0.00000001
            QUAL=str(int(round(-10*math.log(Pvalue, 10),0)))
            DP=str(int(vals[4])+int(vals[5])) #need to check if this is the correct calculation of DP
            infoList=range(4,8)+[13]+range(19,23)
            additionalInfo=", ".join(map(lambda i: header[i]+":"+vals[i], infoList))
            g_germline.write("\t".join([chrom, start, end, ref, alt, zygocity, QUAL, DP, "comments: VariantType:SNV, "+additionalInfo])+'\n')
        #how to deal with the LOHs?
        elif vals[12]=="Somatic":
            GT=vals[11]
            zygocity="homo" if GT==[3] else "het"
            Pvalue=float(vals[14]) if vals[14] != '0.0' else 0.00000001
            QUAL=str(int(round(-10*math.log(Pvalue, 10),0)))
            DP=str(int(vals[8])+int(vals[9])) #need to check if this is the correct calculation of DP
            infoList=range(8,12)+range(14,19)
            additionalInfo=", ".join(map(lambda i: header[i]+":"+vals[i], infoList))
            g_somatic.write("\t".join([chrom, start, end, ref, alt, zygocity, QUAL, DP, "comments: VariantType:SNV, "+additionalInfo])+'\n')
        else: continue
    f.close()
    g_germline.close()
    g_somatic.close()
    
def VarscanINDELToAnnovarInput(inputName):
    ''' This function converts an INDEL file generated by Varscan to two Annovar input files, one for 
    the germline variants and one for somatic variants. '''
    f=open(inputName, 'r')
    outputGermline=inputName+'.avinput.germline'
    outputSomatic=inputName+'.avinput.somatic'
    g_germline=open(outputGermline, 'w')
    g_somatic=open(outputSomatic, 'w')
    header=["CHROM", "POS", "REF", "ALT", "RD_normal", "AD_normal", "FREQ_normal", "GT_normal", "RD_tumor", "AD_tumor", "FREQ_tumor", "GT_tumor", "SS", "GPV", "SPV", "RD_p_tumor", "RD_n_tumor", "AD_p_tumor", "AD_n_tumor", "RD_p_normal", "RD_n_normal", "AD_p_normal", "AD_n_normal"]
    for counter, line in enumerate(f):
        if counter==0: continue #skip the header line in the input file  
        vals=line.strip().split("\t")
        chrom=vals[0]
        start=str(int(vals[1])+1)
        end=str(int(vals[1])+len(vals[3])-1) 
        ref, alt=vals[2], vals[3]
        variant_type=alt[0]
        if variant_type=="-": ref, alt, var_type=alt[1:], "-", "deletion" #change the ref and alt format to be consistent with Annovar input format
        else: ref, alt, var_type='-', alt[1:], "insertion"
        
        if vals[12]=="Germline":
            GT=vals[7].split("/")
            zygocity="homo" if GT[0]==GT[1] else "het"
            Pvalue=float(vals[13]) if vals[13] != '0.0' else 0.00000001
            QUAL=str(int(round(-10*math.log(Pvalue, 10),0)))
            DP=str(int(vals[4])+int(vals[5])) #need to check if this is the correct calculation of DP
            infoList=range(4,8)+[13]+range(19,23)
            additionalInfo=", ".join(map(lambda i: header[i]+":"+vals[i], infoList))
            g_germline.write("\t".join([chrom, start, end, ref, alt, zygocity, QUAL, DP, "comments: VariantType:"+var_type+", "+additionalInfo])+'\n')
        #how to deal with the LOHs?
        elif vals[12]=="Somatic":
            GT=vals[11].split("/")
            zygocity="homo" if GT[0]==GT[1] else "het"
            Pvalue=float(vals[14]) if vals[14] != '0.0' else 0.00000001
            QUAL=str(int(round(-10*math.log(Pvalue, 10),0))) 
            DP=str(int(vals[8])+int(vals[9])) #need to check if this is the correct calculation of DP
            infoList=range(8,12)+range(14,19)
            additionalInfo=", ".join(map(lambda i: header[i]+":"+vals[i], infoList))
            g_somatic.write("\t".join([chrom, start, end, ref, alt, zygocity, QUAL, DP, "comments: VariantType:"+var_type+", "+additionalInfo])+'\n')
        else: continue
    f.close()
    g_germline.close()
    g_somatic.close()

def vcfToAnnovarInput(inputName, annovar_path): 
    ''' This function converts a vcf file to an Annovar input file'''
    if inputName[-4:] == '.vcf':
        inputName=inputName[:-4]
    os.system('perl ' + annovar_path + '/convert2annovar.pl '+inputName+'.vcf --includeinfo --withzyg -format vcf4 -outfile '+inputName+'.avinput')
             
# ExomeCNVToAnnovarInput("/Users/cmelton/Documents/Aptana Studio 3 Workspace/CODES/TestData/CNV_Variant_Files/cp2.exomecnv.output")
#CrestToAnnovarInput("/Users/cmelton/Documents/Aptana Studio 3 Workspace/CODES/TestData/Structural_Variant_Files/cp2.crest.output")
# ExomeCNVToAnnovarInput("/Users/wwang/Documents/BMI212/data/TestData/CNV_Variant_Files/cp2.exomecnv.output")
#VarscanSNPToAnnovarInput("/Users/wwang/Documents/Aptana Studio 3 Workspace/Cancer/TestData/Somatic_SNVs_Indels/cp2.varscan.snp")
#VarscanINDELToAnnovarInput("/Users/wwang/Documents/Aptana Studio 3 Workspace/Cancer/TestData/Somatic_SNVs_Indels/cp2.varscan.indel")

#CrestInsertions("/Users/wwang/Documents/BMI212/data/TestData/Structural_Variant_Files/cp2.crest.output")
#CrestToAnnovarInput("/Users/wwang/Documents/BMI212/data/TestData/Structural_Variant_Files/cp2.crest.output")

# annovar_path="/Users/wwang/Documents/SVs/annovar"
# inputName="/Users/wwang/Documents/Aptana\ Studio\ 3\ Workspace/Cancer/TestData/Germline_Variant_Files/cp2.UnifiedGenotyper.truncated.vcf"
# vcfToAnnovarInput(inputName, annovar_path)


