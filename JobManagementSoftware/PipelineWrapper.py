
import PatientWrapperClass
from optparse import OptionParser ## this is code to help with option parsing
import os

def getOptions():
    parser = OptionParser()
    parser.add_option("--P", dest = "patientID", help = "",
                      metavar = "STRING", type = "string", default = "")
    parser.add_option("--O", dest = "outputDirectory", help = "",
                      metavar = "STRING", default = "", type = "string")
    parser.add_option("--C", dest = "csvFile", help = "",
                      metavar = "FILE", default = "", type = "string")
    parser.add_option("--GZFQDIR", dest = "gzippedFastQDir", help = "",
                      metavar = "PATH", default = "", type = "string")
    parser.add_option("--CFQ1", dest = "cancerFASTQ1", help = "",
                      metavar = "FILE", default = "", type = "string")
    parser.add_option("--CFQ2", dest = "cancerFASTQ2", help = "",
                      metavar = "FILE", default = "", type = "string")
    parser.add_option("--NFQ1", dest = "normalFASTQ1", help = "",
                      metavar = "FILE", default = "", type = "string")
    parser.add_option("--NFQ2", dest = "normalFASTQ1", help = "",
                      metavar = "FILE", default = "", type = "string")   
    (options, args) = parser.parse_args()
    return options

def run():
    options=getOptions()
    if not os.path.exists(options.outputDirectory):
        os.makedirs(options.outputDirectory)
    if not os.path.exists(os.path.join(options.outputDirectory, "scripts")):    
        os.makedirs(os.path.join(options.outputDirectory, "scripts"))
    if not os.path.exists(os.path.join(options.outputDirectory, "error_and_outputs")):
        os.makedirs(os.path.join(options.outputDirectory, "error_and_outputs"))
    newPatient=PatientWrapperClass.Patient(options.patientID, options.outputDirectory, options.csvFile)    
    
    # ungzip and name fastqs properly
#     gzippedFastqs=[]
#     if options.gzippedFastQDir!="":
#         for x in os.listdir(options.gzippedFastQDir):
#             if ".fastq.gz" in x: gzippedFastqs.append(x)
#     print gzippedFastqs
#     multVarNorm, multVarCancer = newPatient.unzipPersonalis(options.gzippedFastQDir, 
#                                                             gzippedFastqs) #["AH86YJADXX_normal-121713_ATTGGCTC_L001_R1_001.fastq.gz","AH86YJADXX_normal-121713_ATTGGCTC_L001_R1_002.fastq.gz","AH86YJADXX_normal-121713_ATTGGCTC_L001_R1_003.fastq.gz","AH86YJADXX_normal-121713_ATTGGCTC_L001_R1_004.fastq.gz","AH86YJADXX_normal-121713_ATTGGCTC_L001_R1_005.fastq.gz","AH86YJADXX_normal-121713_ATTGGCTC_L001_R2_001.fastq.gz","AH86YJADXX_normal-121713_ATTGGCTC_L001_R2_002.fastq.gz","AH86YJADXX_normal-121713_ATTGGCTC_L001_R2_003.fastq.gz","AH86YJADXX_normal-121713_ATTGGCTC_L001_R2_004.fastq.gz","AH86YJADXX_normal-121713_ATTGGCTC_L001_R2_005.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R1_001.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R1_002.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R1_003.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R1_004.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R1_005.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R1_006.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R2_001.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R2_002.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R2_003.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R2_004.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R2_005.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L001_R2_006.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R1_001.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R1_002.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R1_003.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R1_004.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R1_005.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R1_006.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R2_001.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R2_002.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R2_003.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R2_004.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R2_005.fastq.gz","AH89AFADXX_normal-121713_ATTGGCTC_L002_R2_006.fastq.gz.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R1_001.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R1_002.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R1_003.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R1_004.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R1_005.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R1_006.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R1_007.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R1_008.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R2_001.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R2_002.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R2_003.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R2_004.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R2_005.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R2_006.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R2_007.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L001_R2_008.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R1_001.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R1_002.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R1_003.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R1_004.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R1_005.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R1_006.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R1_007.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R1_008.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R2_001.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R2_002.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R2_003.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R2_004.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R2_005.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R2_006.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R2_007.fastq.gz","AH89AFADXX_tumor-111513_CTAAGGTC_L002_R2_008.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L001_R1_001.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L001_R1_002.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L001_R1_003.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L001_R1_004.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L001_R1_005.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L001_R2_001.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L001_R2_002.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L001_R2_003.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L001_R2_004.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L001_R2_005.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L002_R1_001.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L002_R1_002.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L002_R1_003.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L002_R1_004.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L002_R1_005.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L002_R2_001.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L002_R2_002.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L002_R2_003.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L002_R2_004.fastq.gz","AH89AFADXX_tumor-111513_GAACAGGC_L002_R2_005.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R1_001.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R1_002.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R1_003.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R1_004.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R1_005.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R1_006.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R2_001.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R2_002.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R2_003.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R2_004.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R2_005.fastq.gz","AH86YJADXX_tumor-111513_CTAAGGTC_L001_R2_006.fastq.gz","AH86YJADXX_tumor-111513_GAACAGGC_L001_R1_001.fastq.gz","AH86YJADXX_tumor-111513_GAACAGGC_L001_R1_002.fastq.gz","AH86YJADXX_tumor-111513_GAACAGGC_L001_R1_003.fastq.gz","AH86YJADXX_tumor-111513_GAACAGGC_L001_R1_004.fastq.gz","AH86YJADXX_tumor-111513_GAACAGGC_L001_R2_001.fastq.gz","AH86YJADXX_tumor-111513_GAACAGGC_L001_R2_002.fastq.gz","AH86YJADXX_tumor-111513_GAACAGGC_L001_R2_003.fastq.gz","AH86YJADXX_tumor-111513_GAACAGGC_L001_R2_004.fastq.gz"])
#     multVarNorm, multVarCancer = newPatient.getMultVars(gzippedFastqs)
    
#     print multVarNorm
#     print multVarCancer
    # run bwa from fastqs and process bam for gatk
    multVarNorm="1"
    multVarCancer="1"
#     newPatient.addMultipleFASTQs(multVarNorm, multVarCancer, isPersonalis=False)
     
    # run local realignment and base recalibration
    newPatient.addRealignmentAndRecal()
     
    # find germline snvs, also run on cancer to do quality check on shared snvs
#     newPatient.addGermlineSNVsAndQualityCheck() 
     
    # find somatic snvs and indels and annotate them
    newPatient.addMutectAndVarscan2()
     
    # find copy number alterations and annotate them
#     newPatient.addBicseq("low")
     
    # find rearrangements and other somatic structural variants
#     newPatient.addCREST()
    
    # find CNVs
#     newPatient.addExomeCNV()
    
    # generate csv file and run pipeline
    newPatient.run()
     
    
run()