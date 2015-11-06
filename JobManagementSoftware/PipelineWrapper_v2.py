
import PatientWrapperClass_v2
from optparse import OptionParser ## this is code to help with option parsing
import os

# example:
# python PipelineWrapper_v2 --P 601453_NA5K5CHE --O /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/TestNewPipeline_cp3 --C ./pancreas_cp3_test.csv --CFQ1 /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/TestNewPipeline_cp3/601453_NA5K5CHE.cancer.output1.fastq --CFQ2 /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/TestNewPipeline_cp3/601453_NA5K5CHE.cancer.output2.fastq --NFQ1 /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/TestNewPipeline_cp3/601453_NA5K5CHE.normal.output1.fastq --NFQ2 /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/TestNewPipeline_cp3/601453_NA5K5CHE.normal.output2.fastq --E cmelton@stanford.edu 

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
    parser.add_option("--CB", dest = "cancerBam", help = "",
                      metavar = "FILE", default = "", type = "string")
    parser.add_option("--NB", dest = "normalBam", help = "",
                      metavar = "FILE", default = "", type = "string")
    parser.add_option("--CFQ1", dest = "cancerFASTQ1", help = "",
                      metavar = "FILE", default = "", type = "string")
    parser.add_option("--CFQ2", dest = "cancerFASTQ2", help = "",
                      metavar = "FILE", default = "", type = "string")
    parser.add_option("--NFQ1", dest = "normalFASTQ1", help = "",
                      metavar = "FILE", default = "", type = "string")
    parser.add_option("--NFQ2", dest = "normalFASTQ2", help = "",
                      metavar = "FILE", default = "", type = "string")   
    parser.add_option("--E", dest = "emailAddress", help = "",
                      metavar = "STRING", default = "cmelton@stanford.edu", type = "string")
    parser.add_option("--CID", dest = "cancerID", help = "",
                      metavar = "STRING", default = "cancer", type = "string")
    parser.add_option("--NID", dest = "normalID", help = "",
                      metavar = "STRING", default = "normal", type = "string")
    parser.add_option("--NL", dest = "normalLoc", help = "",
                      metavar = "STRING", default = "", type = "string")
    parser.add_option("--RN", dest = "runNormals", help = "",
                      metavar = "STRING", default = "T", type = "string")
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

    newPatient=PatientWrapperClass_v2.Patient(options.patientID, options.outputDirectory, options.csvFile, options.emailAddress)    
    
    # ungzip and name fastqs properly
#     FASTQDirs=options.gzippedFastQDir.split(",")
#     print FASTQDirs
#     Fastqs=[]
#     if options.gzippedFastQDir!="":
#         for dirname in FASTQDirs: 
#             for x in os.listdir(dirname):
#                 if ".fastq" in x: Fastqs.append(x)
#     cid, nid =options.cancerID, options.normalID
#     print cid, nid
#     multVarNorm, multVarCancer, dep = newPatient.unzipPersonalis(FASTQDirs, gzippedFastqs, cid, nid, run=False) 
    multVarNorm = "" #"1|2"
    multVarCancer  = "" #"1|2"
    
    multVarNorm = "|".join(map(str, range(100))) #"1|2"
    multVarCancer  = "|".join(map(str, range(100))) #"1|2"
    
    print "multVarNorm:", multVarNorm
    print "multVarCancer:", multVarCancer
    # run bwa from fastqs and process bam for gatk
#     multVarNorm="1"
#     multVarCancer="1"
#     print "dependencies:", dep
    cancerFASTQ1=options.cancerFASTQ1
    cancerFASTQ2=options.cancerFASTQ2
    normalFASTQ1=options.normalFASTQ1
    normalFASTQ2=options.normalFASTQ2
#     previousNormalDep, previousCancerDep=newPatient.addMultipleFASTQs(multVarNorm, multVarCancer, isPersonalis=False, dependencies=[])
    
    
    # convert bam to fastq
    
    # convert bams to fastq
#     previousNormalDep, previousCancerDep = newPatient.bamToFastq(options.cancerBam, cancerFASTQ1, cancerFASTQ2, options.normalBam, normalFASTQ1, normalFASTQ2, dependencies=[])
    
    # convert fastq to bam
    previousNormalDep, previousCancerDep=[],[]
    previousNormalDep, previousCancerDep=newPatient.addMultipleFASTQs(multVarNorm, multVarCancer, cancerFASTQ1=cancerFASTQ1, 
                                                                      cancerFASTQ2=cancerFASTQ2, normalFASTQ1=normalFASTQ1, 
                                                                      normalFASTQ2=normalFASTQ2, isPersonalis=False, 
                                                                      dependencies=previousNormalDep+previousCancerDep, 
                                                                      runNormals=(options.runNormals=="T"), normalLoc=options.normalLoc)
    
#     # start here if data come aligned
#     
   
    # run local realignment and base recalibration
    previousNormalDep, previousCancerDep=newPatient.addRealignmentAndRecal(previousNormalDep=previousNormalDep, previousCancerDep=previousCancerDep)
           
    # find germline snvs, also run on cancer to do quality check on shared snvs
    newPatient.addGermlineSNVsAndQualityCheck(previousNormalDep=previousNormalDep, previousCancerDep=previousCancerDep)
    
          
    # find somatic snvs and indels and annotate them
    newPatient.addMutectAndVarscan2(previousNormalDep=previousNormalDep, previousCancerDep=previousCancerDep)
         
#     # find copy number alterations and annotate them
#     newPatient.addBicseq("low", previousNormalDep=previousNormalDep, previousCancerDep=previousCancerDep)
         
    # find rearrangements and other somatic structural variants
    newPatient.addCREST(previousNormalDep=previousNormalDep, previousCancerDep=previousCancerDep)
       
    # find CNVs
    newPatient.addExomeCNV(previousNormalDep=previousNormalDep, previousCancerDep=previousCancerDep)
#     newPatient.addExomeCNV(previousNormalDep=[], previousCancerDep=[])

    # generate csv file and run pipeline
    newPatient.run()
     
    
run()