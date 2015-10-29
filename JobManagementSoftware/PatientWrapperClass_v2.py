import os, re
import csv
import subprocess

class Patient():
    def __init__(self, patientID, outputDirectory, csvFile, emailAddress):
        self.patientID=patientID
        self.outputDirectory=outputDirectory
        self.csvFile=csvFile
        self.csvFileHandle=open(self.csvFile, 'w')
        self.csvWriter=csv.writer(self.csvFileHandle, delimiter=",")
        self.nextJobNumber=1
        self.csvWriter.writerow(["run","notes","patientID","scriptName","outputPath","multiplicity","scriptPath","scriptTime",
                                 "scriptOutputFileDirectory","scriptErrorFileDirectory","scriptCustomizations","scriptMemory","scriptEmailAddress","scriptCommand","inputs", "dependencies"])
        self.emailAddress=emailAddress
    
    def getMultVars(self, filenames):
        multVarNorm=set([])
        multVarCancer=set([])
        for filename in filenames:
            sampletype="normal"
            if "cancer" in filename or "tumor" in filename or "Tumor" in filename: sampletype="cancer"
            R="R1"
            if "R1" not in filename: R="R2"
            multvar=filename.replace(".fastq", "").replace(".gz", "").replace("_"+R, "").replace("."+R, "")
            if sampletype=="cancer": multVarCancer.add(multvar)
            if sampletype=="normal": multVarNorm.add(multvar)
        return "|".join(list(multVarNorm)), "|".join(list(multVarCancer))
    
    # writes a job to csv file
    def writeJob(self, scriptName, scriptTime, scriptMemory, scriptCommand, inputs, run="TRUE",
                 notes="", outputPath="", multiplicity="", scriptPath="", scriptOutputFileDirectory="",
                 scriptErrorFileDirectory="", scriptCustomizations="#", scriptEmailAddress="", dependencies=[]):
        if outputPath=="": outputPath=self.outputDirectory
        if scriptPath=="": scriptPath=self.outputDirectory+"/scripts/"#+scriptName+".sh"
        if scriptOutputFileDirectory=="": scriptOutputFileDirectory= self.outputDirectory+"/error_and_outputs"
        if scriptErrorFileDirectory=="": scriptErrorFileDirectory= self.outputDirectory+"/error_and_outputs"
        if scriptEmailAddress=="": scriptEmailAddress=self.emailAddress
        self.csvWriter.writerow([run,notes, self.patientID, scriptName, outputPath, multiplicity,
                                 scriptPath, scriptTime, scriptOutputFileDirectory, scriptErrorFileDirectory,
                                 scriptCustomizations, scriptMemory, scriptEmailAddress, scriptCommand,
                                 "|".join(inputs), "|".join(dependencies)])
    
    # unzip fastq files
    def unzipPersonalis(self, zippedDirs, filenames, cancerID, normalID, run=True):
        multVarNorm=set([])
        multVarCancer=set([])
        cancerpattern=re.compile(cancerID)
        normalpattern=re.compile(normalID)
        scriptNames=[]
        for filename in filenames:
            sampletype="none"
            if cancerpattern.search(filename)!=None: sampletype="cancer"
            elif normalpattern.search(filename)!=None: sampletype="normal"
            R="R1"
            if "R1" not in filename: R="R2"
            multvar=filename.strip(".fastq.gz").replace("_"+R, "")
            scriptName="gunzip"+"_"+R+"_"+multvar+"_"+sampletype
            scriptNames.append(scriptName)
            for zd in zippedDirs:
                filepath=os.path.join(zd, filename)
                if os.path.exists(filepath):
                    if run: self.writeJob(scriptName, "6:00:00", "4", "gunzip -c $1 > $2", [filepath,"$OUTPUTPATH/$PATIENTID_"+R+"_"+multvar+"."+sampletype+".fastq"])
            if sampletype=="cancer": multVarCancer.add(multvar)
            if sampletype=="normal": multVarNorm.add(multvar)
        self.nextJobNumber+=1
        if run: return "|".join(list(multVarNorm)), "|".join(list(multVarCancer)), scriptNames
        else: return "|".join(list(multVarNorm)), "|".join(list(multVarCancer)), []
    
    # convert bam to fastq
    def bamToFastq(self, cancerBam, cancerFastqR1, cancerFastqR2, normalBam, normalFastqR1, normalFastqR2, dependencies=[]):
        self.writeJob("tofastq_cancer", "150:00:00", "12", "JAVAPATH -Xmx8g -XX:-UseGCOverheadLimit -jar /srv/gsfs0/software/picard-tools/1.92/SamToFastq.jar INPUT=$1 FASTQ=$2 SECOND_END_FASTQ=$3", 
                      [cancerBam, cancerFastqR1, cancerFastqR2], dependencies=dependencies)
        self.writeJob("tofastq_normal", "150:00:00", "12", "JAVAPATH -Xmx8g -XX:-UseGCOverheadLimit -jar /srv/gsfs0/software/picard-tools/1.92/SamToFastq.jar INPUT=$1 FASTQ=$2 SECOND_END_FASTQ=$3", 
                      [normalBam, normalFastqR1, normalFastqR2], dependencies=dependencies)
        # return names of last jobs to be run
        return ["tofastq_cancer"], ["tofastq_normal"]
     
    # align fastqs using bwa
    def addMultipleFASTQs(self, multVarNorm, multVarCancer, 
                          cancerFASTQ1="$OUTPUTPATH/$PATIENTID_R1_$MULTIPLICITYVAR_FORFILE.cancer.fastq", 
                          cancerFASTQ2="$OUTPUTPATH/$PATIENTID_R2_$MULTIPLICITYVAR_FORFILE.cancer.fastq", 
                          normalFASTQ1="$OUTPUTPATH/$PATIENTID_R1_$MULTIPLICITYVAR_FORFILE.normal.fastq", 
                          normalFASTQ2="$OUTPUTPATH/$PATIENTID_R2_$MULTIPLICITYVAR_FORFILE.normal.fastq", 
                          isPersonalis=False, dependencies=[], runNormals=True, normalLoc=""):
        
        newdependencies = []
        if runNormals: fastqs = [cancerFASTQ1, cancerFASTQ2, normalFASTQ1, normalFASTQ2]
        else: fastqs = [cancerFASTQ1, cancerFASTQ2]
        for fastq in fastqs:
            print fastq
            if fastq[-3:]==".gz": # or ".gz" in fastq:
                self.writeJob("unzip_"+fastq.split("/")[-1].split(".")[0], "6:00:00", "4", "gunzip -c $1 > $2", [fastq, fastq[:-3]], dependencies=dependencies)
                newdependencies.append("unzip_"+fastq.split("/")[-1].split(".")[0])
        for fq in [cancerFASTQ1, cancerFASTQ2, normalFASTQ1, normalFASTQ2]:
            if fq[-3:]==".gz": fq = fq[:-3]
        # split fastqs
        if runNormals:
#             print ["$OUTPUTPATH/$PATIENTID.normal", multVarNorm, normalFASTQ1, normalFASTQ2]
            self.writeJob("split_normal_fastq", "48:00:00", "2", "python /srv/gsfs0/clinical/cancerPatientAnno/SCGCancerPipeline/OtherScripts/SplitFastqs.py --P $1 --M "+multVarNorm+" --FQ1 $2 --FQ2 $3",
                          ["$OUTPUTPATH/$PATIENTID.normal", normalFASTQ1, normalFASTQ2], 
                          dependencies=newdependencies, multiplicity="")
#         print ["$OUTPUTPATH/$PATIENTID.cancer", multVarCancer, cancerFASTQ1, cancerFASTQ2]
        self.writeJob("split_cancer_fastq", "48:00:00", "2", "python /srv/gsfs0/clinical/cancerPatientAnno/SCGCancerPipeline/OtherScripts/SplitFastqs.py --P $1 --M "+multVarCancer+" --FQ1 $2 --FQ2 $3",
                      ["$OUTPUTPATH/$PATIENTID.cancer", cancerFASTQ1, cancerFASTQ2], 
                      dependencies=newdependencies, multiplicity="")
        return ["split_cancer_fastq"], ["split_normal_fastq"]
        # align
        if runNormals:
            self.writeJob("bwa_normal", "150:00:00", "6", "BWAPATH mem -M REFERENCEPATH $1 $2 -t $3 | SAMTOOLSPATH view -Sbt REFERENCEINDEX -o $4 -",
                          ["$OUTPUTPATH/$PATIENTID.normal_$MULTIPLICITYVAR_FORFILE_R1.fastq", "$OUTPUTPATH/$PATIENTID.normal_$MULTIPLICITYVAR_FORFILE_R2.fastq", "1", "$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.bam"], 
                          dependencies=["split_normal_fastq"], multiplicity=multVarNorm)
        self.writeJob("bwa_cancer", "150:00:00", "6", "BWAPATH mem -M REFERENCEPATH $1 $2 -t $3 | SAMTOOLSPATH view -Sbt REFERENCEINDEX -o $4 -",
                      ["$OUTPUTPATH/$PATIENTID.cancer_$MULTIPLICITYVAR_FORFILE_R1.fastq", "$OUTPUTPATH/$PATIENTID.cancer_$MULTIPLICITYVAR_FORFILE_R2.fastq", "1", "$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.bam"], 
                      dependencies=["split_cancer_fastq"], multiplicity=multVarCancer)
        if isPersonalis:
            if runNormals:
                normalScriptNames=[]
                for m in multVarNorm.split("|"):
                    vals=m.split("_")
                    rgid, rglb, rgpl, rgpu, rgsm= vals[0], vals[2]+"_"+vals[3], "illumina", vals[2], vals[1]  # id, library (barcode?), platform, platform unit (barcode),   
                    scriptName="addReadGroup"+"_normal_"+m
                    normalScriptNames.append(scriptName)
                    self.writeJob(scriptName, "100:00:00", "8", "JAVAPATH -Xmx2g -jar PICARDPATH/AddOrReplaceReadGroups.jar INPUT=$1 OUTPUT=$2 RGID=$3 RGLB=$4 RGPL=$5 RGPU=$6 RGSM=$7 VALIDATION_STRINGENCY=LENIENT",
                          ["$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.bam","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.bam",rgid,rglb,rgpl,rgpu,rgsm], 
                          dependencies=["bwa_normal"], multiplicity=m)
                     
            cancerScriptNames=[]
            for m in multVarCancer.split("|"):
                vals=m.split("_")
                rgid, rglb, rgpl, rgpu, rgsm= vals[0], vals[2]+"_"+vals[3], "illumina", vals[2], vals[1]  # id, library (barcode?), platform, platform unit (barcode),  
                scriptName="addReadGroup"+"_normal_"+m
                cancerScriptNames.append(scriptName)
                self.writeJob(scriptName, "100:00:00", "8", "JAVAPATH -Xmx2g -jar PICARDPATH/AddOrReplaceReadGroups.jar INPUT=$1 OUTPUT=$2 RGID=$3 RGLB=$4 RGPL=$5 RGPU=$6 RGSM=$7 VALIDATION_STRINGENCY=LENIENT",
                      ["$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.bam","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.bam",rgid,rglb,rgpl,rgpu,rgsm], 
                      dependencies=["bwa_cancer"], multiplicity=m)
        else:
            normalScriptNames=["addReadGroup_normal"]
            cancerScriptNames=["addReadGroup_cancer"]
            if runNormals:
                self.writeJob("addReadGroup_normal", "100:00:00", "8", "JAVAPATH -Xmx2g -jar PICARDPATH/AddOrReplaceReadGroups.jar INPUT=$1 OUTPUT=$2 RGID=$3 RGLB=$4 RGPL=$5 RGPU=$6 RGSM=$7 VALIDATION_STRINGENCY=LENIENT",
                          ["$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.bam","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.bam","$PATIENTID","$MULTIPLICITYVAR_FORFILE", "Illumina","$PATIENTID","normal"], 
                          dependencies=["bwa_normal"], multiplicity=multVarNorm)
            self.writeJob("addReadGroup_cancer", "100:00:00", "8", "JAVAPATH -Xmx2g -jar PICARDPATH/AddOrReplaceReadGroups.jar INPUT=$1 OUTPUT=$2 RGID=$3 RGLB=$4 RGPL=$5 RGPU=$6 RGSM=$7 VALIDATION_STRINGENCY=LENIENT",
                      ["$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.bam","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.bam","$PATIENTID|$MULTIPLICITYVAR_FORFILE","Illumina","$PATIENTID","cancer"], 
                      dependencies=["bwa_cancer"], multiplicity=multVarCancer)
        
        # sort normal and cancer
        if runNormals:
            self.writeJob("sort_normal_1", "48:00:00", "14", "JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate",
                          ["$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.bam","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.sorted.bam"], 
                          dependencies=normalScriptNames, multiplicity=multVarNorm)
        self.writeJob("sort_cancer_1", "48:00:00", "14", "JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate",
                      ["$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.bam","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.sorted.bam"], 
                      dependencies=cancerScriptNames, multiplicity=multVarCancer)
         
        # index normal and cancer
        if runNormals:
            self.writeJob("index_normal_1", "48:00:00", "14", "JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT",
                          ["$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.sorted.bam"], 
                          dependencies=["sort_normal_1"], multiplicity=multVarNorm)
        self.writeJob("index_cancer_1", "48:00:00", "14", "JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT",
                      ["$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.sorted.bam"], 
                      dependencies=["sort_cancer_1"], multiplicity=multVarCancer)
                 
        # merge
        if runNormals:
            self.writeJob("merge_normal_1", "48:00:00", "12", 
                          "JAVAPATH -Xmx2g -jar PICARDPATH/MergeSamFiles.jar "+" ".join(map(lambda x: "INPUT=${"+str(x+1)+"}", range(len(multVarNorm.split("|")))))+" OUTPUT=${"+str(len(multVarNorm.split("|"))+1)+"} VALIDATION_STRINGENCY=LENIENT USE_THREADING=true  ",
                          ("|".join(map(lambda x: "$OUTPUTPATH/$PATIENTID."+x+".normal.RG.sorted.bam", multVarNorm.split("|")))+"|$OUTPUTPATH/$PATIENTID.merged.normal.RG.bam").split("|"), 
                          dependencies=["index_normal_1"])
        self.writeJob("merge_cancer_1", "48:00:00", "12", 
                      "JAVAPATH -Xmx2g -jar PICARDPATH/MergeSamFiles.jar "+" ".join(map(lambda x: "INPUT=${"+str(x+1)+"}", range(len(multVarCancer.split("|")))))+" OUTPUT=${"+str(len(multVarCancer.split("|"))+1)+"} VALIDATION_STRINGENCY=LENIENT USE_THREADING=true  ",
                      ("|".join(map(lambda x: "$OUTPUTPATH/$PATIENTID."+x+".cancer.RG.sorted.bam", multVarCancer.split("|")))+"|$OUTPUTPATH/$PATIENTID.merged.cancer.RG.bam").split("|"), 
                      dependencies=["index_cancer_1"])
         
        # sort normal and cancer
        if runNormals:
            self.writeJob("sort_normal_2", "48:00:00", "14", 
                          "JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate",
                          ["$OUTPUTPATH/$PATIENTID.merged.normal.RG.bam","$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.bam"], 
                          dependencies=["merge_normal_1"])
        self.writeJob("sort_cancer_2", "48:00:00", "14", 
                      "JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate",
                      ["$OUTPUTPATH/$PATIENTID.merged.cancer.RG.bam","$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.bam"], 
                      dependencies=["merge_cancer_1"])
         
        # index normal and cancer
        if runNormals:
            self.writeJob("index_normal_2", "48:00:00", "14", "JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT",
                          ["$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.bam"], 
                          dependencies=["sort_normal_2"])
        self.writeJob("index_cancer_2", "48:00:00", "14", "JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT",
                      ["$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.bam"], 
                    dependencies=["sort_cancer_2"])
        
        # mark duplicates for cancer and normal
        if runNormals:
            self.writeJob("mark_dup_normal_1", "150:00:00", "32", 
                          "JAVAPATH -Xmx24g -jar PICARDPATH/MarkDuplicates.jar INPUT=$1 OUTPUT=$2 METRICS_FILE=$3 VALIDATION_STRINGENCY=LENIENT",
                          ["$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.bam","$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.dedup.bam","$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.dedup.metrics.txt"], 
                          dependencies=["index_normal_2"])
        self.writeJob("mark_dup_cancer_1", "150:00:00", "32", 
                      "JAVAPATH -Xmx24g -jar PICARDPATH/MarkDuplicates.jar INPUT=$1 OUTPUT=$2 METRICS_FILE=$3 VALIDATION_STRINGENCY=LENIENT",
                      ["$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.bam","$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.bam","$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.metrics.txt"], 
                      dependencies=["index_cancer_2"])        

        # if not run normals tranfer normal file over
        if not runNormals and normalLoc!="":
            self.writeJob("copyNormal", "12:00:00", "2", "cp $1 > $2",
                      [normalLoc, "$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.bam"], 
                      dependencies=[])
        
        # index normal and cancer
        self.writeJob("index_normal_3", "48:00:00", "14", "JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT",
                      ["$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.dedup.bam"], 
                      dependencies=["mark_dup_normal_1"])
        self.writeJob("index_cancer_3", "48:00:00", "14", "JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT",
                      ["$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.bam"], 
                      dependencies=["mark_dup_cancer_1"])
        
        # return names of last jobs to be run
        return ["index_normal_3"], ["index_cancer_3"]

   
    # run local realignment and base recalibration
    def addRealignmentAndRecal(self, previousNormalDep=[], previousCancerDep=[]):
        
        # split by chromosomes
        self.writeJob("split_chrom_normal_1", "48:00:00", "14", "SAMTOOLSPATH view $1 $2 -b > $3",
                      ["$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.dedup.bam","$MULTIPLICITYVAR","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.bam"], 
                      dependencies=previousNormalDep, multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")
        self.writeJob("split_chrom_cancer_1", "48:00:00", "14", "SAMTOOLSPATH view $1 $2 -b > $3",
                      ["$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.bam","$MULTIPLICITYVAR","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.bam"], 
                      dependencies=previousCancerDep, multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")
         
        # sort normal and cancer
        self.writeJob("sort_normal_4", "48:00:00", "14", 
                      "JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate",
                      ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam"], 
                      dependencies=["split_chrom_normal_1"],multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")
        self.writeJob("sort_cancer_4", "48:00:00", "14", 
                      "JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam"], 
                      dependencies=["split_chrom_cancer_1"],multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")
         
        # index normal and cancer
        self.writeJob("index_normal_4", "48:00:00", "14", "JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT",
                      ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam"], 
                      dependencies=["sort_normal_4"], multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")
        self.writeJob("index_cancer_4", "48:00:00", "14", "JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam"], 
                      dependencies=["sort_cancer_4"], multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")
         
        # get realign targets
        self.writeJob("realign_targets", "24:00:00", "8", "JAVAPATH -Xmx4g -jar GATKPATH -T RealignerTargetCreator -R REFERENCEPATH -I $1 -I $2 -known G1KVCF -known MILLSINDELS -o $3 -nt 4",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam","$OUTPUTPATH/$PATIENTID.GATKRealignTargets.$MULTIPLICITYVAR_FORFILE.intervals"], 
                      dependencies=["index_normal_4", "index_cancer_4"], multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")
                 
        # realign bams
        self.writeJob("nwayout", "1:00:00", "2", "python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/OtherScripts/MakeNWayOutFile.py --I $1 --J $2 --K $3 --L $4 --F $5",
                      ["$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam","$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.nwayout.map"], 
                      dependencies=["realign_targets"], multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")
        self.writeJob("realigner", "150:00:00", "8", "JAVAPATH -Xmx4g -jar GATKPATH -T IndelRealigner -R REFERENCEPATH -I $1 -I $2 -known G1KVCF -known MILLSINDELS -targetIntervals $3 -nWayOut $4",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam","$OUTPUTPATH/$PATIENTID.GATKRealignTargets.$MULTIPLICITYVAR_FORFILE.intervals","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.nwayout.map"], 
                      dependencies=["nwayout"], multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")

        # get recalibration data
        self.writeJob("recal_data_normal", "48:00:00", "8", "JAVAPATH -Xmx4g -jar GATKPATH -T BaseRecalibrator -R REFERENCEPATH -I $1 -knownSites DBSNPSITES -knownSites HAPMAPSITES -o $2 -nct 8",
                      ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.grp"], 
                      dependencies=["realigner"],
#                       dependencies=[], 
                      multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")

        self.writeJob("recal_data_cancer", "48:00:00", "8", "JAVAPATH -Xmx4g -jar GATKPATH -T BaseRecalibrator -R REFERENCEPATH -I $1 -knownSites DBSNPSITES -knownSites HAPMAPSITES -o $2 -nct 8",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.grp"], 
                      dependencies=["realigner"],
#                       dependencies=[], 
                      multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")

        # recalibrate
        self.writeJob("recal_normal", "100:00:00", "12", "JAVAPATH -Xmx8g -jar GATKPATH -T PrintReads -R REFERENCEPATH -I $1 -BQSR $2 -o $3 -nct 8 -rf BadCigar",
                      ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.grp","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recal.bam"], 
                      dependencies=["recal_data_normal"], multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")
        self.writeJob("recal_cancer", "100:00:00", "12", "JAVAPATH -Xmx8g -jar GATKPATH -T PrintReads -R REFERENCEPATH -I $1 -BQSR $2 -o $3 -nct 8 -rf BadCigar",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.grp","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recal.bam"], 
                      dependencies=["recal_data_cancer"], multiplicity="1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT")

        # merge bams
        self.writeJob("merge_recal_normal", "48:00:00", "12", "JAVAPATH -Xmx2g -jar PICARDPATH/MergeSamFiles.jar INPUT=$1 INPUT=$2 INPUT=$3 INPUT=$4 INPUT=$5 INPUT=$6 INPUT=$7 INPUT=$8 INPUT=$9 INPUT=${10} INPUT=${11} INPUT=${12} INPUT=${13} OUTPUT=${14} VALIDATION_STRINGENCY=LENIENT USE_THREADING=true",
                      ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.1.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.2.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.3_22.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.4_21.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.5_19.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.6_20.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.7_18.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.8_17.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.9_16.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.10_15.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.11_14.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.12_13.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.X_Y_MT.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.bam"], 
                      dependencies=["recal_normal"])
        self.writeJob("merge_recal_cancer", "48:00:00", "12", "JAVAPATH -Xmx2g -jar PICARDPATH/MergeSamFiles.jar INPUT=$1 INPUT=$2 INPUT=$3 INPUT=$4 INPUT=$5 INPUT=$6 INPUT=$7 INPUT=$8 INPUT=$9 INPUT=${10} INPUT=${11} INPUT=${12} INPUT=${13} OUTPUT=${14} VALIDATION_STRINGENCY=LENIENT USE_THREADING=true",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.1.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.2.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.3_22.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.4_21.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.5_19.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.6_20.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.7_18.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.8_17.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.9_16.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.10_15.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.11_14.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.12_13.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.X_Y_MT.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.bam"], 
                      dependencies=["recal_cancer"])

        # sort normal and cancer
        self.writeJob("sort_normal_5", "48:00:00", "14", 
                      "JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate",
                      ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam"], 
                      dependencies=["merge_recal_normal"])
        self.writeJob("sort_cancer_5", "48:00:00", "14", 
                      "JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam"], 
                      dependencies=["merge_recal_cancer"])
        
        # index normal and cancer
        self.writeJob("index_normal_5", "48:00:00", "14", "JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT",
                      ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam"], 
                      dependencies=["sort_normal_5"])
        self.writeJob("index_cancer_5", "48:00:00", "14", "JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam"], 
                      dependencies=["sort_cancer_5"])

        return ["index_normal_5"], ["index_cancer_5"]
        
    # find germline snvs, also run on cancer to do quality check on shared snvs
    def addGermlineSNVsAndQualityCheck(self, previousNormalDep=[], previousCancerDep=[]): 
        # run haplotype caller
        self.writeJob("haplotype_caller_normal", "150:00:00", "7", "JAVAPATH -Xmx4g -jar GATKPATH -T HaplotypeCaller -R REFERENCEPATH --dbsnp DBSNPSITES -I $1 -o $2 -nct 4 -L $3",
                      ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR.HaplotypeCaller.vcf","$MULTIPLICITYVAR"], 
                      dependencies=previousNormalDep, multiplicity="1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X")
        self.writeJob("haplotype_caller_cancer", "150:00:00", "7", "JAVAPATH -Xmx4g -jar GATKPATH -T HaplotypeCaller -R REFERENCEPATH --dbsnp DBSNPSITES -I $1 -o $2 -nct 4 -L $3",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR.HaplotypeCaller.cancer.vcf","$MULTIPLICITYVAR"], 
                      dependencies=previousCancerDep, multiplicity="1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X")
        # merge unified genotyper results        
        self.writeJob("merge_haplotype_caller_normal", "2:00:00", "8", "module add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/ConcatenateVCFs.R $1 $2 $3",
                      ["$OUTPUTPATH","$PATIENTID.*.HaplotypeCaller.vcf","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.vcf"], 
                      dependencies=["haplotype_caller_normal"])
#                       dependencies=[])
        self.writeJob("merge_haplotype_caller_cancer", "2:00:00", "8", "module add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/ConcatenateVCFs.R $1 $2 $3",
                      ["$OUTPUTPATH","$PATIENTID.*.HaplotypeCaller.cancer.vcf","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.vcf"], 
#                       dependencies=[])
                    dependencies=["haplotype_caller_cancer"])
        
        # variant recalibrate normal variants
        self.writeJob("haplotype_caller_vr_normal_snps", "150:00:00", "12", "JAVAPATH -Xmx4g -jar GATKPATH -T VariantRecalibrator -R REFERENCEPATH -input $1 "+
                      "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 HAPMAPSITES "+
                      "-resource:omni,known=false,training=true,truth=true,prior=12.0 OMNISITES "+ 
                      "-resource:1000G,known=false,training=true,truth=true,prior=10.0 G1KSNPS "+
                      "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 DBSNPSITES "+
                      "-an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum "+
                      "-mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $2 -tranchesFile $3",
                      ["$OUTPUTPATH/$PATIENTID.HaplotypeCaller.vcf","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.snp.recal",
                       "$OUTPUTPATH/$PATIENTID.HaplotypeCaller.snp.tranches"], 
                      dependencies=["merge_haplotype_caller_normal"], multiplicity="")
        self.writeJob("haplotype_caller_ar_normal_snps", "150:00:00", "12", "JAVAPATH -Xmx4g -jar GATKPATH -T ApplyRecalibration -R REFERENCEPATH -input $1 "+
                      "-mode SNP -ts_filter_level 99.9 -recalFile $2 -tranchesFile $3 -o $4",
                      ["$OUTPUTPATH/$PATIENTID.HaplotypeCaller.vcf","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.snp.recal",
                       "$OUTPUTPATH/$PATIENTID.HaplotypeCaller.snp.tranches", "$OUTPUTPATH/$PATIENTID.HaplotypeCaller.snp.recal.indel.raw.vcf"], 
                      dependencies=["haplotype_caller_vr_normal_snps"], multiplicity="")
        self.writeJob("haplotype_caller_vr_normal_indels", "150:00:00", "12", "JAVAPATH -Xmx4g -jar GATKPATH -T VariantRecalibrator -R REFERENCEPATH -input $1 "+
                      "-resource:mills,known=true,training=true,truth=true,prior=12.0 MILLSINDELS "+
                      "-an DP -an FS -an MQRankSum -an ReadPosRankSum "+
                      "-mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $2 -tranchesFile $3",
                      ["$OUTPUTPATH/$PATIENTID.HaplotypeCaller.snp.recal.indel.raw.vcf","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.indel.recal","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.indel.tranches"], 
                      dependencies=["haplotype_caller_ar_normal_snps"], multiplicity="")
        self.writeJob("haplotype_caller_ar_normal_indels", "150:00:00", "12", "JAVAPATH -Xmx4g -jar GATKPATH -T ApplyRecalibration -R REFERENCEPATH -input $1 "+
                      "-mode INDEL -ts_filter_level 99.9 -recalFile $2 -tranchesFile $3 -o $4",
                      ["$OUTPUTPATH/$PATIENTID.HaplotypeCaller.snp.recal.indel.raw.vcf","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.indel.recal",
                       "$OUTPUTPATH/$PATIENTID.HaplotypeCaller.indel.tranches", "$OUTPUTPATH/$PATIENTID.HaplotypeCaller.recal.vcf"], 
                      dependencies=["haplotype_caller_vr_normal_indels"], multiplicity="")

        # variant recalibrate cancer variants
        self.writeJob("haplotype_caller_vr_cancer_snps", "150:00:00", "12", "JAVAPATH -Xmx4g -jar GATKPATH -T VariantRecalibrator -R REFERENCEPATH -input $1 "+
                      "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 HAPMAPSITES "+
                      "-resource:omni,known=false,training=true,truth=true,prior=12.0 OMNISITES "+ 
                      "-resource:1000G,known=false,training=true,truth=true,prior=10.0 G1KSNPS "+
                      "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 DBSNPSITES "+
                      "-an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum "+
                      "-mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $2 -tranchesFile $3",
                      ["$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.vcf","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.snp.recal",
                       "$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.snp.tranches"], 
                      dependencies=["merge_haplotype_caller_cancer"], multiplicity="")
        self.writeJob("haplotype_caller_ar_cancer_snps", "150:00:00", "12", "JAVAPATH -Xmx4g -jar GATKPATH -T ApplyRecalibration -R REFERENCEPATH -input $1 "+
                      "-mode SNP -ts_filter_level 99.9 -recalFile $2 -tranchesFile $3 -o $4",
                      ["$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.vcf","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.snp.recal",
                       "$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.snp.tranches", "$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.snp.recal.indel.raw.vcf"], 
                      dependencies=["haplotype_caller_vr_cancer_snps"], multiplicity="")
        self.writeJob("haplotype_caller_vr_cancer_indels", "150:00:00", "12", "JAVAPATH -Xmx4g -jar GATKPATH -T VariantRecalibrator -R REFERENCEPATH -input $1 "+
                      "-resource:mills,known=true,training=true,truth=true,prior=12.0 MILLSINDELS "+
                      "-an DP -an FS -an MQRankSum -an ReadPosRankSum "+
                      "-mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $2 -tranchesFile $3",
                      ["$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.snp.recal.indel.raw.vcf","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.indel.recal","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.indel.tranches"], 
                      dependencies=["haplotype_caller_ar_cancer_snps"], multiplicity="")
        self.writeJob("haplotype_caller_ar_cancer_indels", "150:00:00", "12", "JAVAPATH -Xmx4g -jar GATKPATH -T ApplyRecalibration -R REFERENCEPATH -input $1 "+
                      "-mode INDEL -ts_filter_level 99.9 -recalFile $2 -tranchesFile $3 -o $4",
                      ["$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.snp.recal.indel.raw.vcf","$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.indel.recal",
                       "$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.indel.tranches", "$OUTPUTPATH/$PATIENTID.HaplotypeCaller.cancer.recal.vcf"], 
                      dependencies=["haplotype_caller_vr_cancer_indels"], multiplicity="")


        # compare haplotype caller results
        self.writeJob("compare_haplotype_caller", "2:00:00", "8", 
                      "source /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/pythonenvs/forDrmaa/bin/activate \\npython /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/QualityControl.py --T $1 --N $2 --O $3",
                      ["$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR.HaplotypeCaller.cancer.recal.vcf","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR.HaplotypeCaller.recal.vcf",
                       "$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR.HaplotypeCallerComp.txt"], 
#                       dependencies=["haplotype_caller_ar_normal_indels", "haplotype_caller_ar_cancer_indels"],
                      dependencies=["haplotype_caller_ar_cancer_indels"], 
                      multiplicity="1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X")
        
        # annotate germline variants
        self.writeJob("germline_stannovar", "12:00:00", "20", "module add stanovar/0.1\\npython /srv/gs1/software/stanovar/stanovar-0.3/annotate_vars.py -i $1 -o $2 -n 1",
                      ["$OUTPUTPATH/$PATIENTID.HaplotypeCaller.recal.vcf","$OUTPUTPATH/$PATIENTID.haplotype.caller"], 
#                       dependencies=[])
                    dependencies=["haplotype_caller_ar_normal_indels"])
        
        # filter relevant variants 
        self.writeJob("filter_germline_stannovar", "6:00:00", "8", "python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/KeepRelevantRegions.py --I $1 --O $2",
                      ["$OUTPUTPATH/$PATIENTID.haplotype.caller.genome_summary.tsv","$OUTPUTPATH/$PATIENTID.haplotype.caller.genome_summary.subset.tsv"], 
                      dependencies=["germline_stannovar"])
        
        # get dgidb info
        self.writeJob("germline_dgidb", "6:00:00", "8", "source /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/pythonenvs/forDrmaa/bin/activate \\nmodule add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/DGIDB_R_for_stannovar.R $1 $2",
                      ["$OUTPUTPATH/$PATIENTID.haplotype.caller.genome_summary.subset.tsv","$OUTPUTPATH/$PATIENTID.haplotype.caller.genome_summary.subset.dgidb.tsv"], 
                      dependencies=["filter_germline_stannovar"])
        
        # add cosmic info
        self.writeJob("germline_cosmic", "6:00:00", "8", "python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/FindCosmicGeneInfo.py --I $1 --O $2 --C $3 --CI $4 --CN $5",
                      ["$OUTPUTPATH/$PATIENTID.haplotype.caller.genome_summary.subset.dgidb.tsv","$OUTPUTPATH/$PATIENTID.haplotype.caller.genome_summary.subset.dgidb.cosmic.tsv","/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv","/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv.index","cosmic_v68_mut_count"], 
                      dependencies=["germline_dgidb"])

        # add cosmic gene info        
        self.writeJob("germline_cosmic_gene", "6:00:00", "8", "python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/FindCosmicGeneMutCountInfo.py --I $1 --O $2 --C $3 --CI $4 --CN $5",
                      ["$OUTPUTPATH/$PATIENTID.haplotype.caller.genome_summary.subset.dgidb.cosmic.tsv","$OUTPUTPATH/$PATIENTID.haplotype.caller.genome_summary.subset.dgidb.cosmic.cosmicgene.tsv","/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv","/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv.index","cosmic_v68_gene_count"], 
                      dependencies=["germline_cosmic"])
        
        # add cancer annotations
        self.writeJob("germline_cancer_gene_annot", "6:00:00", "8", "module add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AddCancerAnnotationToStannovar.R $1 $2 $3",
                      ["$OUTPUTPATH/$PATIENTID.haplotype.caller.genome_summary.subset.dgidb.cosmic.cosmicgene.tsv","$OUTPUTPATH/$PATIENTID.haplotype.caller.genome_summary.subset.dgidb.cosmic.driver.tsv","/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/"], 
                      dependencies=["germline_cosmic_gene"])
        
        return ["germline_cancer_gene_annot"], ["germline_cancer_gene_annot"]
    
    # find somatic snvs and indels and annotate them
    def addMutectAndVarscan2(self, previousNormalDep=[], previousCancerDep=[]):
         
        # run mutect
        self.writeJob("mutect", "150:00:00", "8", "MUTECTJAVAPATH -Xmx4g -jar MUTECTPATH --analysis_type MuTect --reference_sequence REFERENCEPATH --cosmic COSMICVCF --dbsnp DBSNPMUTECTVCF --input_file:normal $1 --input_file:tumor $2 --out $3 --coverage_file $4 -vcf $5",
                      ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","$OUTPUTPATH/$PATIENTID.mutect.callstats.out","$OUTPUTPATH/$PATIENTID.mutect.coverage.wig","$OUTPUTPATH/$PATIENTID.mutect.vcf"],
                      dependencies=(previousNormalDep+previousCancerDep))
  
        # get pileups        
        self.writeJob("pileup_normal", "100:00:00", "8", "SAMTOOLSPATH mpileup -A -q 1 -f REFERENCEPATH $1 > $2",
                      ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.pileup"],
                      dependencies=(previousNormalDep+previousCancerDep))
          
        self.writeJob("pileup_cancer", "100:00:00", "8", "SAMTOOLSPATH mpileup -A -q 1 -f REFERENCEPATH $1 > $2",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.pileup"],
                      dependencies=(previousNormalDep+previousCancerDep))
          
        # run varscan2
        self.writeJob("varscan2", "150:00:00", "8", "JAVAPATH -Xmx4g -jar VARSCANPATH somatic $1 $2 $3",
                      ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.pileup","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.pileup","$OUTPUTPATH/$PATIENTID.varscan"],
                      dependencies=["pileup_cancer", "pileup_normal"])
        
        # filter varscan output
        self.writeJob("varscan2_processSomatic", "150:00:00", "8", "JAVAPATH -Xmx4g -jar VARSCANPATH processSomatic $1",
                      ["$OUTPUTPATH/$PATIENTID.varscan.snp"],
                      dependencies=["varscan2"])
    
        # merge mutations
        self.writeJob("merge_somatic_muts", "6:00:00", "8", "python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/MergeVarscanAndMutect.py --M $1 --V $2 --VI $3 --T $4 --N $5 --O $6",
                      ["$OUTPUTPATH/$PATIENTID.mutect.callstats.out","$OUTPUTPATH/$PATIENTID.varscan.snp","$OUTPUTPATH/$PATIENTID.varscan.indel","cancer","normal","$OUTPUTPATH/$PATIENTID.unannotated.merge"],
                      dependencies=["varscan2", "mutect"])

        # convert merged mutations to vcf
        self.writeJob("merge_somatic_muts2", "6:00:00", "8", "python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/ConvertMergeToVCF.py --I $1 --O $2",
                      ["$OUTPUTPATH/$PATIENTID.unannotated.merge","$OUTPUTPATH/$PATIENTID.unannotated.merge.vcf"],
                      dependencies=["merge_somatic_muts"])
        
        # annotate somatic variants
        self.writeJob("somatic_stannovar", "12:00:00", "20", "module add stanovar/0.1\\npython /srv/gs1/software/stanovar/stanovar-0.3/annotate_vars.py -i $1 -o $2 -n 1",
                      ["$OUTPUTPATH/$PATIENTID.unannotated.merge.vcf","$OUTPUTPATH/$PATIENTID.unannotated.merge"], 
                      dependencies=["merge_somatic_muts2"])
        
        # filter relevant variants 
        self.writeJob("filter_somatic_stannovar", "6:00:00", "8", "python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/KeepRelevantRegions.py --I $1 --O $2",
                      ["$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.tsv","$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.tsv"], 
                      dependencies=["somatic_stannovar"])

        # get dgidb info
        self.writeJob("somatic_dgidb", "6:00:00", "8", "source /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/pythonenvs/forDrmaa/bin/activate \\nmodule add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/DGIDB_R_for_stannovar.R $1 $2",
                      ["$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.tsv","$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.tsv"], 
                      dependencies=["filter_somatic_stannovar"])

        # add cosmic info
        self.writeJob("somatic_cosmic", "6:00:00", "8", "python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/FindCosmicGeneInfo.py --I $1 --O $2 --C $3 --CI $4 --CN $5",
                      ["$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.tsv","$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.cosmic.tsv","/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv","/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv.index","cosmic_v68_mut_count"], 
                      dependencies=["somatic_dgidb"])

        # add cosmic gene info        
        self.writeJob("somatic_cosmic_gene", "6:00:00", "8", "python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/FindCosmicGeneMutCountInfo.py --I $1 --O $2 --C $3 --CI $4 --CN $5",
                      ["$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.cosmic.tsv","$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.cosmic.cosmicgene.tsv","/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv","/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv.index","cosmic_v68_gene_count"], 
                      dependencies=["somatic_cosmic"])
        
        # add cancer annotations
        self.writeJob("somatic_cancer_gene_annot", "6:00:00", "8", "module add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AddCancerAnnotationToStannovar.R $1 $2 $3",
                      ["$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.cosmic.cosmicgene.tsv","$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.cosmic.driver.tsv","/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/"], 
                      dependencies=["somatic_cosmic_gene"])
    
        return ["somatic_cancer_gene_annot"], ["somatic_cancer_gene_annot"]

        
    # find copy number alterations and annotate them
    def addBicseq(self, coverage, previousNormalDep=[], previousCancerDep=[]):
        
        # run bicseq
        self.writeJob("bicseq", "72:00:00", "12", 
                      "module add r/R-2.15.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/OtherScripts/RunBicseQ.R $1 $2 $3 $4 $5",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_X_Y_MT","$OUTPUTPATH/$PATIENTID.BICSEQ.summary.csv", str(coverage)], 
                      dependencies=(previousNormalDep+previousCancerDep))
        
        # analyze bicseq results
        self.writeJob("bicseq_analysis", "6:00:00", "4", 
                      "module add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/CancerDriverBicseqAnalysis.R $1 $2 $3",
                      ["/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/Cancer_Mutation_List.bed","$OUTPUTPATH/$PATIENTID.BICSEQ.summary.csv","$OUTPUTPATH/$PATIENTID.BICSEQ.copy_number.cancer.genes.tsv"], 
                      dependencies=["bicseq"])
        
        return ["bicseq_analysis"], ["bicseq_analysis"]
    
    # find copy number variants
    def addExomeCNV(self, previousNormalDep=[], previousCancerDep=[]):
        # produce coverage files for exome cnv
#         self.writeJob("exome_cnv_coverage_normal", "72:00:00", "12", 
#                       "JAVAPATH -Xmx4g -jar GATKPATH -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R REFERENCEPATH -I $1 -L $2 -o $3",
#                       ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","PERSONALISEXOMEINTERVALLIST","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.coverage"], 
#                       dependencies=(previousNormalDep+previousCancerDep))
#          
#         self.writeJob("exome_cnv_coverage_cancer", "72:00:00", "12", 
#                       "JAVAPATH -Xmx4g -jar GATKPATH -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R REFERENCEPATH -I $1 -L $2 -o $3",
#                       ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","PERSONALISEXOMEINTERVALLIST","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.coverage"], 
#                       dependencies=(previousNormalDep+previousCancerDep))
        # run exome cnv 
        self.writeJob("exome_cnv", "48:00:00", "16", 
                     "module add r/2.15.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/OtherScripts/ExomeCNV.R $1 $2 $3 $4",
                     ["$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted","$OUTPUTPATH/$PATIENTID.exomecnv.output","100"], 
                     dependencies=(previousNormalDep+previousCancerDep))
#                      dependencies=[])
#                      dependencies=["exome_cnv_coverage_cancer", "exome_cnv_coverage_normal"])
        
        # annotate exome cnv result
        self.writeJob("annotate_exome_cnv", "48:00:00", "16", 
             "module add r/3.1.0\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/CancerDriverExomeCNVAnalysis.R $1 $2 $3 $4",
             ["/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/AllGenesBED.collapsed.tsv",
              "$OUTPUTPATH/$PATIENTID.exomecnv.output.cnv.txt",
              "$OUTPUTPATH/$PATIENTID.exomecnv.output.cnv",
              "/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/PancancerCNVDrivers.tsv"], 
#              dependencies=[])
            dependencies=["exome_cnv"])
        
        return ["annotate_exome_cnv"], ["annotate_exome_cnv"]
        
    # find rearrangements and other somatic structural variants
    def addCREST(self, previousNormalDep=[], previousCancerDep=[]):
        
        # copy bam index files with new extensions
        self.writeJob("rename_bam_index", "4:00:00", "4", 
                      "cp $1 $2\\ncp $3 $4",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bai",
                       "$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam.bai",
                       "$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bai",
                       "$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam.bai"], 
                      dependencies=(previousNormalDep+previousCancerDep), multiplicity="")
        
        # extract sclip
        self.writeJob("extract_sclip", "72:00:00", "4", 
                      "module add perl-scg\\nmodule add samtools\\nmodule add crest/1.0\\nextractSClip.pl -i $1 --ref_genome $2 -r $3 -o $4",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","REFERENCEPATH","$MULTIPLICITYVAR","$OUTPUTPATH"], 
                      dependencies=["rename_bam_index"], multiplicity="1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|MT")

        # concat covers
        self.writeJob("concat_cover", "4:00:00", "12", 
                      "echo '$1 > $2'\\ncat $1 > $2",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam.*.cover","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam.cover"], 
                      dependencies=["extract_sclip"])
        
        # run crest
        self.writeJob("crest", "72:00:00", "12", 
                      "module add cap3/1.0\\nmodule add perl-scg\\nmodule add samtools\\nmodule add crest/1.0\\nmodule add blat\\n/srv/gs1/software/blat/3.5/bin/x86_64/gfServer start localhost 8888 $5 &\\nsleep 240\\nCREST.pl -f $1 -d $2 -g $3 --ref_genome $4 -t $5 -r $6 --blatserver localhost --blatport 8888 -o $7 -p $6",
                      ["$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam.cover","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam","REFERENCEPATH","REFERENCETWOBIT","$MULTIPLICITYVAR","$OUTPUTPATH"], 
                      dependencies=["concat_cover"], multiplicity="1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|MT")
        
        # combine crest results
        self.writeJob("concat_crest_results", "4:00:00", "4", 
                      "echo '$1 > $2'\\ncat $1 > $2",
                      ["$OUTPUTPATH/*.predSV.txt", "$OUTPUTPATH/mergedCrestPredSV.txt"], 
                      dependencies=["crest"])
        
    # generate csv file
    def writeCSV(self):
        #self.csvWriter.close()
        self.csvFileHandle.close()
        
    # run pipeline
    def run(self):
        self.writeCSV()
        print os.path.isfile(self.csvFile)
        subprocess.call("python /srv/gsfs0/clinical/cancerPatientAnno/SCGCancerPipeline/ClusterJobManager/Main.py --C "+self.csvFile+" --V /srv/gsfs0/clinical/cancerPatientAnno/SCGCancerPipeline/JobManagementSoftware/Variables.csv"+" --L /srv/gsfs0/clinical/cancerPatientAnno/SCGCancerPipeline/JobManagementSoftware/logfile.txt"+ " --IF ''"+" --OF ''"+" --SF ''"+" --JEF ''"+" --JOF ''"+" --JSF ''", shell=True)
#         subprocess.call("python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/JobManagementSoftware/Main_v2.py --I "+self.csvFile, shell=True)
        
