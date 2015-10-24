import os
import csv
import subprocess

class Patient():
    def __init__(self, patientID, outputDirectory, csvFile):
        self.patientID=patientID
        self.outputDirectory=outputDirectory
        self.csvFile=csvFile
        self.csvFileHandle=open(self.csvFile, 'w')
        self.csvWriter=csv.writer(self.csvFileHandle, delimiter=",")
        self.nextJobNumber=1
        self.csvWriter.writerow(["run","notes","patientID","subTaskGroup","scriptName","outputPath","multiplicity","scriptPath","scriptTime","scriptOutputFileDirectory","scriptErrorFileDirectory","scriptCustomizations","scriptMemory","scriptEmailAddress","scriptCommand","inputs"])
    
    def writeJob(self):
        return ""
    
    def getMultVars(self, filenames):
        multVarNorm=set([])
        multVarCancer=set([])
        for filename in filenames:
            type="normal"
            if "cancer" in filename or "tumor" in filename: type="cancer"
            R="R1"
            if "R1" not in filename: R="R2"
            multvar=filename.strip(".fastq.gz").replace("_"+R, "")
            if type=="cancer": multVarCancer.add(multvar)
            if type=="normal": multVarNorm.add(multvar)
        return "|".join(list(multVarNorm)), "|".join(list(multVarCancer))
    
    def unzipPersonalis(self, zippedDir, filenames):
        multVarNorm=set([])
        multVarCancer=set([])
        for filename in filenames:
            type="normal"
            if "cancer" in filename or "tumor" in filename: type="cancer"
            R="R1"
            if "R1" not in filename: R="R2"
            multvar=filename.strip(".fastq.gz").replace("_"+R, "")
            self.csvWriter.writerow(["TRUE","gunzip",self.patientID,str(self.nextJobNumber),"gunzip",self.outputDirectory,"",
                                 self.outputDirectory+"/scripts/gunzip.sh","6:00:00",self.outputDirectory+"/error_and_outputs",
                                 self.outputDirectory+"/error_and_outputs","#","4","cmelton@stanford.edu",
                                 "gunzip -c $1 > $2",
                                 os.path.join(zippedDir,filename)+"|"+"$OUTPUTPATH/$PATIENTID_"+R+"_"+multvar+"."+type+".fastq"])
            if type=="cancer": multVarCancer.add(multvar)
            if type=="normal": multVarNorm.add(multvar)
        self.nextJobNumber+=1
        return "|".join(list(multVarNorm)), "|".join(list(multVarCancer))
     
    
    def addMultipleFASTQs(self, multVarNorm, multVarCancer, 
                          cancerFASTQ1="$OUTPUTPATH/$PATIENTID_R1_$MULTIPLICITYVAR_FORFILE.cancer.fastq", 
                          cancerFASTQ2="$OUTPUTPATH/$PATIENTID_R2_$MULTIPLICITYVAR_FORFILE.cancer.fastq", 
                          normalFASTQ1="$OUTPUTPATH/$PATIENTID_R1_$MULTIPLICITYVAR_FORFILE.normal.fastq", 
                          normalFASTQ2="$OUTPUTPATH/$PATIENTID_R2_$MULTIPLICITYVAR_FORFILE.normal.fastq", isPersonalis=False):
        self.csvWriter.writerow(["TRUE","bwa",self.patientID,str(self.nextJobNumber),"BWA",self.outputDirectory,multVarNorm,
                                 self.outputDirectory+"/scripts/BWA.sh","150:00:00",self.outputDirectory+"/error_and_outputs",
                                 self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu",
                                 "BWAPATH mem -M REFERENCEPATH $1 $2 -t $3 | SAMTOOLSPATH view -Sbt REFERENCEINDEX -o $4 -",
                                 normalFASTQ1+"|"+normalFASTQ2+"|1|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.bam"])
        self.csvWriter.writerow(["TRUE","bwa",self.patientID,str(self.nextJobNumber),"BWA",self.outputDirectory,multVarCancer,
                                 self.outputDirectory+"/scripts/BWA.sh","150:00:00",self.outputDirectory+"/error_and_outputs",
                                 self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu",
                                 "BWAPATH mem -M REFERENCEPATH $1 $2 -t $3 | SAMTOOLSPATH view -Sbt REFERENCEINDEX -o $4 -",
                                 cancerFASTQ1+"|"+cancerFASTQ2+"|1|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.bam"])
        self.nextJobNumber+=1
        if isPersonalis:
            for m in multVarNorm.split("|"):
                vals=m.split("_")
                rgid, rglb, rgpl, rgpu, rgsm= vals[0], vals[2]+"_"+vals[3], "illumina", vals[2], vals[1]  # id, library (barcode?), platform, platform unit (barcode),   
                self.csvWriter.writerow(["TRUE","add read groups?",self.patientID,str(self.nextJobNumber),"AddReadGroups",
                                         self.outputDirectory,m,self.outputDirectory+"/scripts/AddReadGroups.sh",
                                         "100:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs",
                                         "#","8","cmelton@stanford.edu","JAVAPATH -Xmx2g -jar PICARDPATH/AddOrReplaceReadGroups.jar INPUT=$1 OUTPUT=$2 RGID=$3 RGLB=$4 RGPL=$5 RGPU=$6 RGSM=$7 VALIDATION_STRINGENCY=LENIENT",
                                         "$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.bam|"
                                         +rgid+"|"+rglb+"|"+rgpl+"|"+rgpu+"|"+rgsm])
            for m in multVarCancer.split("|"):
                vals=m.split("_")
                rgid, rglb, rgpl, rgpu, rgsm= vals[0], vals[2]+"_"+vals[3], "illumina", vals[2], vals[1]  # id, library (barcode?), platform, platform unit (barcode),  
                self.csvWriter.writerow(["TRUE","add read groups?",self.patientID,str(self.nextJobNumber),"AddReadGroups",
                                     self.outputDirectory,m,self.outputDirectory+"/scripts/AddReadGroups.sh",
                                     "100:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs",
                                     "#","8","cmelton@stanford.edu","JAVAPATH -Xmx2g -jar PICARDPATH/AddOrReplaceReadGroups.jar INPUT=$1 OUTPUT=$2 RGID=$3 RGLB=$4 RGPL=$5 RGPU=$6 RGSM=$7 VALIDATION_STRINGENCY=LENIENT",
                                     "$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.bam|"
                                     +rgid+"|"+rglb+"|"+rgpl+"|"+rgpu+"|"+rgsm])
        else:
            self.csvWriter.writerow(["TRUE","add read groups?",self.patientID,str(self.nextJobNumber),"AddReadGroups",self.outputDirectory,multVarNorm,self.outputDirectory+"/scripts/AddReadGroups.sh","100:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","JAVAPATH -Xmx2g -jar PICARDPATH/AddOrReplaceReadGroups.jar INPUT=$1 OUTPUT=$2 RGID=$3 RGLB=$4 RGPL=$5 RGPU=$6 RGSM=$7 VALIDATION_STRINGENCY=LENIENT","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.bam|$PATIENTID|$MULTIPLICITYVAR_FORFILE|Illumina|$PATIENTID|normal"])
            self.csvWriter.writerow(["TRUE","add read groups?",self.patientID,str(self.nextJobNumber),"AddReadGroups",self.outputDirectory,multVarCancer,self.outputDirectory+"/scripts/AddReadGroups.sh","100:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","JAVAPATH -Xmx2g -jar PICARDPATH/AddOrReplaceReadGroups.jar INPUT=$1 OUTPUT=$2 RGID=$3 RGLB=$4 RGPL=$5 RGPU=$6 RGSM=$7 VALIDATION_STRINGENCY=LENIENT","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.bam|$PATIENTID|$MULTIPLICITYVAR_FORFILE|Illumina|$PATIENTID|cancer"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","sort normal",self.patientID,str(self.nextJobNumber),"SortSam",self.outputDirectory,multVarNorm,self.outputDirectory+"/scripts/SortSam.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.sorted.bam"])
        self.csvWriter.writerow(["TRUE","sort cancer",self.patientID,str(self.nextJobNumber),"SortSam",self.outputDirectory,multVarCancer,self.outputDirectory+"/scripts/SortSam.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.sorted.bam"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","index normal",self.patientID,str(self.nextJobNumber),"BuildBamIndex",self.outputDirectory,multVarNorm,self.outputDirectory+"/scripts/BuildBamIndex.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.sorted.bam"])
        self.csvWriter.writerow(["TRUE","index cancer",self.patientID,str(self.nextJobNumber),"BuildBamIndex",self.outputDirectory,multVarCancer,self.outputDirectory+"/scripts/BuildBamIndex.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.sorted.bam"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","mark duplicates?",self.patientID,str(self.nextJobNumber),"MarkDuplicates",self.outputDirectory,multVarNorm,self.outputDirectory+"/scripts/MarkDuplicates.sh","100:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","32","cmelton@stanford.edu","JAVAPATH -Xmx24g -jar PICARDPATH/MarkDuplicates.jar INPUT=$1 OUTPUT=$2 METRICS_FILE=$3 VALIDATION_STRINGENCY=LENIENT","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.sorted.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.sorted.dedup.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.sorted.dedup.metrics.txt"])
        self.csvWriter.writerow(["TRUE","mark duplicates?",self.patientID,str(self.nextJobNumber),"MarkDuplicates",self.outputDirectory,multVarCancer,self.outputDirectory+"/scripts/MarkDuplicates.sh","100:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","32","cmelton@stanford.edu","JAVAPATH -Xmx24g -jar PICARDPATH/MarkDuplicates.jar INPUT=$1 OUTPUT=$2 METRICS_FILE=$3 VALIDATION_STRINGENCY=LENIENT","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.sorted.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.sorted.dedup.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.sorted.dedup.metrics.txt"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","index normal",self.patientID,str(self.nextJobNumber),"BuildBamIndex",self.outputDirectory,multVarNorm,self.outputDirectory+"/scripts/BuildBamIndex.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.normal.RG.sorted.dedup.bam"])
        self.csvWriter.writerow(["TRUE","index cancer",self.patientID,str(self.nextJobNumber),"BuildBamIndex",self.outputDirectory,multVarCancer,self.outputDirectory+"/scripts/BuildBamIndex.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.cancer.RG.sorted.dedup.bam"])
        self.nextJobNumber+=1
        # merge
        self.csvWriter.writerow(["TRUE","merge normal bams",self.patientID,str(self.nextJobNumber),"MergeBams",self.outputDirectory,"",self.outputDirectory+"/scripts/MergeBams.sh",
                                 "48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu",
                                 "JAVAPATH -Xmx2g -jar PICARDPATH/MergeSamFiles.jar "+" ".join(map(lambda x: "INPUT=${"+str(x+1)+"}", range(len(multVarNorm.split("|")))))+" OUTPUT=${"+str(len(multVarNorm.split("|"))+1)+"} VALIDATION_STRINGENCY=LENIENT USE_THREADING=true  ",
                                 "|".join(map(lambda x: "$OUTPUTPATH/$PATIENTID."+x+".normal.RG.sorted.dedup.bam", multVarNorm.split("|")))+"|$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.dedup.bam"])
        self.csvWriter.writerow(["TRUE","merge normal bams",self.patientID,str(self.nextJobNumber),"MergeBams",self.outputDirectory,"",self.outputDirectory+"/scripts/MergeBams.sh",
                                 "48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu",
                                 "JAVAPATH -Xmx2g -jar PICARDPATH/MergeSamFiles.jar "+" ".join(map(lambda x: "INPUT=${"+str(x+1)+"}", range(len(multVarCancer.split("|")))))+" OUTPUT=${"+str(len(multVarCancer.split("|"))+1)+"} VALIDATION_STRINGENCY=LENIENT USE_THREADING=true  ",
                                 "|".join(map(lambda x: "$OUTPUTPATH/$PATIENTID."+x+".cancer.RG.sorted.dedup.bam", multVarCancer.split("|")))+"|$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.bam"])
        self.nextJobNumber+=1
        # sort and index
        self.csvWriter.writerow(["TRUE","sort normal",self.patientID,str(self.nextJobNumber),"SortSam",self.outputDirectory,"",self.outputDirectory+"/scripts/SortSam.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate","$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.dedup.bam|$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.dedup.sorted.bam"])
        self.csvWriter.writerow(["TRUE","sort tumor",self.patientID,str(self.nextJobNumber),"SortSam",self.outputDirectory,"",self.outputDirectory+"/scripts/SortSam.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate","$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.bam|$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.sorted.bam"]) 
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","index bam?",self.patientID,str(self.nextJobNumber),"BuildBamIndex",self.outputDirectory,"",self.outputDirectory+"/scripts/BuildBamIndex.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT\\ncp $2 $3","$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.dedup.sorted.bam|$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.dedup.sorted.bai|$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.dedup.sorted.bam.bai"])
        self.csvWriter.writerow(["TRUE","index bam?",self.patientID,str(self.nextJobNumber),"BuildBamIndex",self.outputDirectory,"",self.outputDirectory+"/scripts/BuildBamIndex.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT\\ncp $2 $3","$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.sorted.bam|$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.sorted.bai|$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.sorted.bam.bai"])
        self.nextJobNumber+=1
#         self.csvWriter.writerow(["TRUE","remove additional bams and fastqs",self.patientID,str(self.nextJobNumber),"RemoveBams",
#                                  self.outputDirectory,multVarNorm+"|"+multVarCancer,self.outputDirectory+"/scripts/RemoveBams.sh","6:00:00",
#                                  self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs",
#                                  "#","2","cmelton@stanford.edu","rm $1","$OUTPUTPATH/*$MULTIPLICITYVAR_FORFILE*.bam"])
        self.nextJobNumber+=1
        
   
    # run local realignment and base recalibration
    def addRealignmentAndRecal(self):
#         self.csvWriter.writerow(["TRUE","break normal up by chrom",self.patientID,str(self.nextJobNumber),"SegregateByChrom",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/SegregateByChrom.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","4","cmelton@stanford.edu","SAMTOOLSPATH view $1 $2 -b > $3","$OUTPUTPATH/$PATIENTID.merged.cancer.RG.sorted.dedup.sorted.bam|$MULTIPLICITYVAR|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.bam"])
#         self.csvWriter.writerow(["TRUE","break normal up by chrom",self.patientID,str(self.nextJobNumber),"SegregateByChrom",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/SegregateByChrom.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","4","cmelton@stanford.edu","SAMTOOLSPATH view $1 $2 -b > $3","$OUTPUTPATH/$PATIENTID.merged.normal.RG.sorted.dedup.sorted.bam|$MULTIPLICITYVAR|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.bam"])
#         self.nextJobNumber+=1
#         self.csvWriter.writerow(["TRUE","sort cancer",self.patientID,str(self.nextJobNumber),"SortSam",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/SortSam.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam"])
#         self.csvWriter.writerow(["TRUE","sort normal",self.patientID,str(self.nextJobNumber),"SortSam",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/SortSam.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam"])
#         self.nextJobNumber+=1
#         self.csvWriter.writerow(["TRUE","index bam?",self.patientID,str(self.nextJobNumber),"BuildBamIndex",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/BuildBamIndex.sh","100:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam"])
#         self.csvWriter.writerow(["TRUE","index bam?",self.patientID,str(self.nextJobNumber),"BuildBamIndex",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/BuildBamIndex.sh","100:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam"])
#         self.nextJobNumber+=1
#         self.csvWriter.writerow(["TRUE","get realign targets for each chrom",self.patientID,str(self.nextJobNumber),"RealignTargetScript",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/RealignTargetScript.sh","24:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar GATKPATH -T RealignerTargetCreator -R REFERENCEPATH -I $1 -I $2 -known G1KVCF -known MILLSINDELS -o $3 -nt 4","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam|$OUTPUTPATH/$PATIENTID.GATKRealignTargets.$MULTIPLICITYVAR_FORFILE.intervals"])
#         self.nextJobNumber+=1
#         self.csvWriter.writerow(["TRUE","make nway out file for realignerscript",self.patientID,str(self.nextJobNumber),"MakeNWayOutFile",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/MakeNWayOutFile.sh","1:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","2","cmelton@stanford.edu","python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/OtherScripts/MakeNWayOutFile.py --I $1 --J $2 --K $3 --L $4 --F $5","$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam|$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.nwayout.map"])
#         self.nextJobNumber+=1
#         self.csvWriter.writerow(["TRUE","realign each tumor/normal chrom pair",self.patientID,str(self.nextJobNumber),"RealignerScript",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/RealignerScript.sh","150:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar GATKPATH -T IndelRealigner -R REFERENCEPATH -I $1 -I $2 -known G1KVCF -known MILLSINDELS -targetIntervals $3 -nWayOut $4","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.bam|$OUTPUTPATH/$PATIENTID.GATKRealignTargets.$MULTIPLICITYVAR_FORFILE.intervals|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR_FORFILE.nwayout.map"])
#         self.nextJobNumber+=1
#         self.csvWriter.writerow(["TRUE","get recalibration data for each tumor chrom",self.patientID,str(self.nextJobNumber),"GetRecalibrationDataScript",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/GetRecalibrationDataScript.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar GATKPATH -T BaseRecalibrator -R REFERENCEPATH -I $1 -knownSites DBSNPSITES -knownSites HAPMAPSITES -o $2 -nct 8 -plots $3 --intermediate_csv_file $4","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.grp|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.pdf|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.csv"])
#         self.csvWriter.writerow(["TRUE","get recalibration data for each normal chrom",self.patientID,str(self.nextJobNumber),"GetRecalibrationDataScript",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/GetRecalibrationDataScript.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar GATKPATH -T BaseRecalibrator -R REFERENCEPATH -I $1 -knownSites DBSNPSITES -knownSites HAPMAPSITES -o $2 -nct 8 -plots $3 --intermediate_csv_file $4","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.grp|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.pdf|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.csv"])
#         self.nextJobNumber+=1
#         self.csvWriter.writerow(["TRUE","recalibrate each tumor chrom ",self.patientID,str(self.nextJobNumber),"RecalibrationBamScript",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/RecalibrationBamScript.sh","100:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu","JAVAPATH -Xmx8g -jar GATKPATH -T PrintReads -R REFERENCEPATH -I $1 -BQSR $2 -o $3 -nct 8 -rf BadCigar","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.grp|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recal.bam"])
#         self.csvWriter.writerow(["TRUE","recalibrate each normal chrom ",self.patientID,str(self.nextJobNumber),"RecalibrationBamScript",self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/RecalibrationBamScript.sh","100:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu","JAVAPATH -Xmx8g -jar GATKPATH -T PrintReads -R REFERENCEPATH -I $1 -BQSR $2 -o $3 -nct 8 -rf BadCigar","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recalData.grp|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.$MULTIPLICITYVAR_FORFILE.sorted.realigned.recal.bam"])
#         self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","merge bams into one prior to running variant calling algorithms",self.patientID,str(self.nextJobNumber),"MergeBams",self.outputDirectory,"",self.outputDirectory+"/scripts/MergeBams.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu","JAVAPATH -Xmx2g -jar PICARDPATH/MergeSamFiles.jar INPUT=$1 INPUT=$2 INPUT=$3 INPUT=$4 INPUT=$5 INPUT=$6 INPUT=$7 INPUT=$8 INPUT=$9 INPUT=${10} INPUT=${11} INPUT=${12} INPUT=${13} OUTPUT=${14} VALIDATION_STRINGENCY=LENIENT USE_THREADING=true","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.1.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.2.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.3_22.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.4_21.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.5_19.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.6_20.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.7_18.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.8_17.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.9_16.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.10_15.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.11_14.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.12_13.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.X_Y_MT.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.bam"])
        self.csvWriter.writerow(["TRUE","merge bams into one prior to running variant calling algorithms",self.patientID,str(self.nextJobNumber),"MergeBams",self.outputDirectory,"",self.outputDirectory+"/scripts/MergeBams.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu","JAVAPATH -Xmx2g -jar PICARDPATH/MergeSamFiles.jar INPUT=$1 INPUT=$2 INPUT=$3 INPUT=$4 INPUT=$5 INPUT=$6 INPUT=$7 INPUT=$8 INPUT=$9 INPUT=${10} INPUT=${11} INPUT=${12} INPUT=${13} OUTPUT=${14} VALIDATION_STRINGENCY=LENIENT USE_THREADING=true","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.1.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.2.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.3_22.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.4_21.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.5_19.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.6_20.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.7_18.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.8_17.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.9_16.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.10_15.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.11_14.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.12_13.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.X_Y_MT.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.bam"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","sort tumor",self.patientID,str(self.nextJobNumber),"SortSam",self.outputDirectory,"",self.outputDirectory+"/scripts/SortSam.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam"])
        self.csvWriter.writerow(["TRUE","sort normal",self.patientID,str(self.nextJobNumber),"SortSam",self.outputDirectory,"",self.outputDirectory+"/scripts/SortSam.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xms8g -Xmx8g -jar PICARDPATH/SortSam.jar INPUT=$1 OUTPUT=$2 MAX_RECORDS_IN_RAM=$((8*250000)) VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","index bam?",self.patientID,str(self.nextJobNumber),"BuildBamIndex",self.outputDirectory,"",self.outputDirectory+"/scripts/BuildBamIndex.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT\\ncp $2 $3","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bai|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam.bai"])
        self.csvWriter.writerow(["TRUE","index bam?",self.patientID,str(self.nextJobNumber),"BuildBamIndex",self.outputDirectory,"",self.outputDirectory+"/scripts/BuildBamIndex.sh","48:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","14","cmelton@stanford.edu","JAVAPATH -Xmx4g -jar PICARDPATH/BuildBamIndex.jar INPUT=$1 VALIDATION_STRINGENCY=LENIENT\\ncp $2 $3","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bai|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam.bai"])
        self.nextJobNumber+=1
#         self.csvWriter.writerow(["TRUE","remove additional bams and fastqs",self.patientID,str(self.nextJobNumber),"RemoveBams",
#                          self.outputDirectory,"1|2|3 22|4 21|5 19|6 20|7 18|8 17|9 16|10 15|11 14|12 13|X Y MT",self.outputDirectory+"/scripts/RemoveBams.sh","6:00:00",
#                          self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs",
#                          "#","2","cmelton@stanford.edu","rm $1","$OUTPUTPATH/*dedup.$MULTIPLICITYVAR_FORFILE.*"])
        self.nextJobNumber+=1
        
    # find germline snvs, also run on cancer to do quality check on shared snvs
    def addGermlineSNVsAndQualityCheck(self): 
        resetVal=self.nextJobNumber
        self.csvWriter.writerow(["TRUE","run unified genotyper",self.patientID,str(self.nextJobNumber),"UnifiedGenotyper",self.outputDirectory,
                                 "1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X",self.outputDirectory+"/scripts/UnifiedGenotyper.sh",
                                 "150:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu",
                                 "JAVAPATH -Xmx4g -jar GATKPATH -T UnifiedGenotyper -R REFERENCEPATH --dbsnp DBSNPSITES -I $1 -o $2 -nct 4 -L $3",
                                 "$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR.UnifiedGenotyper.vcf|$MULTIPLICITYVAR"])
        self.csvWriter.writerow(["TRUE","run unified genotyper",self.patientID,str(self.nextJobNumber),"UnifiedGenotyper",self.outputDirectory,
                                 "1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X",self.outputDirectory+"/scripts/UnifiedGenotyper.sh",
                                 "150:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu",
                                 "JAVAPATH -Xmx4g -jar GATKPATH -T UnifiedGenotyper -R REFERENCEPATH --dbsnp DBSNPSITES -I $1 -o $2 -nct 4 -L $3",
                                 "$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR.UnifiedGenotyper.cancer.vcf|$MULTIPLICITYVAR"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","merge unified gennotyper results",self.patientID,str(self.nextJobNumber),"MergeVCFs",self.outputDirectory,"",self.outputDirectory+"/scripts/MergeVCFs.sh",
                                 "2:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu",
                                 "module add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/ConcatenateVCFs.R $1 $2 $3",
                                 "$OUTPUTPATH|$PATIENTID.*.UnifiedGenotyper.vcf|$OUTPUTPATH/$PATIENTID.UnifiedGenotyper.vcf"])
        self.csvWriter.writerow(["TRUE","merge unified gennotyper results",self.patientID,str(self.nextJobNumber),"MergeVCFs",self.outputDirectory,"",self.outputDirectory+"/scripts/MergeVCFs.sh",
                                 "2:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu",
                                 "module add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/ConcatenateVCFs.R $1 $2 $3",
                                 "$OUTPUTPATH|$PATIENTID.*.UnifiedGenotyper.cancer.vcf|$OUTPUTPATH/$PATIENTID.UnifiedGenotyper.cancer.vcf"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","CompareVCFs",self.patientID,str(self.nextJobNumber),"CompareVCFs",self.outputDirectory,
                                 "1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X",self.outputDirectory+"/scripts/CompareVCFs.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","source /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/pythonenvs/forDrmaa/bin/activate \\npython /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/QualityControl.py --T $1 --N $2 --O $3","$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR.UnifiedGenotyper.cancer.vcf|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR.UnifiedGenotyper.vcf|$OUTPUTPATH/$PATIENTID.$MULTIPLICITYVAR.UnifiedGenotyperComp.txt"])
        self.csvWriter.writerow(["TRUE","RunStannovar",self.patientID,str(self.nextJobNumber),"RunStannovar",self.outputDirectory,"",self.outputDirectory+"/scripts/RunStannovar.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","20","cmelton@stanford.edu","module add stanovar/0.1\\npython /srv/gs1/software/stanovar/stanovar-0.3/annotate_vars.py -i $1 -o $2 -n 1","$OUTPUTPATH/$PATIENTID.UnifiedGenotyper.vcf|$OUTPUTPATH/$PATIENTID.unified.genotyper"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","KeepRelevantRegions",self.patientID,str(self.nextJobNumber),"KeepRelevantRegions",self.outputDirectory,"",self.outputDirectory+"/scripts/KeepRelevantRegions.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/KeepRelevantRegions.py --I $1 --O $2","$OUTPUTPATH/$PATIENTID.unified.genotyper.genome_summary.tsv|$OUTPUTPATH/$PATIENTID.unified.genotyper.genome_summary.subset.tsv"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","GetDGIDB",self.patientID,str(self.nextJobNumber),"GetDGIDB",self.outputDirectory,"",self.outputDirectory+"/scripts/GetDGIDB.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","source /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/pythonenvs/forDrmaa/bin/activate \\nmodule add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/DGIDB_R_for_stannovar.R $1 $2","$OUTPUTPATH/$PATIENTID.unified.genotyper.genome_summary.subset.tsv|$OUTPUTPATH/$PATIENTID.unified.genotyper.genome_summary.subset.dgidb.tsv"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","AddCosmic",self.patientID,str(self.nextJobNumber),"AddCosmic",self.outputDirectory,"",self.outputDirectory+"/scripts/AddCosmic.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/FindCosmicGeneInfo.py --I $1 --O $2 --C $3 --CI $4 --CN $5","$OUTPUTPATH/$PATIENTID.unified.genotyper.genome_summary.subset.dgidb.tsv|$OUTPUTPATH/$PATIENTID.unified.genotyper.genome_summary.subset.dgidb.cosmic.tsv|/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv|/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv.index|cosmic_v68_mut_count"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","AddCosmicGene",self.patientID,str(self.nextJobNumber),"AddCosmicGene",self.outputDirectory,"",self.outputDirectory+"/scripts/AddCosmicGene.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/FindCosmicGeneMutCountInfo.py --I $1 --O $2 --C $3 --CI $4 --CN $5","$OUTPUTPATH/$PATIENTID.unified.genotyper.genome_summary.subset.dgidb.cosmic.tsv|$OUTPUTPATH/$PATIENTID.unified.genotyper.genome_summary.subset.dgidb.cosmic.cosmicgene.tsv|/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv|/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv.index|cosmic_v68_gene_count"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","AddCancerAnnotationToStannovar",self.patientID,str(self.nextJobNumber),"AddCancerAnnotationToStannovar",self.outputDirectory,"",self.outputDirectory+"/scripts/AddCancerAnnotationToStannovar.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","module add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AddCancerAnnotationToStannovar.R $1 $2 $3","$OUTPUTPATH/$PATIENTID.unified.genotyper.genome_summary.subset.dgidb.cosmic.cosmicgene.tsv|$OUTPUTPATH/$PATIENTID.unified.genotyper.genome_summary.subset.dgidb.cosmic.driver.tsv|/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/"])
        self.nextJobNumber+=1
        self.nextJobNumber=resetVal
    
    # find somatic snvs and indels and annotate them
    def addMutectAndVarscan2(self):
        resetVal=self.nextJobNumber
        self.csvWriter.writerow(["TRUE","run mutect on each tumor normal pair",self.patientID,str(self.nextJobNumber),"writeMutectScript",self.outputDirectory,"",self.outputDirectory+"/scripts/writeMutectScript.sh","150:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","MUTECTJAVAPATH -Xmx4g -jar MUTECTPATH --analysis_type MuTect --reference_sequence REFERENCEPATH --cosmic COSMICVCF --dbsnp DBSNPMUTECTVCF --input_file:normal $1 --input_file:tumor $2 --out $3 --coverage_file $4 -vcf $5","$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|$OUTPUTPATH/$PATIENTID.mutect.callstats.out|$OUTPUTPATH/$PATIENTID.mutect.coverage.wig|$OUTPUTPATH/$PATIENTID.mutect.vcf"])
        self.csvWriter.writerow(["TRUE","get pileup file for varscan2",self.patientID,str(self.nextJobNumber),
                                 "writePileupScript",self.outputDirectory,"",self.outputDirectory+"/scripts/writePileupScript.sh","72:00:00",
                                 self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu",
                                 "SAMTOOLSPATH mpileup -A -q 1 -f REFERENCEPATH $1 > $2",
                                 "$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.pileup"])
             
        self.csvWriter.writerow(["TRUE","get pileup file for varscan2",self.patientID,str(self.nextJobNumber),
                                 "writePileupScript",self.outputDirectory,"",self.outputDirectory+"/scripts/writePileupScript.sh","72:00:00",
                                 self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu",
                                 "SAMTOOLSPATH mpileup -A -q 1 -f REFERENCEPATH $1 > $2",
                                 "$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|"+
                                 "$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.pileup"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","run varscan2 on each tumor/normal pair",self.patientID,str(self.nextJobNumber),
                                 "VarscanSomaticScript",self.outputDirectory,"",self.outputDirectory+"/scripts/VarscanSomaticScript.sh","150:00:00",
                                 self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu",
                                 "JAVAPATH -Xmx4g -jar VARSCANPATH somatic $1 $2 $3",
                                 "$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.pileup|"+
                                 "$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.pileup|"+
                                 "$OUTPUTPATH/$PATIENTID.varscan"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","MergeMutationsScript.sh",self.patientID,str(self.nextJobNumber),"MergeMutationsScript",self.outputDirectory,"",self.outputDirectory+"/scripts/MergeMutationsScript.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/MergeVarscanAndMutect.py --M $1 --V $2 --VI $3 --T $4 --N $5 --O $6","$OUTPUTPATH/$PATIENTID.mutect.callstats.out|$OUTPUTPATH/$PATIENTID.varscan.snp|$OUTPUTPATH/$PATIENTID.varscan.indel|cancer|normal|$OUTPUTPATH/$PATIENTID.unannotated.merge"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","ConvertMergeToVCF",self.patientID,str(self.nextJobNumber),"MergeMutationsScript",self.outputDirectory,"",self.outputDirectory+"/scripts/ConvertMergeToVCF.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/ConvertMergeToVCF.py --I $1 --O $2","$OUTPUTPATH/$PATIENTID.unannotated.merge|$OUTPUTPATH/$PATIENTID.unannotated.merge.vcf"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","RunStannovar",self.patientID,str(self.nextJobNumber),"RunStannovar",self.outputDirectory,"",self.outputDirectory+"/scripts/RunStannovar.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","module add stanovar/0.1\\npython /srv/gs1/software/stanovar/stanovar-0.3/annotate_vars.py -i $1 -o $2 -n 1","$OUTPUTPATH/$PATIENTID.unannotated.merge.vcf|$OUTPUTPATH/$PATIENTID.unannotated.merge"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","KeepRelevantRegions",self.patientID,str(self.nextJobNumber),"KeepRelevantRegions",self.outputDirectory,"",self.outputDirectory+"/scripts/KeepRelevantRegions.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/KeepRelevantRegions.py --I $1 --O $2","$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.tsv|$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.tsv"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","GetDGIDB",self.patientID,str(self.nextJobNumber),"GetDGIDB",self.outputDirectory,"",self.outputDirectory+"/scripts/GetDGIDB.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","source /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/pythonenvs/forDrmaa/bin/activate \\nmodule add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/DGIDB_R_for_stannovar.R $1 $2","$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.tsv|$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.tsv"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","AddCosmic",self.patientID,str(self.nextJobNumber),"AddCosmic",self.outputDirectory,"",self.outputDirectory+"/scripts/AddCosmic.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/FindCosmicGeneInfo.py --I $1 --O $2 --C $3 --CI $4 --CN $5","$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.tsv|$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.cosmic.tsv|/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv|/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv.index|cosmic_v68_mut_count"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","AddCosmicGene",self.patientID,str(self.nextJobNumber),"AddCosmicGene",self.outputDirectory,"",self.outputDirectory+"/scripts/AddCosmicGene.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/FindCosmicGeneMutCountInfo.py --I $1 --O $2 --C $3 --CI $4 --CN $5","$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.cosmic.tsv|$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.cosmic.cosmicgene.tsv|/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv|/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/CosmicMutantExportIncFus_v68.tsv.index|cosmic_v68_gene_count"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","AddCancerAnnotationToStannovar",self.patientID,str(self.nextJobNumber),"AddCancerAnnotationToStannovar",self.outputDirectory,"",self.outputDirectory+"/scripts/AddCancerAnnotationToStannovar.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","8","cmelton@stanford.edu","module add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AddCancerAnnotationToStannovar.R $1 $2 $3","$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.cosmic.cosmicgene.tsv|$OUTPUTPATH/$PATIENTID.unannotated.merge.genome_summary.subset.dgidb.cosmic.driver.tsv| /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/"])
        self.nextJobNumber+=1
        self.nextJobNumber=resetVal
        
    # find copy number alterations and annotate them
    def addBicseq(self, coverage):
        resetVal=self.nextJobNumber
        self.csvWriter.writerow(["TRUE","BicseQ",self.patientID,str(self.nextJobNumber),"BicseQ",self.outputDirectory,"",self.outputDirectory+"/scripts/BicseQ.sh","72:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu","module add r/R-2.15.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/OtherScripts/RunBicseQ.R $1 $2 $3 $4 $5","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_X_Y_MT|$OUTPUTPATH/$PATIENTID.BICSEQ.summary.csv|"+coverage])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","BicseQAnalysis",self.patientID,str(self.nextJobNumber),"BicseQAnalysis",self.outputDirectory,"",self.outputDirectory+"/scripts/BicseQAnalysis.sh","6:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","4","cmelton@stanford.edu","module add r/3.0.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/CancerDriverBicseqAnalysis.R $1 $2 $3","/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/Annotation/AnnotationFiles/Cancer_Mutation_List.bed|$OUTPUTPATH/$PATIENTID.BICSEQ.summary.csv|$OUTPUTPATH/$PATIENTID.BICSEQ.copy_number.cancer.genes.tsv"])
        self.nextJobNumber+=1
        self.nextJobNumber=resetVal
    
    # find copy number variants
    def addExomeCNV(self):
        resetVal=self.nextJobNumber
        self.csvWriter.writerow(["TRUE","produceCoverage",self.patientID,str(self.nextJobNumber),"produceCoverage",self.outputDirectory,"",
                                 self.outputDirectory+"/scripts/produceCoverage.sh","72:00:00",self.outputDirectory+"/error_and_outputs",
                                 self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu",
                                 "JAVAPATH -Xmx4g -jar GATKPATH -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R REFERENCEPATH -I $1 -L $2 -o $3",
                                 "$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|PERSONALISEXOMEINTERVALLIST|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.coverage"])
        self.csvWriter.writerow(["TRUE","produceCoverage",self.patientID,str(self.nextJobNumber),"produceCoverage",self.outputDirectory,"",
                                 self.outputDirectory+"/scripts/produceCoverage.sh","72:00:00",self.outputDirectory+"/error_and_outputs",
                                 self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu",
                                 "JAVAPATH -Xmx4g -jar GATKPATH -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R REFERENCEPATH -I $1 -L $2 -o $3",
                                 "$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|PERSONALISEXOMEINTERVALLIST|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.coverage"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","ExomeCNV",self.patientID,str(self.nextJobNumber),"ExomeCNV",self.outputDirectory,"",self.outputDirectory+"/scripts/ExomeCNV.sh","48:00:00",
                                 self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","16","cmelton@stanford.edu",
                                 "module add r/R-2.15.1\\nRscript /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/OtherScripts/ExomeCNV.R $1 $2 $3 $4",
                                 "$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted|$OUTPUTPATH/$PATIENTID.exomecnv.output|100"])
        self.nextJobNumber=resetVal
        
    # find rearrangements and other somatic structural variants
    def addCREST(self):
        resetVal=self.nextJobNumber
        self.csvWriter.writerow(["TRUE","extractSClip",self.patientID,str(self.nextJobNumber),"extractSClip",self.outputDirectory,"1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|MT",self.outputDirectory+"/scripts/extractSClip.sh","72:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","4","cmelton@stanford.edu","module add perl-scg\\nmodule add samtools\\nmodule add crest/1.0\\nextractSClip.pl -i $1 --ref_genome $2 -r $3 -o $4","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|REFERENCEPATH|$MULTIPLICITYVAR|$OUTPUTPATH"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","concatCovers",self.patientID,str(self.nextJobNumber),"catCover",self.outputDirectory,"",self.outputDirectory+"/scripts/catCover.sh","4:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu","echo '$1 > $2'\\ncat $1 > $2","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam.*.cover|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam.cover"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","CREST",self.patientID,str(self.nextJobNumber),"CREST",self.outputDirectory,"1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|MT",self.outputDirectory+"/scripts/CREST.sh","72:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","12","cmelton@stanford.edu","module add cap3/1.0\\nmodule add perl-scg\\nmodule add samtools\\nmodule add crest/1.0\\nmodule add blat\\nscreen -dm /srv/gs1/software/blat/3.5/bin/x86_64/gfServer start localhost 8888 $5\\n/srv/gs1/projects/snyder/collinmelton/bundle/2.3/b37/d5/hs37d5.2bit\\nsleep 240\\nCREST.pl -f $1 -d $2 -g $3 --ref_genome $4 -t $5 -r $6 --blatserver localhost --blatport 8888 -o $7 -p $6","$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam.cover|$OUTPUTPATH/$PATIENTID.cancer.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|$OUTPUTPATH/$PATIENTID.normal.RG.sorted.dedup.merged.sorted.realigned.recal.sorted.bam|REFERENCEPATH|REFERENCETWOBIT|$MULTIPLICITYVAR|$OUTPUTPATH"])
        self.nextJobNumber+=1
        self.csvWriter.writerow(["TRUE","concatCrestResults",self.patientID,str(self.nextJobNumber),"concatCrestResults",self.outputDirectory,"",self.outputDirectory+"/scripts/concatCrestResults.sh","4:00:00",self.outputDirectory+"/error_and_outputs",self.outputDirectory+"/error_and_outputs","#","4","cmelton@stanford.edu","echo '$1 > $2'\\ncat $1 > $2","$OUTPUTPATH/*.predSV.txt|$OUTPUTPATH/mergedCrestPredSV.txt"])
        self.nextJobNumber+=1
        self.nextJobNumber=resetVal
        
    # generate csv file
    def writeCSV(self):
        #self.csvWriter.close()
        self.csvFileHandle.close()
        
    # run pipeline
    def run(self):
        self.writeCSV()
        print os.path.isfile(self.csvFile)
        subprocess.call("python /srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/JobManagementSoftware/Main.py --I "+self.csvFile, shell=True)