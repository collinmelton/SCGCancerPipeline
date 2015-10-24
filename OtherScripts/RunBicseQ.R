
library("BICseq")

args<-commandArgs(TRUE)
cancerBam<-args[1] #"/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/PairedAnalysis/601453_NA5K5CHE.cancer.RG.sorted.dedup.X_Y_MT.sorted.realigned.recal.bam"
normalBam<-args[2] #"/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/PairedAnalysis/601453_NA5K5CHE.normal.RG.sorted.dedup.X_Y_MT.sorted.realigned.recal.bam"
chroms<-strsplit(args[3], split="_", fixed=T)[[1]]#c(1:22, "X", "Y")
outputFile<-args[4]
coverage<-args[5]
print(chroms)
print(outputFile)

bicseq<-BICseq(sample=cancerBam, reference=normalBam, seqNames=chroms)

# lambda=4 good for high coverage (5X-30X), 2 for medium (2-5X), and 1 or 1.2 for low coverage
if (coverage=="medium") {
  lambda=2
} else if (coverage=="low") {
  lambda=1.2
} else if (coverage=="high") {
  lambda=4
} else {
  lambda=2
}
segs <- getBICseg(object = bicseq, bin = 100, lambda = lambda, winSize = 200,quant = 0.95, mult = 1)

# optional plot of copy number
# plot(segs, sampleName =outputImage, save = T, plotBin=TRUE, chromOrder=chroms)

# get segmentation summary and write to output
seg.summary <- BICseq:::getSummary(segs,correction=TRUE)
write.csv(seg.summary, file=outputFile, quote=F, row.names=F)