library(ExomeCNV)

# read in command line args
args<-commandArgs(TRUE)
normalPrefix=args[1]
cancerPrefix=args[2]
outputFile=args[3]
readlength=as.numeric(args[4])

# define list of chromosomes
#chr.list=paste("chr",c(as.character(1:22), "X"),sep="")

# read in coverage files
suffix = ".coverage.sample_interval_summary"
normal = read.coverage.gatk(paste(normalPrefix, suffix, sep=""))
tumor = read.coverage.gatk(paste(cancerPrefix, suffix, sep=""))

# calculate log coverage ratio
logR = calculate.logR(normal, tumor)

# call CNV results for exons
chr.list<-as.character(unique(normal$chr))
print(chr.list)
# eCNVResults = c()
# for (i in length(chr.list)) {
#   idx = (normal$chr == chr.list[i])
#   ecnv = classify.eCNV(normal=normal[idx,], tumor=tumor[idx,], logR=logR[idx], min.spec=0.9999, min.sens=0.9999, option="spec", c=0.5, l=readlength)
#   eCNVResults = rbind(eCNVResults, ecnv)
# }
eCNVResults<-do.call(rbind, lapply(chr.list, function(chrom) 
  classify.eCNV(normal=normal[normal$chr==chrom,], tumor=tumor[normal$chr==chrom,], logR=logR[normal$chr==chrom],
                min.spec=0.9999, min.sens=0.9999, option="spec", c=0.5, l=readlength)))

# combine exon cnv results into segments
demo.cnv = multi.CNV.analyze(normal, tumor, logR=logR, all.cnv.ls=list(eCNVResults), coverage.cutoff=5, min.spec=0.99, min.sens=0.99, option="auc", c=0.5)

do.plot.eCNV(demo.cnv, lim.quantile=0.99, style="bp", bg.cnv=eCNVResults, line.plot=T)
write.output(eCNVResults, demo.cnv, outputFile)