'''
Created on Oct 29, 2015

@author: cmelton
'''
from optparse import OptionParser ## this is code to help with option parsing

def getOptions():
    parser = OptionParser()
    parser.add_option("--P", dest = "prefix", help = "",
                      metavar = "STRING", type = "string", default = "./test/test")
    parser.add_option("--M", dest = "multiplicityVar", help = "",
                      metavar = "STRING", type = "string", default = "1|2|3")
    parser.add_option("--FQ1", dest = "FASTQ1", help = "",
                      metavar = "FILE", default = "", type = "string")
    parser.add_option("--FQ2", dest = "FASTQ2", help = "",
                      metavar = "FILE", default = "", type = "string")
    parser.add_option("--T", dest = "test", help = "",
                      metavar = "FILE", default = "F", type = "string")
    (options, args) = parser.parse_args()
    return options


def readChunk(f, lines=4):
    result = ""
    for i in range(lines):
        result += f.readline()
    return result

def readChunks(fq1, fq2, lines=4):
    chunk1 = readChunk(fq1, lines=lines)
    chunk2 = readChunk(fq2, lines=lines)
    finished = (len(chunk1.split("\n"))<lines and len(chunk1.split("\n"))<lines)
    return finished, chunk1, chunk2

def splitFiles(fastq1, fastq2, filePrefix, multiplicity):
    # open original fastq files
    fq1 = open(fastq1, 'r')
    fq2 = open(fastq2, 'r')
    # generate dict of new fastqs
    newFastqs = [{"R1": open(filePrefix+"_"+m+"_R1.fastq", 'w'), "R2": open(filePrefix+"_"+m+"_R2.fastq", 'w')} for m in multiplicity]
    # while files still not over iterate over new fastqs and add reads
    finished = False
    while not finished:
        for fqDict in newFastqs:
            if not finished:
                finished, chunk1, chunk2 = readChunks(fq1, fq2)
            if not finished:
                fqDict["R1"].write(chunk1)
                fqDict["R2"].write(chunk2)
    # close all files
    fq1.close()
    fq2.close()
    for fqDict in newFastqs:
        fqDict["R1"].close()
        fqDict["R2"].close()

def writeTestFastqs():
    fnames = ["./test/fastq_test_R1.fastq", "./test/fastq_test_R2.fastq"]
    for fname in fnames:
        f = open(fname, 'w')
        f.write("\n".join([str(i) for i in range(100)]))
        f.close()
        print "wrote file", fname
    return fnames[0], fnames[1]

def run():
    options = getOptions()
    filePrefix = options.prefix
    multiplicity = options.multiplicityVar.split("|")
    
    if options.test=="T": 
        fastq1, fastq2 = writeTestFastqs()
    else:
        fastq1 = options.FASTQ1
        fastq2 = options.FASTQ2
        
    splitFiles(fastq1, fastq2, filePrefix, multiplicity)


if __name__ == '__main__':
    run()