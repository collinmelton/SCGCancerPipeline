
import os


def parseFile(filename, outfile, header, stopline="\n"):
    f=open(filename, 'r')
    g=open(outfile, 'w')
    line=f.readline()
    start=False
    while line!="":
        if header in line: 
            g.write("\t".join(line.strip("\n").split()))
            start=True
        elif line==stopline: 
            start=False
        elif start:
            g.write("\n"+"\t".join(line.strip("\n").split()))
        line=f.readline()
    f.close()
    g.close()
    
dir_name="/Users/cmelton/Documents/Lab/SnyderLab/CancerPipeline/RecalDataAnalysis/"
for file_name in os.listdir(dir_name):
    if ".grp" in file_name:
        parseFile(os.path.join(dir_name, file_name), 
                  os.path.join(dir_name, file_name+".RecalTable2"),
                  "#:GATKTable:RecalTable2:")