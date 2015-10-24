import csv
from collections import OrderedDict

def getVariablesDict(variablesFile):
    f = open(variablesFile, "rU")
    reader = csv.reader(f)
    i = 1
    result = {}
    header=reader.next()
    for row in reader:
        result[row[0]]=row[1]
    #print result
    return result

# eventually this should parse inputs from the JobInfoFile, to get information for each job,
# right now I just return a dictionary object directly as a placeholder
def parsePipelineInfoFile(filename, variablesFile):
    variableReplacements=getVariablesDict(variablesFile)
    f = open(filename, "rU")
    data=f.read()
    f.close()
    for key in variableReplacements.keys():
        data=data.replace(key, variableReplacements[key])
    f = open(filename+".temp", "w")
    f.write(data)
    f.close()
    f = open(filename+".temp", "rU")
    reader = csv.reader(f)
    i = 1
    result = OrderedDict({})
    header=reader.next()
    #print header
    for row in reader:
        newSubTask={}
        for i in range(0, min(len(header), len(row))):
            newSubTask[header[i]]=row[i]
        if newSubTask["run"]=="TRUE":
            newSubtasks=[]
            multiplicity=newSubTask["multiplicity"].split("|")
            for m in multiplicity:
                newMSubtask={}
                for key in newSubTask.keys():
                    newMSubtask[key]=newSubTask[key].replace("$MULTIPLICITYVAR_FORFILE", m.replace(" ", "_")).replace("$MULTIPLICITYVAR", m).replace("$OUTPUTPATH", newSubTask["outputPath"]).replace("$PATIENTID", newSubTask["patientID"]).replace("\\n", "\n")
                newSubtasks.append(newMSubtask)
            if newSubTask["patientID"] not in result.keys():
                result[newSubTask["patientID"]]=OrderedDict({})
            if newSubTask["subTaskGroup"] not in result[newSubTask["patientID"]].keys():
                result[newSubTask["patientID"]][newSubTask["subTaskGroup"]]=[]
            for newSubTask in newSubtasks:
                result[newSubTask["patientID"]][newSubTask["subTaskGroup"]].append(newSubTask)
    f.close()
    #print "done"
    return result

if __name__ == "__main__":
    getVariablesDict("/Users/cmelton/Documents/Aptana Studio 3 Workspace/CancerPipeline/Variables.csv")
    x= parsePipelineInfoFile("/Users/cmelton/Documents/Aptana Studio 3 Workspace/CancerPipeline/PersonalisBrca121713_111513_ThroughRecal_2.csv", "/Users/cmelton/Documents/Aptana Studio 3 Workspace/CancerPipeline/Variables.csv")
    f = open("outputTest.tsv", "w")
    for group in x['1'].keys():
        for task in x["1"][group]:
            f.write("\n"+task["scriptCommand"]+"\t"+task["inputs"])
    f.close()