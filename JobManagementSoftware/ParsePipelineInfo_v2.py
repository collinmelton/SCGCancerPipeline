import csv
from Job_v2 import pipelineJob
from GridEngine import GridEngine

def getVariablesDict(variablesFile):
    f = open(variablesFile, "rU")
    reader = csv.reader(f)
    result = {}
    header=reader.next()
    for row in reader:
        result[row[0]]=row[1]
    #print result
    return result

# eventually this should parse inputs from the JobInfoFile, to get information for each job,
# right now I just return a dictionary object directly as a placeholder
def parsePipelineInfoFile(filename, variablesFile, grid):
    # replace global variable names
    variableReplacements=getVariablesDict(variablesFile)
    f = open(filename, "rU")
    data=f.read()
    f.close()
    for key in variableReplacements.keys():
        data=data.replace(key, variableReplacements[key])
    f = open(filename+".temp", "w")
    f.write(data)
    f.close()
    
    # read filename with variables replaced
    f = open(filename+".temp", "rU")
    reader = csv.reader(f)
    header=reader.next()
    
    # read in jobs
    Jobs = {} # make job dict
    dependencyDict={} # keep track of dependencies after multiplicity expansion
    for row in reader:
        # read job info into newSubTask dict
        newSubTask={}
        for i in range(0, min(len(header), len(row))):
            newSubTask[header[i]]=row[i]
        newSubTask["scriptName"]=("pat_"+newSubTask["patientID"]+"_"+newSubTask["scriptName"]).replace(" ", "_") # add patient id to scriptname to uniquely identify
        dependencyDict[newSubTask["scriptName"]]=[] # make list for dependency names after splitting 
        
        # get jobs for multiplicity var
        if newSubTask["run"]=="TRUE":
            newSubtasks=[]
            multiplicity=newSubTask["multiplicity"].split("|")
            for m in multiplicity:
                # create new subtask for each multiplicity
                newMSubtask={}
                for key in newSubTask.keys():
                    newMSubtask[key]=newSubTask[key].replace("$MULTIPLICITYVAR_FORFILE", m.replace(" ", "_")).replace("$MULTIPLICITYVAR", m).replace("$OUTPUTPATH", newSubTask["outputPath"]).replace("$PATIENTID", newSubTask["patientID"]).replace("\\n", "\n")
                newMSubtask["scriptName"]=(newMSubtask["scriptName"]+"_"+m).replace(" ", "_")
                Jobs[newMSubtask["scriptName"]]=pipelineJob(grid,newMSubtask) 
                dependencyDict[newSubTask["scriptName"]].append(newMSubtask["scriptName"])
    for job in Jobs:
        Jobs[job].setDependencies(dependencyDict)
    f.close()
    #print "done"
    return Jobs

if __name__ == "__main__":
    grid = GridEngine("test")
    getVariablesDict("/Users/cmelton/Documents/Aptana Studio 3 Workspace/CancerPipeline/Variables.csv")
    Jobs= parsePipelineInfoFile("/Users/cmelton/Documents/Aptana Studio 3 Workspace/CancerPipelineSharedCode/JobManagementSoftware/v2jobfiletest.csv", "/Users/cmelton/Documents/Aptana Studio 3 Workspace/CancerPipelineSharedCode/JobManagementSoftware/Variables.csv", grid)
    for job in Jobs:
        print Jobs[job].name