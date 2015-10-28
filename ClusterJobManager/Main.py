

#
# This is potential start code for the Cancer Pipeline.
#

import time, drmaa ## time has a sleep function, Queue is a Queue for multithreaded applications, threading has a Thread superclass that is utilized here, drmaa has the code for interacting with the grid engine
from optparse import OptionParser ## this is code to help with option parsing
import Job ## this is a subclass of a Job superclass that is used to run any Job
import LogFile ## this is a class to write to a log file
from GridEngine import GridEngine
import os
from ParsePipelineInfo import parsePipelineInfoFile
import pickle

# determine if any jobs remain
def JobsRemain(JobsDict):
    for jobName in JobsDict:
        job=JobsDict[jobName]
        if (not job.started) and job.dependenciesOkay(JobsDict):
            print job.name, job.started
            return True
    return False

# # update jobs dict with all finished jobs
# def UpdateJobsDict(JobsDict):
#     for jobName in JobsDict:
#         job=JobsDict[jobName]
#         job.updateStatus()

def CheckForNewJobs(JobsDict, grid, csvFolder, jobinfofiles, variablesFile, JobVariables):
    files =  os.listdir(csvFolder)
    newfiles = [f for f in files if (".csv" in f and f not in jobinfofiles and ".temp" not in f)]
    jobinfofiles=list(set(jobinfofiles+newfiles))
    for f in newfiles:
        try: 
            newDict=parsePipelineInfoFile(os.path.join(csvFolder, f), variablesFile, grid, JobVariables=JobVariables)
            for key in newDict: JobsDict[key]=newDict[key]
        except: print "error parsing file:", f
    return JobsDict, jobinfofiles

# writes status to dictionary in pickle file
def WriteStatus(JobsDict, statusFolder):
    StatusDict={}
    for job in JobsDict:
        StatusDict[job]=JobsDict[job].status
    pickle.dump(StatusDict, open(os.path.join(statusFolder, "status.pickle"), "w"))

# run all jobs until none remain
def RunJobs(JobsDict, grid, csvFolder, statusFolder, variablesFile, JobVariables):
    jobinfofiles=[]
    # forever check if jobs are ready and start them
    while True:
        # add new jobs if they exist
        JobsDict, jobinfofiles = CheckForNewJobs(JobsDict, grid, csvFolder, jobinfofiles, variablesFile, JobVariables)
        # update Jobs Dict with newly completed Jobs
        grid.updateJobDict(JobsDict)
        # start jobs that are ready
        for jobName in JobsDict:
            job=JobsDict[jobName]
            if job.readyToRun(JobsDict) and not job.started:
                job.start()
                print job.name, job.started
        # write status to pickle file
        WriteStatus(JobsDict, statusFolder)
        # wait before checking on jobs to run
        time.sleep(30)
        
# this functions gets the command line options for running the program
def getOptions():
    parser = OptionParser()
    parser.add_option("--V", dest="VariablesFile", help = "this file contains key value pairs for variable to be replaced in the JobInfoFile",
                      metavar = "FILE", type = "string", default = "/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/JobManagementSoftware/Variables.csv")
    parser.add_option("--L", dest = "logFileFilePath", help = "path/name of log file to be output"+
                      "concurrently", metavar = "FILE", default = "logFileTest.txt", type = "string")
    parser.add_option("--C", dest = "csvFolder", help = "", metavar = "FILE", default = "/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/CODESPatientData/JobCSVFiles", type = "string")
    parser.add_option("--IF", dest = "inputFolder", help = "", metavar = "FILE", default = "/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/CODESPatientData/InputFiles", type = "string")
    parser.add_option("--OF", dest = "outputFolder", help = "", metavar = "FILE", default = "/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/CODESPatientData/OutputFiles", type = "string")
    parser.add_option("--SF", dest = "statusFolder", help = "", metavar = "FILE", default = "/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/CODESPatientData/StatusFiles", type = "string")
    parser.add_option("--JEF", dest = "jobErrorFolder", help = "", metavar = "FILE", default = "/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/CODESPatientData/JobErrorFiles", type = "string")
    parser.add_option("--JOF", dest = "jobOutputFolder", help = "", metavar = "FILE", default = "/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/CODESPatientData/JobOutputFiles", type = "string")
    parser.add_option("--JSF", dest = "jobScriptFolder", help = "", metavar = "FILE", default = "/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/CODESPatientData/JobScriptFiles", type = "string")
    
    (options, args) = parser.parse_args()
    return options

# this is the main function
def run():
    # get options and defaults
    options = getOptions()
    
    grid = GridEngine(options.logFileFilePath)
    
    try:
        # create Grid Engine class
        print "reading job info"
        
        # get a dictionary of class jobs for all jobs to be run
        JobsDict={}
    
        JobVariables={"$InputsFolder":options.inputFolder, "$OutputsFolder": options.outputFolder,  
                      "$JobOutputFolder":options.jobOutputFolder, "$JobErrorFolder":options.jobErrorFolder, 
                      "$ScriptPath": options.jobScriptFolder}
    
        # run all jobs        
        RunJobs(JobsDict, grid, options.csvFolder, options.statusFolder, options.VariablesFile, JobVariables)

    finally:
        # exit drmaa session
        if grid:
            grid.exit()

run()