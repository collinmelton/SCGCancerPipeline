

#
# This is potential start code for the Cancer Pipeline.
#

import time, drmaa ## time has a sleep function, Queue is a Queue for multithreaded applications, threading has a Thread superclass that is utilized here, drmaa has the code for interacting with the grid engine
from optparse import OptionParser ## this is code to help with option parsing
import Job_v2 ## this is a subclass of a Job superclass that is used to run any Job
import LogFile ## this is a class to write to a log file
from ParsePipelineInfo_v2 import parsePipelineInfoFile ## parses a text input to get parameters for the jobs
from GridEngine import GridEngine

# determine if any jobs remain
def JobsRemain(JobsDict):
    for jobName in JobsDict:
        job=JobsDict[jobName]
        if (not job.started) and job.dependenciesOkay(JobsDict):
            print job.name, job.started
            return True
    return False

# update jobs dict with all finished jobs
def UpdateJobsDict(JobsDict):
    for jobName in JobsDict:
        job=JobsDict[jobName]
        job.updateStatus()

# run all jobs until none remain
def RunJobs(JobsDict, grid):
    # while jobs remain check if jobs are ready and start them
    while (JobsRemain(JobsDict)):
        # update Jobs Dict with newly completed Jobs
        grid.updateJobDict(JobsDict)
        # start jobs that are ready
        for jobName in JobsDict:
            job=JobsDict[jobName]
            if job.readyToRun(JobsDict) and not job.started:
                job.start()
                print job.name, job.started
        # wait before checking on jobs to run
        time.sleep(60)
        
# this functions gets the command line options for running the program
def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest = "JobInfoFile", help = "this should be a file with info about the job",
                      metavar = "FILE",
                       type = "string", default = "v2jobfiletest.csv")
    parser.add_option("--V", dest="VariablesFile", help = "this file contains key value pairs for variable to be replaced in the JobInfoFile",
                      metavar = "FILE", type = "string", default = "/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/JobManagementSoftware/Variables.csv")
    parser.add_option("--L", dest = "logFileFilePath", help = "path/name of log file to be output"+
                      "concurrently", metavar = "FILE", default = "logFileTest.txt", type = "string")
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
        JobsDict=parsePipelineInfoFile(options.JobInfoFile, options.VariablesFile, grid)

        # run all jobs        
        RunJobs(JobsDict, grid)

    finally:
        # exit drmaa session
        if grid:
            grid.exit()

run()