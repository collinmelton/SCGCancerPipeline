

#
# This is potential start code for the Cancer Pipeline.
#

import time, Queue, threading, drmaa ## time has a sleep function, Queue is a Queue for multithreaded applications, threading has a Thread superclass that is utilized here, drmaa has the code for interacting with the grid engine
from optparse import OptionParser ## this is code to help with option parsing
import Job ## this is a subclass of a Job superclass that is used to run any Job
import LogFile ## this is a class to write to a log file
import ParsePipelineInfo ## this will be a class that parses a text input to get parameters for the run

# this class is where we define the pipeline, an example (jobA) is currently run to demonstrate that this works
class ThreadRunPipeline(threading.Thread):
    def __init__(self, queue, drmaaSession, logFileWriter):
        threading.Thread.__init__(self)
        self.queue = queue
        self.drmaaSession = drmaaSession 
        self.logFileWriter = logFileWriter
        self.parallelDicts=[]
    
    # this method runs a job for a particular key in the pipelineDict
    def runThreadAndWait(self, jobThread, jobEvent):
        # start job thread
        jobThread.setDaemon(True)
        jobThread.start()

        # wait for thread to finish
        jobEvent.wait()
        return (not jobThread.checkError())
    
    # this method runs through the pipeline for each job dictionary in the queue
    def run(self):
        while True:
            patientAllTasksDict = self.queue.get()
            patientError=False
            subtaskList=map(lambda x: str(x), sorted(map(lambda x: int(x), patientAllTasksDict.keys())))
            for subTaskGroup in subtaskList:
                if patientError: break
                subTaskGroupDict=patientAllTasksDict[subTaskGroup]
                for JobDict in subTaskGroupDict:
                    self.parallelDicts.append(JobDict)
                jobEvent = threading.Event()
                jobThread = Job.pipelineJob(jobEvent, self.drmaaSession, self.logFileWriter, self.parallelDicts)
                if not self.runThreadAndWait(jobThread, jobEvent): 
                    patientError=True
                    continue # break out of current iteration of for loop because error occurred in job
                self.parallelDicts=[]
            self.queue.task_done()

# this functions gets the command line options for running the program
def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest = "JobInfoFile", help = "this should be a file with info about the job",
                      metavar = "FILE", type = "string", default = "PersonalisBrca121713_111513_ThroughRecal_2.csv")#"BRCATest.csv")#"/home/cmelton/apps/CancerPipeline/601453_NA5K5CHE_Pipeline.csv")#PersonalisBrca121713_111513_ThroughRecal_2.csv")#/home/cmelton/apps/CancerPipeline/601453_NA5K5CHE_Pipeline.csv")#601453_NA5K5CHE_Pipeline.csv")
    parser.add_option("--T", dest = "numberOfThreads", help = "number of tumor/normal pairs to process "+
                      "concurrently",
                      metavar = "INTEGER", default = 3, type = "int")
    parser.add_option("--V", dest="VariablesFile", help = "this file contains key value pairs for variable to be replaced in the JobInfoFile",
                      metavar = "FILE", type = "string", default = "/srv/gsfs0/projects/gbsc/Clinical/cancerPatientAnno/SharedSoftware/JobManagementSoftware/Variables.csv")
    parser.add_option("--L", dest = "logFileFilePath", help = "path/name of log file to be output"+
                      "concurrently", metavar = "FILE", default = "/home/cmelton/scripts/logFileTest.txt", type = "string")
    (options, args) = parser.parse_args()
    return options

# this is the main function
def run():
    # get options and defaults
    options = getOptions()
    
    # create drmaa session
    drmaaSession = drmaa.Session()
    try:
        # initialize drmaa session (this is session object is our means of interacting with
        # the grid engine
        drmaaSession.initialize()

        # get list of dicts with parameters for our different runs of the pipeline
        print "reading job info"
        pipelineDicts = ParsePipelineInfo.parsePipelineInfoFile(options.JobInfoFile, options.VariablesFile)

        #initialize queue to process one tumor normal set at a time
        print "initializing queue"
        queue = Queue.Queue();
        
        # put all samples in the queue, by passing in the dictionary of parameters for the each
        # particular run of the pipeline
        for key in pipelineDicts:
            queue.put(pipelineDicts[key])
            
        # initialize logfile writer, this is a common instance shared by all the threads
        logFileWriter = LogFile.LogFile(options.logFileFilePath)    
        
        #make a single thread
        for i in range(options.numberOfThreads):
            t = ThreadRunPipeline(queue, drmaaSession, logFileWriter)
            t.setDaemon(True)
            t.start()
            time.sleep(10)
        
        #wait for queue to be done processing
        queue.join()

    finally:
        # exit drmaa session
        if drmaaSession:
            drmaaSession.exit()

run()