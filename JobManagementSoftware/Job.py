import threading, drmaa, time
import GridEngineJobThread

# this is a class for running a single job in the pipeline, we would write many of these,
# one for each distinct job/component of the pipeline, it is subclassed off of Job which
# has some nice functions that allow you to write a shell script interpretable by the grid engine,
# each of us could write one of these for each component of the pipeline we are responsible for
class pipelineJob(GridEngineJobThread.Job):
    def __init__(self, jobEvent, drmaaSession, logFileWriter, jobInfoDicts):
        super(pipelineJob, self).__init__(jobEvent, drmaaSession, logFileWriter, jobInfoDicts)
    
    def start(self):
        for jobInfoDict in self.jobInfoDicts:
            # write a shell script to run
            self.writeShellScript(jobInfoDict["scriptPath"], jobInfoDict["scriptName"], 
                                  jobInfoDict["scriptTime"], jobInfoDict["scriptErrorFileDirectory"],
                                  jobInfoDict["scriptOutputFileDirectory"], jobInfoDict["scriptCustomizations"],
                                  jobInfoDict["scriptCommand"], jobInfoDict["scriptMemory"], 
                                  jobInfoDict["scriptEmailAddress"])
            # create job
            self.createJob(jobInfoDict["scriptPath"], jobInfoDict["inputs"].split("|"))
            
            # run job
            self.runJob()
        
        # synchronize jobs
        self.synchronize()
        # mark thread complete
        for jobInfoDict in self.jobInfoDicts:
            self.writeToLogFile(jobInfoDict["patientID"]+"\t"+jobInfoDict["subTaskGroup"]+
                            "\t"+jobInfoDict["scriptName"]+"\t"+"job complete, error status: "+str(self.error))
        self.markComplete()
        print "finished job"