import threading, drmaa

# This is the common super class to create and start jobs
class Job(threading.Thread):
    def __init__(self, jobEvent, drmaaSession, logFileWriter, jobInfoDicts):
        threading.Thread.__init__(self)
        self.jobEvent = jobEvent
        self.drmaaSession = drmaaSession
        self.logFileWriter = logFileWriter
        self.jobInfoDicts = jobInfoDicts
        self.jobList = []
        self.error=False
        self.jobTemplate=None
        
    # write to the common log file
    def writeToLogFile(self, textToWrite):
        self.logFileWriter.write(textToWrite)
        
    # this function checks for an error
    def checkError(self):
        return self.error
    
    # this function helps write a generic shell script for qsub
    def writeShellScript(self, filePath, scriptName, runTime, error, output, customizations, executionLine, memory, emailAddress):
        newline = "\n"
        f = open(filePath, "w")
        f.write("#!/bin/sh" + newline)
        # set the name of the job
        f.write("#$ -N " + scriptName + newline)
        # set max memory usage per slot
        f.write("#$ -l h_vmem=" + str(memory) + "G" + newline)
        # set max run time
        f.write("#$ -l h_rt="+runTime+ newline)    
        # send mail when job ends or aborts
        f.write("#$ -m ea" + newline)
        # specify email address
        f.write("#$ -M " + emailAddress + newline)
        # check for errors in job submission
        f.write("#$ -w e" + newline)
        # specify output and error files= names
        f.write("#$ -o " + output + newline)
        f.write("#$ -e " + error + newline)
        # add custom qsub options
        f.write(customizations + newline)
        # add code to be run
        f.write(executionLine)
        f.close()     
    
    # this method generates a jobtemplate object to be run
    # script path is the path to the script .sh .py or whatever, that is to be run, 
    # args is a list of arguments for the script
    def createJob(self, scriptPath, args):
        self.jobTemplate = self.drmaaSession.createJobTemplate()
        self.jobTemplate.remoteCommand = scriptPath
        self.jobTemplate.args = args
        # the -b no option tells sge to not read the command as a binary file,
        # this means all the sge commented out options will be read, default is -b yes
        self.jobTemplate.nativeSpecification = "-b no"
        
    # this method rus the jobtemplate object
    def runJob(self):
#         print self.jobTemplate.args
#         print self.jobTemplate.remoteCommand
#         print self.jobTemplate.nativeSpecification
        jobid = self.drmaaSession.runJob(self.jobTemplate)
        self.jobList.append(jobid)

    # this method waits for all running jobs to finish
    def synchronize(self):
        self.drmaaSession.synchronize(self.jobList, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
        for curjob in self.jobList:
            retval = self.drmaaSession.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
            print "return value", retval.exitStatus
            if not self.error: self.error=(retval.exitStatus!=0)
        
            
    # this method marks the jobevent complete
    def markComplete(self):
        self.drmaaSession.deleteJobTemplate(self.jobTemplate)
        self.jobEvent.set()