from optparse import OptionParser

# this functions gets the command line options for running the program
def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest = "input1", help = "this is the input cancer file",
                      metavar = "FILE", type = "string")
    parser.add_option("--J", dest = "output1", help = "this is the output cancer file",
                      metavar = "FILE", type = "string")
    parser.add_option("--K", dest = "input2", help = "this is the input normal file",
                      metavar = "FILE", type = "string")
    parser.add_option("--L", dest = "output2", help = "this is the output normal file",
                      metavar = "FILE", type = "string")
    parser.add_option("--F", dest = "outputFileName", help = "this is the output map file name",
                      metavar = "FILE", type = "string")
    (options, args) = parser.parse_args()
    return options

# this is the main function
def run():
    # get options and defaults
    options = getOptions()
    f =open(options.outputFileName, "w")
    f.write(options.input1+"\t"+options.output1+"\n"+options.input2+"\t"+options.output2)
    f.close()

run()