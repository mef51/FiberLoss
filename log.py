
INFO    = 0
WARNING = 1
ERROR   = 2
OFF     = 3

logLevel = INFO

def info(s):
    if logLevel <= INFO:
        print '[INFO]: ' + s

def warning(s):
    if logLevel <= WARNING:
        print '[WARNING]: ' + s

def error(s):
    if logLevel <= ERROR:
        print '[ERROR]: ' + s

def infoVar(v, name):
    info(name + ": " + str(v))
