
INFO    = 'info'
WARNING = 'warning'
ERROR   = 'error'
OFF     = 'off'

logLevel = 'info'

def info(s):
    if logLevel == 'info':
        print '[INFO]: ' + s

def warning(s):
    if logLevel == 'warning':
        print '[WARNING]: ' + s

def error(s):
    if logLevel == 'error':
        print '[ERROR]: ' + s

def infoVar(v, name):
    info(name + ": " + str(v))
