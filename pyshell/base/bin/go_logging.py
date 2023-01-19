
"""
GO_Logging - logging tools
"""

def defaultLogger( name=None ) :

    """
    Return simple logger to standard output.
    """
    
    # external:
    import logging
    import sys
    
    # new logger instance, or receive the logger with the provided name:
    logger = logging.getLogger(name)
    
    # no handlers yet ?
    if len(logger.handlers) == 0 :
        # add handler to standard output:
        logger.addHandler(logging.StreamHandler(sys.stdout))
        # let logger print all messages:
        logger.setLevel(0)
    #endif
    
    # ok
    return logger
    
#enddef

