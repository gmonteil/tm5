
# ------------------------------------------------
# help
# ------------------------------------------------

"""
GO SubProcess - tools to call a subprocess in combination with logging

Routines

    p = call( *popenargs )

        Call subprocess without logging.
        Returns an object with fields:
          retcode   (integer)
          stdout    (list of strings)
          stderr    (list of strings)
        Raises a CallingError or StatusError on failure.

    p = log_call( *popenargs, logger=None, loglevel=0 )

        Call subprocess. After process is terminated, send the standard output 
        and standard error to the logger. If no logger instance is applied,
        a standard logger to created. The standard output is written with
        loglevel logging.INFO or the level specified as argument ;
        the standard error is written with loggin.DEBUG.
        Returns an object with fields:
          retcode   (integer)
          stdout    (list of strings)
          stderr    (list of strings)
        Raises a CallingError or StatusError on failure.

    p = watch_call( *popenargs, logger=None, loglevel=0 )

        Call subprocess and send standard output to logger while running.
        Standard error is redirected to standard output to maintain the
        order in which lines are written as good as possible.
        See arguments of 'log_call' for the logging arguments.
        Returns an object with fields:
          retcode   (integer)
          stdout    (list of strings)
          stderr    (list of strings; empty)
        Raises a CallingError or StatusError on failure.
        
Exceptions

    CallingError
        Raised when a command could not be started at all.
        Attribute: command

    StatusError
        Raised when a non-zero exit status is received from a command.
        Attributes: command, returncode, stdout, stderr

Example

    # external:
    import logging

    # run command, log outputs:
    p = log_call( ['/bin/ls','-l'] )
    
    # idem, but log outside call:
    try :
        p = call( ['/bin/ls','-l'] )
    except CallingError, err :
        logging.error( err )
        raise Exception
    except StatusError, err :
        for line in err.stderr : logging.error(line)
        logging.error( err )
        raise Exception
    except Exception, err :
        logging.error( 'unknown error : %s' % err )
        raise Exception
    #endtry
    # display output:
    for line in p.stdout : logging.info(line)
    for line in p.stderr : logging.error(line)

"""


# ------------------------------------------------
# exceptions
# ------------------------------------------------

class Error( Exception ):
    """
    Base class for exceptions in this module.
    """
    pass
#endclass

# *

class CallingError( Error ) :
    """
    This exception is raised when an attempt to call a command fails.

    Attributes:
        command       # str or str list
    """

    def __init__( self, command, strerror ) :
        """
        Store exception attributes.
        """
        # store:
        self.command = command
        self.strerror = strerror
    #enddef
    
    def __str__( self ) :
        """
        Format error message.
        """
        # error message:
        return 'Error from calling "%s" : %s' % (self.command,self.strerror)
    #enddef
    
#endclass

# *

class StatusError( Error ) :
    """
    This exception is raised when a called command returns a non-zero exit status.

    Attributes:
        command       # str or str list
        returncode    # integer return status
        stdout        # str list
        stderr        # str list
    """

    def __init__( self, command, returncode, stdout, stderr ) :
        """
        Store exception attributes.
        """
        # store:
        self.command = command
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
    #enddef
    
    def __str__( self ) :
        """
        Format error message.
        """
        # error message:
        return 'Call "%s" returned non-zero status %i' % (str(self.command),self.returncode)
    #enddef
    
#endclass


# ------------------------------------------------
# routines
# ------------------------------------------------

class call( object ) :

    def __init__( self, command, **kwargs ) :

        """
        Call subprocess and return object with returncode, output, and error.
        """

        # external:
        import subprocess
        
        # call suprocess, catch calling errors:
        try :
            p = subprocess.Popen( command, 
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  **kwargs )
        except Exception, err :
            raise CallingError( command, err )
        #endtry

        # wait for process to terminate, receive standard output and error:
        stdout_str,stderr_str = p.communicate()
        
        # remove final end-of-lines:
        stdout_str = stdout_str.strip('\n')
        stderr_str = stderr_str.strip('\n')
        
        # store lists:
        if len(stdout_str) > 0 :
            self.stdout = stdout_str.split('\n')
        else :
            self.stdout = []
        #endif
        if len(stderr_str) > 0 :
            self.stderr = stderr_str.split('\n')
        else :
            self.stderr = []
        #endif
        
        # store return code:
        self.returncode = p.returncode
        
        # error ?
        if self.returncode != 0 : raise StatusError(command,self.returncode,self.stdout,self.stderr)

        # ok
        return

    #enddef

#endclass

# *

class log_call( object ) :

    def __init__( self, command, logger=None, loglevel=None, **kwargs ) :

        """
        Call subprocess and send standard output and error to logger.
        """

        # external:
        import logging
        import sys
        import subprocess
        
        # tools:
        import go_logging

        # setup logger to standard output if not provided as argument:
        if logger == None : logger = go_logging.defaultLogger()

        # logging levels:
        errlevel = logging.ERROR
        outlevel = logging.INFO
        if loglevel != None : outlevel = loglevel
        
        # init outputs:
        self.stdout = []
        self.stderr = []

        # call suprocess, catch calling errors:
        try :
            p = subprocess.Popen( command, 
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  **kwargs )
        except Exception, err :
            logger.error( str(err) )
            raise CallingError( command, err )
        #endtry

        # wait for process to terminate, receive standard output and error:
        stdout_str,stderr_str = p.communicate()
        
        # remove final end-of-lines:
        stdout_str = stdout_str.strip('\n')
        stderr_str = stderr_str.strip('\n')
        
        # store lists:
        if len(stdout_str) > 0 :
            self.stdout = stdout_str.split('\n')
        else :
            self.stdout = []
        #endif
        if len(stderr_str) > 0 :
            self.stderr = stderr_str.split('\n')
        else :
            self.stderr = []
        #endif
        
        # write standard output and error to logger:
        for line in self.stdout : logger.log( outlevel, line )
        for line in self.stderr : logger.log( errlevel, line )

        # store return code:
        self.returncode = p.returncode
        
        # error ?
        if self.returncode != 0 : raise StatusError(command,self.returncode,self.stdout,self.stderr)

        # ok
        return

    #enddef

#endclass

# *

class watch_call( object ) :

    def __init__( self, command, logger=None, loglevel=0, **kwargs ) :

        """
        Call subprocess and send standard output to logger while running.
        Standard error is redirected to standard output to maintain the
        order in which lines are written as good as possible.
        """

        # external:
        import logging
        import sys
        import subprocess
        
        # tools:
        import go_logging

        # setup logger to standard output if not provided as argument:
        if logger == None : logger = go_logging.defaultLogger()

        # logging levels:
        errlevel = logging.ERROR
        outlevel = logging.INFO
        if loglevel != 0 : outlevel = loglevel
        
        # init outputs:
        self.stdout = []
        self.stderr = []

        # call suprocess, catch calling errors:
        try :
            p = subprocess.Popen( command, 
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT,
                                  **kwargs )
        except Exception, err :
            logger.error( str(err) )
            raise CallingError( command, err )
        #endtry

        # wait for process to terminate:
        while p.poll() == None :
            # read next line from standard output:
            line = p.stdout.readline()
            # empty means no data ...
            if len(line) == 0 : continue
            # store:
            self.stderr.append(line.strip('\n'))
            # write:
            logger.log( outlevel, line.strip('\n') )
        #endwhile
        
        # read remaing lines:
        for line in p.stdout.readlines() :
            # store:
            self.stderr.append(line.strip('\n'))
            # write:
            logger.log( outlevel, line.strip('\n') )
        #endwhile

        # store return code:
        self.returncode = p.returncode
        
        # error ?
        if self.returncode != 0 : raise StatusError(command,self.returncode,self.stdout,self.stderr)

        # ok
        return

    #enddef

#endclass


# ------------------------------------------------
# end
# ------------------------------------------------
