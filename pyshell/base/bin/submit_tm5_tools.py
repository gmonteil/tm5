#
# TM5 run tools
#


# ***


def Command_Line( rcf, exe, args, in_debugger ) :

    """
    Return command line.

    ARGUMENTS
        rcf
                Rcfile with settings.
        exe
                Name of executable.
        args
                Arguments to be passed to executable.
        in_debugger
                Set to 'True' if the job should be run in a debugger.

    RETURN VALUES
        cmndline
                Command line to be executed.
    """

    # external
    import socket
    import logging

    # mpi run ?
    if rcf.get('par.mpi','bool') :

        # number of mpi tasks:
        ntask = rcf.get('par.ntask','int')

        # get command line:
        cmnd_exec = rcf.get('mpirun.command')
        cmnd_args = rcf.get('mpirun.args'   )

        # write command file ?
        cmdfile = rcf.get('mpirun.cmdfile',default='')
        if len(cmdfile) > 0 :
            # write command line for each task:
            f = open(cmdfile,'w')
            for i in range(ntask) : f.write( '%s %s\n' % (exe,args) )
            f.close()
        else :
            # otherwise, add the executable and its arguments:
            cmnd_args = '%s %s %s' % (cmnd_args,exe,args)
        #endif

        # write host file ?
        hostfile = rcf.get('mpirun.hostfile',default='')
        if len(hostfile) > 0 :
            # get hostname:
            hname = socket.gethostname()
            # write hostname for each task:
            f = open(hostfile,'w')
            for i in range(ntask) : f.write( '%s\n' % hname )
            f.close()
        #endif

    else :

        # standard run:
        cmnd_exec = exe
        cmnd_args = args

    #endif

    # run in debugger ?
    if in_debugger :

        # debugger type:
        debugger = rcf.get( 'debugger' )
        # get debugger command:
        debugger_call = rcf.get( 'debugger.command' )
        # large differences ...
        if debugger == 'totalview' :
            # syntaxis: totalview <executable> [-a <arguments>]
            # pass executable:
            cmndline = '%s %s' % (debugger_call,cmnd_exec)
            # add arguments ?
            if len(cmnd_args) > 0 :
                cmndline = '%s -a %s' % (cmndline,cmnd_args)
            #endif
        elif debugger == 'idb' :
            # syntaxis: idb [-args <executable> <arguments>]
            # fill executable and arguments:
            cmndline = '%s -args %s %s' % (debugger_call,cmnd_exec,cmnd_args)
        else :
            logging.error('unsupported debugger : %s' % debugger )
            raise Exception
        #endif

    else :

        # standard line:
        cmndline = '%s %s' % (cmnd_exec,cmnd_args)

    #endif


    #endif

    # ok
    return cmndline

#endif


# ***


def WriteAndSubmitNewJob( rcfile, bindir ) :

    """
    Write first or next rcfile and job files(s) in the job chain;
    if chain is not finished yet, submit a new job.

    The argument could be the name of an rcfile or an rcfile object itself,
    since the submit scrip might have changed some values given
    the provided command line arguments.

    The following function is used:

      submit_tm5_setup_rcfile.WriteRcfile         # writes the new rcfile

    This is placed in a seperate file since users might need to
    change these routines for their specific projects.
    """

    # external:
    import sys
    import logging
    from pyshell.base.bin import rc

    # import setup module:
    from pyshell.base.bin import submit_tm5_setup_rcfile

    # name provided ?
    if type(rcfile) == str :
        # load:
        rcf = rc.RcFile( rcfile )
    else :
        # just copy ..
        rcf = rcfile
    #endif

    # write next rfile, return name:
    try :
        rcfile_next = submit_tm5_setup_rcfile.WriteRcfile( rcf )
    except :
        logging.error( sys.exc_info()[1] )
        logging.error( 'exception from WriteRcfile' )
        raise Exception
    #endtry

    # finished ?
    if rcfile_next == None :
        # info ...
        logging.info( '  end of job chain !' )
    else :
        # write job file(s) for this period and return the (first) name;
        # last command in a file should submit the next job if necessary:
        logging.info( '  write jobfile for %s ...' % rcfile_next )
        jobfile_next = WriteJob( rcfile_next, bindir )
        logging.info( '  submit next job : %s' % jobfile_next )
        try :
            SubmitJob( jobfile_next, rcfile_next )
        except :
            logging.error( sys.exc_info()[1] )
            logging.error( 'exception from SubmitJob' )
            raise Exception
        #endtry
    #endif

    # ok
    return

#endif


# ***


def WriteJob( rcfile, bindir ) :

    """
    jobfile = WriteJob(rcfile)
    Write job file given the settings in rcfile.
    The name of the jobfile is based on the name of the rcfile.
    The last command in the job should submit the next job,
    and the script is therefore written in python.
    """

    # external:
    import os
    from pyshell.base.bin import rc

    # load settings:
    rcf = rc.RcFile( rcfile )

    # basename for scripts etc is name of rcfile minus extension:
    bname, ext = os.path.splitext(rcfile)

    # which shell ?
    job_shell = '/usr/bin/env python'

    steps = rcf.get('job.steps').split(' ')
    nstep = len(steps)

    # loop over job steps:
    for istep, step in enumerate(steps) :

        # next:
        if istep < nstep-1 : 
            step_next = steps[istep+1]

        # list with queue option lines for this step:
        qopt_step = QueueOptions( bname, rcf, step )

        # call to acutal script:
        if step == 'run' :
            # get command line:
            exe  = os.path.join( os.curdir, rcf.get('job.step.%s.exe' % step) )
            args = rcfile
            indb = rcf.get('submit.debugger','bool')
            cmndline = Command_Line( rcf, exe, args, indb )
            # <script> <commandline>
            step_command = '["submit_tm5_step_%s","%s"]' % (step, cmndline)
        else :
            # <script> <rcfile>
            step_command = '["submit_tm5_step_%s","%s","--bindir=%s"]' % (step, rcfile, bindir)

        # name of step job to be written:
        step_job_template = bname+'_%s.jb'

        # actual name:
        step_job = step_job_template % step
        # open:
        f = open( step_job, 'w' )

        # write header:
        f.write( '#! %s\n' % job_shell )
        f.write( '\n' )

        # add queue options:
        for line in qopt_step : 
            f.write(line)

        # add lines to call the actual script:
        f.write( 'import sys\n' )
        f.write( 'import os\n' )
        f.write( 'import logging\n' )
        f.write( 'import subprocess\n' )
        f.write( 'from pyshell.base.bin import rc\n' )
        f.write( 'from pyshell.base.bin import submit_tm5_tools\n' )
        f.write( 'os.chdir("%s")\n' % os.getcwd() )
        f.write( 'retcode = subprocess.call( %s )\n' % step_command )
        f.write( 'if retcode != 0:\n' )
        f.write( '    logging.error( sys.exc_info()[1] )\n' )
        f.write( '    logging.error( \'exception from subprocess call to : %s\' )\n' % step_command )
        f.write( '    sys.exit(1)\n' )

        # add submission of next step?
        if istep < nstep-1 :
            # job script of next step:
            step_job_next = step_job_template % step_next
            # add submission command:
            f.write( '# submit next step:\n' )
            f.write( 'try :\n' )
            f.write( '    submit_tm5_tools.SubmitJob( "%s", "%s" )\n' % (step_job_next,rcfile)  )
            f.write( 'except:\n' )
            f.write( '    logging.error( sys.exc_info()[1] )\n' )
            f.write( '    logging.error( \'exception from SubmitJob( "%s", "%s" )\' )\n' % (step_job_next,rcfile) )
            f.write( '    sys.exit(1)\n' )
            f.write( '#endtry\n' )
            f.write( '\n' )
        else :
            # last step; might be necessary to submit a new job:
            f.write( '# write and submit next job if necessary:\n' )
            f.write( 'submit_tm5_tools.WriteAndSubmitNewJob( "%s", "%s" )\n' % (rcfile,bindir) )
            f.write( '\n' )
        #endif
        f.write( '# info ...\n' )
        f.write( 'logging.info( "end" )\n' )
        f.write( '\n' )

        # close:
        f.close()

        # make it executable and readible for all, writable for user only:
        #                   u+r    u+w    u+x    g+r    g-w    g+x    o+r    o-w    o+x
        os.chmod( step_job, 2**8 + 2**7 + 2**6 + 2**5 + 0    + 2**3 + 2**2 + 0    + 2**0 )

        # fill return value:
        if istep == 0 : 
            jobfile = step_job

    return jobfile


def SubmitCommand( jobfile, rcfile ) :

    """
    Return submit command.
    """

    # external:
    from pyshell.base.bin import rc

    # read settings:
    rcf = rc.RcFile(rcfile)

    # submit to where ?
    if rcf.get('submit.to') == 'foreground' :

        # fill ...
        submitcommand = SubmitCommand_Foreground( jobfile, rcfile )

    elif rcf.get('submit.to') == 'background' :

        # fill ...
        submitcommand = SubmitCommand_Background( jobfile, rcfile )

    elif rcf.get('submit.to') == 'queue' :

        # queue type:
        queue = rcf.get('queue')

        # different options and commands:
        if queue == 'loadleveler' :
            submitcommand = SubmitCommand_LoadLeveler( jobfile, rcfile )
        elif queue == 'bsub' :
            submitcommand = SubmitCommand_BSub( jobfile, rcfile )
        elif queue == 'qsub' :
            submitcommand = SubmitCommand_QSub( jobfile, rcfile )
        else :
            # not supported ...
            logging.error( 'unsupported queue : %s' % queue )
            raise Exception
        #endif

    else :

        # not supported ...
        logging.error( 'unsupported run environment : %s' % rcf.get('submit.to') )
        raise Exception

    #endif

    # ok
    return submitcmnd

#enddef


# ***


def QueueOptions( bname, rcf, step ) :

    """
    Return list with queue option lines.
    """

    # submit to queue ?
    if rcf.get('submit.to') == 'queue' :

        # queue type:
        queue = rcf.get('queue')

        # different options and commands:
        if queue == 'loadleveler' :
            qopt = QueueOptions_LoadLeveler( bname, rcf, step )
        elif queue == 'bsub' :
            qopt = QueueOptions_BSub( bname, rcf, step )
        elif queue == 'qsub' :
            qopt = QueueOptions_QSub( bname, rcf, step )
        else :
            # not supported ...
            logging.error( 'unsupported queue : %s' % queue )
            raise Exception
        #endif

    else :

        # sumitting to shell (background or foreground):
        qopt = ShellOptions( bname, rcf, step )

    #endif

    # ok
    return qopt

#enddef


# ***


def SubmitJob( job_script, rcfile ) :

    """
    Submit jobscript.
    Where to submit to (foreground,background, queue) is read from rcfile settings.
    """

    # external:
    import sys
    import logging
    from pyshell.base.bin import rc

    # read settings:
    rcf = rc.RcFile( rcfile )

    # where to ?
    submit_to = rcf.get('submit.to')
    # info ...
    logging.info( 'submit %s to %s ...' % (job_script,submit_to) )

    # call specific submit routines:
    if submit_to == 'foreground' :

        # call run script, catch errors:
        try :
            Run_Job_In_Foreground( job_script )
        except :
            logging.error( sys.exc_info()[1] )
            logging.error( 'from Run_Job_In_Foreground for %s' % job_script )
            raise Exception
        #endtry

    elif submit_to == 'background' :

        # call run script, catch errors:
        try :
            Submit_Job_To_Background( job_script, rcf )
        except :
            logging.error( 'from Submit_Job_To_Background for %s' % job_script )
            raise Exception
        #endtry

    elif submit_to == 'queue' :

        # queue type:
        queue = rcf.get('queue')

        # different options and commands:
        if queue == 'loadleveler' :
            Submit_Job_To_LoadLeveler( job_script, rcf )
        elif queue == 'bsub' :
            Submit_Job_To_BSub( job_script, rcf )
        elif queue == 'qsub' :
            Submit_Job_To_QSub( job_script, rcf )
        else :
            # not supported ...
            logging.error( 'unsupported queue : %s' % queue )
            raise Exception
        #endif

    else :

        # not supported ...
        logging.error( 'unsupported run environment : %s' % submit_to )
        sys.exit(1)

    #endif

    # ok
    return

#endif


# ======================================================================
# ===
# === foreground
# ===
# ======================================================================


def Run_Job_In_Foreground( job_script ) :

    """
    Run job script in foreground.
    """

    # external:
    import sys
    import os
    import logging
    import subprocess

    # setup command line, e.g. './myscript.jb' :
    command = os.path.join(os.curdir,job_script)

    # execute:
    retcode = subprocess.call( command )
    if retcode != 0 :
        logging.error( sys.exc_info()[1] )
        logging.error( 'from subprocess call to : %s' % command )
        raise Exception
    #endif

    # ok
    return

#enddef


# ======================================================================
# ===
# === background
# ===
# ======================================================================


def  Submit_Job_To_Background( job_script, rcf ) :

    """
    Submit job to background.
    """

    # external:
    import sys
    import os
    import logging
    import subprocess

    # basename for scripts etc is name of rcfile minus extension:
    bname,ext = os.path.splitext(job_script)

    # output files:
    job_stdout = bname+'.out'
    job_stderr = bname+'.err'
    job_info   = bname+'.info'

    # setup command line, e.g. './myscript.jb' :
    command = os.path.join(os.curdir,job_script)

    # re-direct standard output:
    command  = command+' > %s' % job_stdout
    # write error messages to seperate file:
    #command  = command+' 2> %s' % job_stderr
    command  = command+' 2>&1'

    # run in background, return process id:
    logging.info( 'run shell command : "%s" ...' % command )
    p = subprocess.Popen( command, shell=True )

    # info ...
    infotext = []
    infotext.append( '\n' )
    infotext.append( 'Summary:\n' )
    infotext.append( '\n' )
    infotext.append( '  job script      : %s\n' % job_script )
    infotext.append( '  standard output : %s\n' % job_stdout )
    infotext.append( '  standard error  : %s\n' % job_stderr )
    infotext.append( '\n' )
    infotext.append( 'Process snapshot:\n')
    infotext.append( '\n')
    p2 = subprocess.Popen( '/bin/ps -f -p %i' % p.pid, shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    for line in p2.stdout.readlines() : infotext.append( line )
    infotext.append( '\n')
    infotext.append( 'To manage this process:\n' )
    infotext.append( '\n' )
    infotext.append( '  # show process snapshot:\n' )
    infotext.append( '  ps -f -p %i\n' % p.pid )
    infotext.append( '  \n' )
    infotext.append( '  # kill process:\n' )
    infotext.append( '  kill %i\n' % p.pid )
    infotext.append( '  \n' )
    infotext.append( '  # follow standard output:\n' )
    infotext.append( '  tail -f %s\n' % job_stdout )
    infotext.append( '\n' )

    # write to log:
    for line in infotext : logging.info( line.strip() )

    # write to file:
    f = open( job_info, 'w' )
    f.writelines(infotext)
    f.close()

    # ok
    return

#enddef


# ======================================================================
# ===
# === Shell options
# ===
# ======================================================================


def ShellOptions( bname, rcf, step ) :

    """
    Return list with shell settings (in python).
    """

    # read shell lines directly from rcfile;
    # seperated by '\n' texts:
    lines = rcf.get( 'shell.options.%s' % step, default='' ).split('\\n')
    # add including newline:
    qopt = []
    for line in lines : qopt.append( '%s\n' % line.strip() )
    qopt.append( '\n' )

    # ok
    return qopt

#enddef


# ======================================================================
# ===
# === LoadLeveler queue
# ===
# ======================================================================


def QueueOptions_LoadLeveler( bname, rcf, step ) :

    """
    Return list with queue options.
    """

    # external:
    import math

    # init result:
    qopt = []

    # which step ?
    if step == 'default' :

        # list with options:
        opts = rcf.get( 'queue.ll.options.%s' % step ).split()
        # default options:
        for opt in opts :
            # get value:
            val = rcf.get( 'queue.ll.option.%s.%s' % (step,opt) )
            # write:
            qopt.append( '#@ %-20s = %s\n' % (opt,val) )
        #endfor
        # layout ...
        qopt.append( '\n' )

    else :

        # list with options:
        opts = rcf.get( 'queue.ll.options.%s' % step ).split()
        # default options:
        for opt in opts :
            # get value:
            val = rcf.get( 'queue.ll.option.%s.%s' % (step,opt) )
            # to be set ?
            if val == '<auto>' :
                # differs per option ...
                if opt == 'output' :
                    val = '%s_%s.out' % (bname,step)
                elif opt == 'error' :
                    val = '%s_%s.err' % (bname,step)
                #endif
            #endif
            # no, empty, or normal value ?
            if val == '<none>' :
                # skip this keyword:
                continue
            elif val == '' :
                # just the keyword:
                qopt.append( '#@ %s\n' % opt )
            else :
                # keyword and value:
                qopt.append( '#@ %-20s = %s\n' % (opt,val) )
            #endif
        #endfor
        # layout ...
        qopt.append( '\n' )

    #endif

    # ok
    return qopt

#enddef


# ***


def Submit_Job_To_LoadLeveler( job_script, rcf ) :

    """
    Submit job to LoadLeveler queue.
    """

    # external:
    import sys
    import os
    import logging
    import subprocess

    # basename for scripts etc is name of rcfile minus extension:
    bname,ext = os.path.splitext(job_script)

    # output files:
    job_info   = bname+'.info'

    # options passed directly to submission command:
    qopts = rcf.get( 'queue.ll.submit.options' )
    # add options passed to submit script:
    qopts = qopts+' '+rcf.get('submit.options')

    # info ...
    logging.info( '    launch ...' )

    # setup command line:
    command = 'llsubmit '+qopts
    # last argument is script:
    command = command+' '+job_script

    # info ...
    logging.info( '      command: %s' % command )

    # init submission info file:
    infotext = []
    infotext.append(  '\n' )

    # call submit command, trap errors:
    try:
        # submit; redirect errors to standard output:
        p = subprocess.Popen( command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    except :
        logging.error( sys.exc_info()[1] )
        logging.error( 'from subprocess.Popen( %s )' % command.split() )
        raise Exception
    #endtry
    # extract:
    outlines = p.stdout.readlines()
    # add to help info message:
    infotext = infotext + outlines
    # extract job id from last line:
    #   llsubmit: The job "c1a0303.4290133" with 3 job steps has been submitted.
    firstwords = 'llsubmit: The job'
    lastline = outlines[-1]
    if lastline.startswith(firstwords) :
        job_id = lastline.lstrip(firstwords).split()[0].replace('"','')
    else :
        job_id = '<job-id>'
    #endif

    # add help text to submission info:
    infotext.append(  '\n' )
    infotext.append(  'To manage LoadLeveler jobs:\n' )
    infotext.append(  '\n' )
    infotext.append(  '  llq [-u ${USER}]         # list [your] current jobs\n' )
    infotext.append(  '  llq %s             # list this job\n' % job_id )
    infotext.append(  '  llcancel %s        # kill this job\n' % job_id )
    infotext.append(  '\n' )

    # write to log:
    for line in infotext : logging.info( line.rstrip() )

    # write to file:
    f = open( job_info, 'w' )
    f.writelines(infotext)
    f.close()

    # ok
    return

#enddef


# ======================================================================
# ===
# === BSUB queue
# ===
# ======================================================================


def QueueOptions_BSub( bname, rcf, step ) :

    """
    Return list with queue options.
    """

    # init result:
    qopt = []

    # specials ...
    if step == 'run' :
        # mpi job ?
        if rcf.get('par.mpi','bool') :
            # number of MPI tasks:
            ntask = rcf.get('par.ntask','int')
            # parallel job:
            bsub_n = ntask
        else:
            # serial job:
            bsub_n = 1
        #endif
    else :
        # serial step:
        bsub_n = 1
    #endif

    # list with options:
    opts = rcf.get( 'queue.bsub.options' ).split()
    # default options:
    for opt in opts :
        # get value:
        val = rcf.get( 'queue.bsub.option.%s' % opt )
        # to be set ?
        if val == '<auto>' :
            # differs per option ...
            if opt == 'o' :
                val = '%s_%s.out' % (bname,step)
            elif opt == 'e' :
                val = '%s_%s.err' % (bname,step)
            elif opt == 'n' :
                val = str(bsub_n)
            #endif
        #endif
        # fill option line:
        qopt.append( '#BSUB -%s %s\n' % (opt,val) )
    #endfor
    # layout ...
    qopt.append( '\n' )

    # ok
    return qopt

#enddef


# ***


def Submit_Job_To_BSub( job_script, rcf ) :

    """
    Submit job to BSUB queue.
    """

    # external:
    import os
    import logging
    import subprocess

    # basename for scripts etc is name of rcfile minus extension:
    bname,ext = os.path.splitext(job_script)

    # output files:
    job_info   = bname+'.info'

    # options passed directly to submission command:
    qopts = rcf.get( 'queue.bsub.submit.options' )
    # add options passed to submit script:
    qopts = qopts+' '+rcf.get('submit.options')

    # info ...
    logging.info( '    launch ...' )

    # setup command line:
    command = 'bsub '+qopts
    # last argument is script:
    command = command+' '+job_script

    # info ...
    logging.info( '      command: %s' % command )

    # prepare for OS errors (file does not exist etc.)
    try:
        # submit; redirect errors to standard output:
        p = subprocess.Popen( command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    except OSError, err :
        logging.error( 'OSError: '+err.strerror )
        logging.error( 'from call : %s' % command )
        raise Exception
    #endtry
    # extract:
    outlines = p.stdout.readlines()
    # display:
    for line in outlines : logging.info( '      %s' % line.rstrip() )
    # standard output is:
    #   <jobname> <jobnr>
    # extract job nr:
    job_nr = outlines[0].split()[1].strip('<>')

    # info ...
    infotext = []
    infotext.append(  '\n' )
    infotext.append(  'Summary:\n' )
    infotext.append(  '\n' )
    infotext.append(  '  current dir     : %s\n' % os.getcwd() )
    infotext.append(  '  job script      : %s\n' % job_script )
    infotext.append(  '\n' )
    infotext.append(  'To manage this job:\n' )
    infotext.append(  '  \n' )
    infotext.append(  '  # kill job:\n' )
    infotext.append(  '  bkill %s\n' % job_nr )
    infotext.append(  '  \n' )
    infotext.append(  'To show all your running jobs:\n' )
    infotext.append(  '  \n' )
    infotext.append(  '  bjobs\n' )
    infotext.append(  '  \n' )

    # write to log:
    for line in infotext : logging.info( line.rstrip() )

    # write to file:
    f = open( job_info, 'w' )
    f.writelines(infotext)
    f.close()

    # ok
    return

#enddef


# ======================================================================
# ===
# === QSUB queue
# ===
# ======================================================================


def QueueOptions_QSub( bname, rcf, step ) :

    """
    Return list with queue options.
    """

    # init result:
    qopt = []

    # list with default options:
    opts = rcf.get( 'queue.qsub.options.default' ).split()
    # loop over options:
    for opt in opts :
        # get value:
        val = rcf.get( 'queue.qsub.option.%s' % opt )
        # fill option line:
        qopt.append( '#PBS -%s %s\n' % (opt,val) )
    #endfor

    # list with step specific options:
    opts = rcf.get( 'queue.qsub.options.%s' % step ).split()
    # loop over options:
    for opt in opts :
        # get value:
        val = rcf.get( 'queue.qsub.option.%s.%s' % (opt,step) ) # such as queue.qsub.option.N.init
        # special '<none>' value ? then skip:
        if val == '<none>' : continue
        # to be filled online ?
        if val == '<auto>' :
            if opt == 'o' :
                val = '%s_%s.out' % (bname,step)
            elif opt == 'e' :
                val = '%s_%s.err' % (bname,step)
            else :
                logging.error( 'could not fill value for "<auto>" for qsub option "%s"' % opt )
                raise ValueError
            #endif
        #endif
        # fill option line:
        qopt.append( '#PBS -%s %s\n' % (opt,val) )
    #endfor

    # layout ...
    qopt.append( '\n' )

    # ok
    return qopt

#enddef


# ***


def Submit_Job_To_QSub( job_script, rcf ) :

    """
    Submit job to QSUB queue.
    """

    # external:
    import os
    import logging
    import subprocess

    # basename for scripts etc is name of rcfile minus extension:
    bname,ext = os.path.splitext(job_script)

    # output files:
    job_info   = bname+'.info'

    # options passed directly to submission command:
    qopts = rcf.get( 'queue.qsub.submit.options' )
    # add options passed to submit script:
    qopts = qopts+' '+rcf.get('submit.options')

    # info ...
    logging.info( '    launch ...' )

    # setup command line:
    command = 'qsub '+qopts
    # last argument is script:
    command = command+' '+job_script

    # info ...
    logging.info( '      command: %s' % command )

    # prepare for OS errors (file does not exist etc.)
    try:
        # submit; redirect errors to standard output:
        p = subprocess.Popen( command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    except OSError, err :
        logging.error( 'OSError: '+err.strerror )
        logging.error( 'from call : %s' % command )
        raise Exception
    #endtry
    # extract:
    outlines = p.stdout.readlines()
    # display:
    for line in outlines : logging.info( '      %s' % line.rstrip() )
    # standard output is:
    #    jobnr
    # extract job nr:
    job_nr = outlines[0].split()[0]

    # info ...
    infotext = []
    infotext.append(  '\n' )
    infotext.append(  'Summary:\n' )
    infotext.append(  '\n' )
    infotext.append(  '  current dir     : %s\n' % os.getcwd() )
    infotext.append(  '  job script      : %s\n' % job_script )
    infotext.append(  '\n' )
    infotext.append(  'To manage this job:\n' )
    infotext.append(  '  \n' )
    infotext.append(  '  # kill job:\n' )
    infotext.append(  '  qdel %s\n' % job_nr )
    infotext.append(  '  \n' )
    infotext.append(  'To show all your running jobs:\n' )
    infotext.append(  '  \n' )
    infotext.append(  '  qstat -n [-u ${USER}]\n' )
    infotext.append(  '  \n' )

    # write to log:
    for line in infotext : logging.info( line.rstrip() )

    # write to file:
    f = open( job_info, 'w' )
    f.writelines(infotext)
    f.close()

    # ok
    return

#enddef


# ======================================================================
# ===
# === end
# ===
# ======================================================================

