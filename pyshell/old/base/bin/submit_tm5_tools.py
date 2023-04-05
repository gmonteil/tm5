#!/usr/bin/env python
import sys
import os
import logging
import subprocess
from pyshell.base.bin import rc
from pyshell.base.bin import submit_tm5_setup_rcfile


def Command_Line(rcf, exe, args, in_debugger) :

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

    cmnd_exec = exe
    cmnd_args = args
    cmndline = '%s %s' % (cmnd_exec,cmnd_args)

    return cmndline


def WriteAndSubmitNewJob(rcfile, bindir) :

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

    if type(rcfile) == str :
        rcf = rc.RcFile( rcfile )
    else :
        rcf = rcfile

    rcfile_next = submit_tm5_setup_rcfile.WriteRcfile( rcf )

    if rcfile_next :
        # write job file(s) for this period and return the (first) name;
        # last command in a file should submit the next job if necessary:
        logging.info( '  write jobfile for %s ...' % rcfile_next )
        jobfile_next = WriteJob( rcfile_next, bindir )
        logging.info( '  submit next job : %s' % jobfile_next )
        SubmitJob( jobfile_next, rcfile_next )


def WriteJob( rcfile, bindir ) :

    """
    jobfile = WriteJob(rcfile)
    Write job file given the settings in rcfile.
    The name of the jobfile is based on the name of the rcfile.
    The last command in the job should submit the next job,
    and the script is therefore written in python.
    """

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


def SubmitCommand( jobfile, rcfile) :
    rcf = rc.RcFile(rcfile)
    return SubmitCommand_Foreground(jobfile, rcfile)


def QueueOptions( bname, rcf, step ) :
    # sumitting to shell (background or foreground):
    return ShellOptions( bname, rcf, step )


def SubmitJob(job_script, rcfile) :
    rcf = rc.RcFile(rcfile)
    command = os.path.join(os.curdir,job_script)
    retcode = subprocess.call(command)
    if retcode != 0 :
        logging.error( sys.exc_info()[1] )
        logging.error( 'from subprocess call to : %s' % command )
        raise Exception


def ShellOptions( bname, rcf, step ) :

    """
    Return list with shell settings (in python).
    """

    # read shell lines directly from rcfile;
    # seperated by '\n' texts:
    lines = rcf.get( 'shell.options.%s' % step, default='' ).split('\\n')
    # add including newline:
    qopt = []
    for line in lines :
        qopt.append( '%s\n' % line.strip() )
    qopt.append( '\n' )

    return qopt
