#! /usr/bin/env python

"""
*            PYCASSO - PYthon Compile And Setup Scripts Organizer            *

INSTALL FOR YOUR APPLICATION

    Assume you have an application named 'yourmodel'.
    To create a compile/setup script for this application,
    copy the template setup script to a name that you and other users
    will recoqnize as the main setup script:

       cp  pycasso/py/pycasso_setup_template  setup_yourmodel

    Edit the 'Settings' section in this new script to set the location 
    of the PYCASSO scripts. Also specify the name of the module file
    with user scripts. If these are in a subdirectory './base/py', use fore example:

       # path to the PYCASSO modules:
       base_py = os.path.join('.','pycasso','py')

       # name of file with PYCASSO user scripts for this setup:
       pycasso_user_scripts = 'pycasso_user_scripts_tm5'

    Copy the template settings file 'pycasso__template.rc'
    to a suitable name:

        cp  pycasso/rc/pycasso_template.rc  yourmodel.rc

    Edit the settings if necessary, and start with:

      ./setup_yourmodel yourmodel.rc


LOGGING

    Logging is implemented via the standard python 'logging' module.
    
    By default, one logger is defined that writes to the standard output;
    usually this is the terminal window, but in a batch environment this
    could be redirected to a file.
    In addition, a second logger could be setup if the 'logfile'
    option in the rc file is set.
    
    To add a new message to one of the logs, follow the examples in the code:

        logging.info( 'rc-file read successfully' )

    This would create a log-message of level INFO.
    Other options for log messages are:

        logging.debug    (msg)     # extra information to debug the scripts
        logging.warning  (msg)     # warnings about undesired behaviour
        logging.info     (msg)     # the standard messages
        logging.exception(msg)     # something wrong, not fatal
        logging.error    (msg)     # something wrong, nearly fatal
        logging.critical (msg)     # something wrong, probably fatal

    The threshold level for messages to be included is DEBUG for the file output.
    The threshold level for messages to be included is INFO for the screen, 
    unless 'options.verbose' is True (set by '-v' or '--verbose' on the 
    command line).

"""

# ------------------------------------------------
# begin
# ------------------------------------------------


def Main( args, pycasso_user_scripts ) :

    """
    Start the compilation and setup of the application.
        args
                List of unparsed arguments, probalby from 'sys.argv[1:]'.
        pycasso_user_scripts
                Name of file with user scripts for this particular setup.
    """
    
    # external:
    import sys
    import os
    import exceptions
    import subprocess
    import shutil
    import logging
    import traceback
    import datetime
    
    # tools:
    from pyshell.base.bin import rc
    from pyshell.base.bin import pycasso_tools
    from pyshell.base.bin import pycasso_user_scripts_tm5 as pycus
    options, rcfile = ParseArguments( args, pycus )
    
    # initialise logging system:
    logger = Start_Logger()
    
    # initialise messages to standard output:
    stdout_handler = Start_StdOut_Log(logger)
    
    # info ...
    logging.info( '' )
    tnow = datetime.datetime.now().isoformat(' ').split('.',1)[0]
    logging.info( 'Started script at %s ...' % tnow )
    logging.info( '' )
    
    # info ...
    logging.info( 'options and arguments ...' )
    logging.info( '  parsed options         : %s' % options )
    logging.info( '  parsed rcfile argument : %s' % rcfile )
    
    # info ...
    logging.info( 'read settings from %s ...' % rcfile )
    
    # read settings:
    rcf = rc.RcFile(rcfile)
    
    # options provided at command line ?
    if len(options) == 0 :
        logging.info( 'no rcfile settings to be added or changed by options' )
    else :
        logging.info( 'change rcfile settings given options ...' )
        for key,val in options.iteritems() :
            if rcf.has_key(key) :
                logging.info( '  reset "%s" to "%s" ...' % (key,str(val)) )
                rcf.replace( key, str(val) )
            else :
                logging.info( '  set "%s" to "%s" ...' % (key,str(val)) )
                rcf.add( key, str(val), comment='new key added after options passed to setup script' )
    
    logging.info( 'setup logging according to settings ...')
    
    flag = rcf.get( 'verbose', 'bool', default=False )
    if flag : 
        logging.info( '  verbose mode for standard output; print all messages ...' )
        stdout_handler.setLevel(logging.DEBUG)
    else :
        logging.info( '  quiet mode for standard output; print info messages only ...' )
    
    # setup logfile ?
    key = 'logfile'
    if rcf.has_key(key) :
        logfile = rcf.get(key)
        logging.info( '  open additional log file : %s' % logfile )
        logfile_handler = Start_File_Log( logger, logfile )
    else :
        logging.info( '  no log file; print to standard output only ...' )
        logfile_handler = None

    # test logging:
    logging.info    ('  test messages ...')
    logging.debug   ('    testing debug    message ...')
    logging.info    ('    testing info     message ...')
    logging.warning ('    testing warning  message ...')
    logging.error   ('    testing error    message ...')
    logging.critical('    testing critical message ...')
    
    # display settings:
    logging.debug( '' )
    logging.debug( 'Content of rc dictionary:' )
    for key in rcf.keys() :
        logging.debug( '  %s : %s' % (key,rcf.get(key)) )
    logging.debug( '[done]' )
    
    # create a source ?
    flag = rcf.get('build.copy','bool')
    if flag :
        logging.info( 'copy source to build directory ...' )
        Build_Copy( rcf, pycus )
    else :
        logging.info( 'no source to be copied ...' )

    # configure a source ?
    flag = rcf.get('build.configure','bool')
    if flag :
        logging.info( 'configure source ...' )
        Build_Configure( rcf, pycus )
    else :
        logging.info( 'no source to be configured ...' )

    # make an executable ?
    flag = rcf.get('build.make','bool')
    if flag :
        logging.info( 'make source ...' )
        Build_Make( rcf, pycus )
    else :
        logging.info( 'no source to be made ...' )

    # change to destination directory ?  
    rundir = rcf.get('rundir')
    if len(rundir) > 0 :
        # create target directory if necessary:
        pycasso_tools.CreateDirs( rundir )
        # info ...
        logging.info( 'change to run directory %s ...' % rundir )
        # goto this directory:
        os.chdir(rundir)
    else :
        # info ...
        logging.info( 'no run directory specified; stay here ...' )
    #endif

    # copy files ?
    ifiles = rcf.get('install.copy')
    if len(ifiles) > 0 :
        logging.info( 'install files ...' )
        for ifile in ifiles.split() :
            # if the file contains a colon ':' ...
            if ':' in ifile :
                # ... the name after the colon is the target name
                sour, targ = ifile.split(':')
            else :
                # ... otherwise, copy to current directory:
                sour,targ = ifile,os.path.curdir
            logging.info( '  copy "%s" to "%s" ...' % (sour,targ) )
            if not os.path.exists(sour) :
                logging.error( 'source file not found : %s' % sour )
                raise IOError
            shutil.copy( sour, targ )
    else :
        logging.info( 'no files to be installed ...' )
    
    # write rcfile ...
    newrc = rcf.get('install.rc')
    if len(newrc) > 0 :
        # info ...
        logging.info( '  install processed rcfile ...' )
        # write pre-processed rcfile:
        rcf.WriteFile( newrc )
    #endif

    #
    # * submit script / info / submitting
    #
    
    # settings for submit script:
    submit_script = rcf.get( 'submit.script' )
    submit_command = rcf.get( 'submit.command' )

    # where to search for scripts ?    
    build_prefix = rcf.get('build.prefix')
    # use paths relative to run directory ?
    if rcf.get('submit.relpaths','bool') :
        # only available in recent versions ...
        if sys.version_info[0]+0.1*sys.version_info[1] >= 2.6 :
            build_prefix = os.path.relpath(build_prefix,start=rundir)
        else :
            logging.warning( "Option `submit.relpaths` requires at least python version 2.6 ; use absolute path instead" )
        #endif
    #endif

    # full path to submit script:
    #submit_script_path = submit_script
    #submit_script_path = os.path.join( rundir, submit_script )
    # should exist ...
    #import pdb; pdb.set_trace()
    #if not os.path.exists(submit_script_path) :
    #    logging.error( 'submit script not found:' )
    #    logging.error( '  %s' % submit_script_path )
    #    logging.error( 'not included in "install.copy" list ?' )
    #    raise Exception
    #endif
    # insert path to submit modules:
    # pycasso_tools.modify_text_file( submit_script_path, 
    #                 "pypath_default = os.path.join( os.pardir, 'build', 'bin' )", 
    #                  "pypath_default = os.path.join('%s','bin')" % build_prefix )
    
    # info ...
    logging.info( '' )
    logging.info( 'To submit a job:' )
    logging.info( '' )
    indent = '  '
    # need to change to run directory ?
    if len(rundir) > 0 :
        logging.info( indent+'# change to the run directory:' )
        logging.info( indent+'cd %s' % rundir )
        logging.info( '' )
    #endif
    # disply usage text:
    logging.info( indent+'# submit a job (use --help for more information):' )
    logging.info( indent+'%s' % submit_command )
    # advertisement ...
    logging.info( '' )
    logging.info( 'For first glance on settings and results:' )
    logging.info( '' )
    logging.info( indent+'# run diadem postprocessor:' )
    logging.info( indent+'./tools/diadem/py/diadem %s' % rcfile )
    
    # submit automatically ?
    flag = rcf.get( 'submit.auto', 'bool' )
    if flag :
        # change to run directory:
        os.chdir( rundir )
        # info ...
        logging.info( '' )
        logging.info( 'submit automatically ...' )
        logging.info( '  %s' % submit_command )
        logging.info( '' )
        # call script:
        retcode = subprocess.call( submit_command.split() )
        if retcode != 0 :
            logging.error( 'from submission of : %s' % submit_command )
            raise Exception
        #endif
    #endif
    
    
    #
    # * end
    #
    
    # info ...
    logging.info( '' )
    tnow = datetime.datetime.now().isoformat(' ').split('.',1)[0]
    logging.info( 'End of script at %s .' % tnow )
    logging.info( '' )
    
    # close logs:
    if logfile_handler != None : logfile_handler.close()
    stdout_handler.close()
    
    # ok
    return
         
#enddef


# ***


def ParseArguments( args, pycus ) :

    """
    Define the arguments accepted by a pycasso run script.

    Usage:

        options,rcfile = ParseArguments( args )

    Arguments:

        args
                Passed from the calling script, probably equal to : sys.argv[1:]
        pycasso_user_scripts
                Name of file with user scripts for this particular setup.

    Return values:

        options       # object with following data fields:
          .verbose    # (boolean) gives extra logging messages to the screen

        rcfile    # name of settings file
      
    """

    # external:
    import logging
    import optparse
    
    # load user scripts:
    # pycus = __import__( pycasso_user_scripts )
    
    # set text for 'usage' help line:
    usage = "\n    %prog [options] rcfile\n    %prog -h|--help"
    
    # descriptive text
    description = "Driver script to compile and setup a model application. "\
        "The 'rcfile' is a textfile with settings read by the scripts or "\
        "the application, a template should be available with this script."
    
    # initialise the option parser:
    parser = optparse.OptionParser( usage=usage, description=description )
    
    # define verbose option:
    parser.add_option( "-v", "--verbose", 
                         help="""Print extra logging messages to standard output.
                            This option will set rcfile key 'verbose' to 'True'.""",
                         dest="verbose", action="store_true", default=False )
    # new build ?
    parser.add_option( "-n", "--new", 
                         help="""Create new build; remove old build directory.
                            This option will set rcfile key 'build.new' to 'True'.""",
                         dest="build_new", action="store_true", default=False )
    # how many jobs used for make etc ?
    parser.add_option( "-j", "--jobs", 
                         help="""Number of jobs (commands) to run simultaneously.
                            Empty value '' indicates unlimitted number.
                            Now only used to (re)set the number of jobs used by the maker ('gmake -j JOBS') .
                            This flag will replace the value of 'build.jobs' in the rcfile.""",
                         dest="jobs", action="store", default=None )
    # submit job ?
    parser.add_option( "-s", "--submit", 
                         help="""Submit the job after setup. 
                            See also the section on 'Submit options' below for 
                            options passed directly to the submit script.
                            This option will set rcfile key 'submit.auto' to 'True'.""",
                         dest="submit_auto", action="store_true", default=False )

    # options for submitting the job:
    group = optparse.OptionGroup( parser, "Submit options",
                         description="Options passed directly to the submit script, see its help text for details." )
    # where to submit to ?
    group.add_option( "-f", "--foreground", 
                         help="Run job in foreground.",
                         dest="submit_to", action="store_const", const='foreground' )
    group.add_option( "-b", "--background", 
                         help="Run job in background.",
                         dest="submit_to", action="store_const", const='background' )
    group.add_option( "-q", "--queue", 
                         help="Submit job to a queue system.",
                         dest="submit_to", action="store_const", const='queue' )
    # when submitted, run in debugger ?
    group.add_option( "-d", "--debugger",
                         help="Run executable in a debugger.",
                         dest="submit_debugger", action="store_true" )
    # add group:
    parser.add_option_group( group )
   
    # add the user options:
    group = optparse.OptionGroup( parser, "User model options",
                         description="These options are defined and handled in the 'pycasso_user_scripts_*' module." )
    pycus.DefineOptions( group )
    parser.add_option_group( group )

    # now parse the actual arguments:
    values,args = parser.parse_args( args=args )
    
    # at least rcfile should be specified as argument,
    # and no other arguments:
    if len(args) != 1 :
        parser.error("incorrect number of arguments\n")
    #endif
    
    # translate options into a dictionairy 
    # if they should replace rcfile values:
    opts = {}
    if values.verbose                 : opts['verbose'        ] = True
    if values.build_new               : opts['build.new'      ] = True
    if values.jobs != None            : opts['build.jobs'     ] = values.jobs
    if values.submit_auto             : opts['submit.auto'    ] = True
    if values.submit_to != None       : opts['submit.to'      ] = values.submit_to
    if values.submit_debugger != None : opts['submit.debugger'] = values.submit_debugger

    # add the parsed user options:
    pycus.StoreOptions( opts, values )
    
    # copy name of rcfile:
    rcfile = args[0]

    # return values:
    return opts,rcfile
    
#enddef


# ***


def Start_Logger() :

    """
    logger = Start_Logger()
    """

    # external:
    import logging

    # create logger
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)
    
    # ok
    return logger

#enddef


# ***


def Start_StdOut_Log(logger) :

    """
    stdout_handler = Start_StdOut_Log(logger)
    """

    # external:
    import sys
    import logging

    # set a format for screen use:
    logformat = '[%(levelname)-8s] %(message)s'
    # create formatter:
    formatter = logging.Formatter(logformat)
    
    # handler for standard output:
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)

    ## first messages:
    #logging.debug('testing debug message after start stdout log ...')
    #logging.info ('testing info  message after start stdout log ...')
    
    # ok
    return stdout_handler
    
#endif


# ***


def Start_File_Log(logger,logfile) :

    """
    logfile_handler = Start_File_Log(logger,logfile)
    
    Create handler for log file.
    Initial level is 'DEBUG', thus all messages will be written to the file.
    """
    
    # external:
    import sys
    import logging
    
    # set format of lines written to logging:
    if sys.version_info < (2, 5):
        logformat   = "%(asctime)s %(name)-10s, line %(lineno)4i  [%(levelname)-10s] %(message)s"
    else:
        logformat   = '%(lineno)4i %(filename)-30s -> %(funcName)-30s [%(levelname)-8s] %(message)s'
    #endif
    # create formatter:
    formatter = logging.Formatter(logformat)
    
    # now create a handler for the log file;
    # mode 'w' will cause the log file to be re-written:
    logfile_handler = logging.FileHandler(filename=logfile,mode='w')
    logfile_handler.setLevel(logging.DEBUG)
    logfile_handler.setFormatter(formatter)
    logger.addHandler(logfile_handler)
    
    ## first messages:
    #logging.debug('testing debug message after start file log ...')
    #logging.info ('testing info  message after start file log ...')
    
    # ok
    return logfile_handler

#enddef


# ***


def Build_Copy( rcf, pycus ) :

    import os
    import sys
    import shutil
    import logging
    from pyshell.base.bin import pycasso_tools
    import re
    
    remove_existing_build = rcf.get('build.new','bool')

    prefix = rcf.get('build.prefix')
    logging.info( '  build prefix : %s ' % prefix )
    
    if rcf.get('build.prefix.extend','bool') :
        prefix_ext = prefix
        flags = pycus.Build_FlagGroups( rcf )
        if len(flags) > 0 :
            for flag in flags : prefix_ext = prefix_ext+'_'+flag
        else :
            prefix_ext = prefix_ext+'_'
        logging.info( '  build prefix extended : %s ' % prefix_ext )
        pycasso_tools.CreateDirs( prefix_ext, forceclean=remove_existing_build )
        if os.path.lexists(prefix) :
            if os.path.islink(prefix) :
                os.remove( prefix )
            else :
                logging.error( 'could not replace "'+prefix+'" by a symbolic link; remove first' )
                raise Exception
        os.symlink( os.path.basename(prefix_ext), prefix )
    else :
        pycasso_tools.CreateDirs( prefix, forceclean=remove_existing_build )
    
    subdirs = rcf.get('build.copy.subdirs').split()
    for subdir in subdirs :
        sdir = os.path.join(prefix,subdir)
        pycasso_tools.CreateDirs(sdir)
    
    # loop over source directories:
    sourcedirs = rcf.get('build.copy.dirs')
    if len(sourcedirs) > 0 :
        # info ...
        logging.info( '  copy files from source directories...' )
        # some flags ...
        flag_remove__part = rcf.get('build.copy.remove.__part', 'bool')
        if flag_remove__part :
            logging.info( '    remove "__<name>" parts from sources files ...' )

        # loop over source directories:
        for sourcedir in sourcedirs.split():
            found_some_files = False
            # info ...
            logging.info( '    scanning %s ...' % sourcedir)
            # should be a directory ...
            if not os.path.isdir(sourcedir) :
                logging.error( 'specified source dir is not an existing directory : %s' % sourcedir )
                raise IOError
            
            # empty ? then add something for current directory:
            if len(subdirs) == 0 : subdirs = os.path.curdir
            # loop over sub directories:
            
            for subdir in  subdirs :
                sourcepath = os.path.join(sourcedir, subdir)
                for sfile in os.listdir(sourcepath) :
                    sourcefile = os.path.join(sourcepath, sfile)
                    if not os.path.isdir(sourcefile):
                        name, ext = os.path.splitext(sfile)
                        skipit = False
                        for pattern in rcf.get('build.copy.skip.ext').split() :
                            if re.search(pattern, ext) :
                                logging.debug( '        %-50s     %-30s [%-8s]' % (sourcefile,'','skip') )
                                skipit = True
                        if not skipit:
                            outfile = sfile
                            if flag_remove__part:
                                if '__' in sfile:
                                    name, ext = os.path.splitext(sfile)
                                    outfile = name.split('__')[0] + ext
                            if outfile not in rcf.get('build.copy.skip.file').split():
                                targetfile = os.path.join(prefix, subdir, outfile)
                                stat = 'new'
                                if os.path.exists(targetfile):
                                    if pycasso_tools.diff_text_files(sourcefile, targetfile):
                                        stat = 'differ'
                                    else:
                                        stat = '...'
                                logging.info('        %-50s  -> %-30s [%-8s]'%(sourcefile, os.path.join(subdir, outfile), stat))
                                # copy source to target, preserve times etc:
                                shutil.copy2(sourcefile, targetfile)
                    
                    found_some_files = True
                        

            if not found_some_files:
                logging.warning('  found no source files in standard subdirs %s of %s.  Mistake in source.dirs?' % (str(subdirs),sourcedir))
    else :
        logging.info( '  no source directories specified ...' )

    # # add a new formed directory to the python path?
    # flag = rcf.get('build.copy.pypath','bool')
    # if flag :
    #     logging.info( '  add subdir <prefix>/py to python path ...' )
    #     newdir = os.path.join(prefix,'py')
    #     sys.path.insert(0,newdir)
    # else :
    #     logging.info( '  no request for extension of the python path ...' )


def Build_Configure( rcf, pycus ) :

    # external:
    import os
    import sys
    import logging
    import subprocess
    
    # tools:
    from pyshell.base.bin import go
    from pyshell.base.bin import pycasso_tools
    
    # change directory:
    current_dir = os.getcwd()
    configure_dir = rcf.get('build.sourcedir')
    logging.debug( '  change to %s ...' % configure_dir )
    os.chdir( configure_dir )

    # call user script to set compilers and linker:
    fc, f77, linker = pycus.Build_Compiler(rcf)
    
    # start without any flags:
    fflags  = ''
    ldflags = ''
    libs    = ''
    
    # very basic flags to compile the too large mdf module ...
    fflags_basic = ''   
    
    # get list with names of compiler flag groups to be used:
    flaggroups = pycus.Build_FlagGroups( rcf )
    # loop over groups:
    for flaggroup in flaggroups :
        # add flags for this group:
        fflags  =  fflags.strip()+' '+rcf.get('compiler.flags.'+flaggroup+'.fflags' )
        ldflags = ldflags.strip()+' '+rcf.get('compiler.flags.'+flaggroup+'.ldflags')
    #endfor
    logging.info( '      fflags       : %s' % fflags )
    logging.info( '      ldflags      : %s' % ldflags )
    
    # idem for basic flags:
    flaggroups_basic = pycus.Build_FlagGroups( rcf, basic=True )
    # loop over groups:
    for flaggroup in flaggroups_basic :
        # add flags for this group:
        fflags_basic = fflags_basic.strip()+' '+rcf.get('compiler.flags.'+flaggroup+'.fflags' )
    #endfor
    logging.info( '      fflags basic : %s' % fflags_basic )
    
    # default defined macro's:    
    macros_def   = rcf.get('build.configure.macro.define' ).split()

    # macro groups:
    groups =  rcf.get('build.configure.macro.groups').split()

    # apply user changes to default list:
    for group in groups :
        # add (or remove!) macro defintions by user script:
        macros_def = pycus.Build_Define( rcf, group, macros_def )
    
    # create a 'clean' list of defined macro's without double definitions:
    macros_defined   = []
    for m in macros_def :
        if m not in macros_defined   : macros_defined   = macros_defined   + macros_def

    # initialize a list to store all supported macro's without double definitions:
    macros_supported = []

    # loop over groups:
    for group in groups :
        # start of the rc keys:
        keybase = 'build.configure.macro.'+group
        # values for this group:
        macros_all   = rcf.get(keybase+'.all'    ).split()
        macros_hfile = rcf.get(keybase+'.hfile'  )
        # write header file ?
        if len(macros_hfile) > 0 :
            # info ...
            logging.info( '  update %s (include file with macro definitions) ...' % macros_hfile )
            # fill text for include file in list of strings,
            # each line should end with the newline expression '\n' :
            src = []
            src.append( '!\n' )
            src.append( '! Include file with macro definitions.\n' )
            src.append( '!\n' )
            # loop over macro's to be defined:
            for mdef in macros_def :
                # vallue assigned ?
                if '=' in mdef :
                    # split in name and value:
                    mname,mval = mdef.split('=')
                    # in this group ?
                    if mname in macros_all :
                        # add line to file:
                        src.append( '#define %s %s\n' % (mname,mval) )
                else :
                    # in this group ?
                    if mdef in macros_all :
                        # add line to file:
                        src.append( '#define %s\n' % mdef )
            # write the source, replace existing file only if it was different:
            pycasso_tools.update_text_file( macros_hfile, src )
        # extend list:
        for m in macros_all :
            if m not in macros_supported : macros_supported = macros_supported + macros_all
    
    # check for macro's that are not supported yet:
    any_error = False
    for macr in macros_defined :
        # not supported ?
        if macr not in macros_supported :
            # any unsupported macro's found yet ?
            if not any_error :
                # initial error message:
                logging.error( "one or more macro's have been defined that are not listed" )
                logging.error( "in any 'build.configure.macro.*.all' lists:" )
            #endif
            # display problematic macro:
            logging.error( "  %s" % macr )
            # reset flag:
            any_error = True
    
    # any error ? then leave:
    if any_error : raise Exception

    # create list of macro definitions as command line arguments, e.g. -Dwith_this_flag etc:
    fc_defs = ''   # for fortran compiler
    mk_defs = ''   # for makedepf90
    # add macro definitions ?
    define_D = rcf.get( 'build.configure.define.D', 'bool' )
    if define_D :
        # compiler depended prefix for macro definition:
        fc_D = rcf.get('compiler.defineflag',default='-D')
        mk_D = '-D'
        # loop over macro's to be defined:
        for mdef in macros_defined :
            # add definition to command line argument list:
            fc_defs = fc_defs.strip()+(' %s%s' % (fc_D,mdef) )
            mk_defs = mk_defs.strip()+(' %s%s' % (mk_D,mdef) )
        #endfor
    #endif
    # add definitions to flags:
    fflags = fflags.strip()+' '+fc_defs
    # idem for basic flags:
    fflags_basic = fflags_basic.strip()+' '+fc_defs
    

    #
    # * libraries
    #
    
    # get list of default library names to be linked:
    libnames = rcf.get( 'build.configure.libs' )
    
    # info ...
    logging.debug( '  libraries to be used:' )
    if len(libnames) > 0 : logging.debug( '    %s    (default)' % libnames )
    
    # convert to list:
    libnames = libnames.split()

    # loop over defined macro's:
    for mdef in macros_defined :
        # read list of libraries that should be added if this macro is defined:
        libnames_ifdef = rcf.get( 'build.configure.libs.ifdef.%s' % mdef, default='' )
        # defined ?
        if len(libnames_ifdef) > 0 :
            # add:
            libnames = libnames+libnames_ifdef.split()
            # info ...
            logging.debug( '    %s    (macro `%s` defined)' % (libnames_ifdef,mdef) )
    
    # get list of all supported library names:
    libnames_all = rcf.get( 'build.configure.libs.all' ).split()
    
    # check if all libraries specified to be used are actually supported ...
    for libname in libnames :
        if libname not in libnames_all :
            logging.error( 'library name `%s` not in `build.configure.libs.all` list ...' % libname )
            raise Exception
    
    # info ...
    logging.debug( '  libraries linked (in this order):' )
    # now add compiler and linker flags ;
    # loop over all supported libraries (this specfifies the linking order!)
    for libname in libnames_all :
        # not in use ? then skip:
        if libname not in libnames : continue
        # info ...
        logging.debug( '    %s' % libname )
        # add include, module, and link flags:
        fflags = fflags.strip()+' '+rcf.get('compiler.lib.'+libname+'.fflags')
        libs   =   libs.strip()+' '+rcf.get('compiler.lib.'+libname+'.libs')
        # idem for basic flags:
        fflags_basic = fflags_basic.strip()+' '+rcf.get('compiler.lib.'+libname+'.fflags')
    

    #
    # * write compiler flags to makefile
    #
    
    # name of include file with dependencies to be written:
    makefile_flags = rcf.get('build.configure.flags.includefile')
    # info ...
    logging.info( '  write %s (compiler and flags) ...' % makefile_flags )
    logging.info( '    compiler fflags : '+fflags  )
    logging.info( '            ldflags : '+ldflags )
    logging.info( '               libs : '+libs    )
    # fill content; each line should end with the newline expression '\n' :
    src = []
    src.append( '#\n' )
    src.append( '# include file with compiler flags for Makefile.\n' )
    src.append( '#\n' )
    src.append( '\n' )
    src.append( '# compiler and linker:\n' )
    src.append( 'FC = %s\n' % fc )
    src.append( 'F77 = %s\n' % f77 )
    src.append( 'LINKER = %s\n' % linker )
    src.append( '\n' )
    src.append( '# compile flags:\n' )
    src.append( 'FFLAGS = %s\n' % fflags )
    src.append( '\n' )
    src.append( '# compile flags without optim:\n' )
    src.append( 'FFLAGS_BASIC = %s\n' % fflags_basic )
    src.append( '\n' )
    src.append( '# linker flags:\n' )
    src.append( 'LDFLAGS = %s\n' % ldflags )
    src.append( '\n' )
    src.append( '# library flags:\n' )
    src.append( 'LIBS = %s\n' % libs )
    src.append( '\n' )
    # write:
    pycasso_tools.write_text_file( makefile_flags, src )
    

    #
    # * cleanup
    #
    
    # remove files ?
    logging.info( '  remove files ...' )
    for sfile in rcf.get('build.configure.remove').split() :
        # info ...
        logging.debug( '    %s ...' % sfile )
        # remove file if present:
        if os.path.exists(sfile) : os.remove( sfile )

    # clean for defined macro's:
    logging.info( '  remove files for defined macro\'s ...' )
    for macr in macros_defined :
        # try to read list with files to be removed if macro is defined:
        sfiles = rcf.get( 'build.configure.remove.ifdef.%s' % macr, default='' )
        # line was found ?
        if len(sfiles) > 0 :
            # loop over files to be removed:
            for sfile in sfiles.split() :
                # info ...
                logging.debug( '    %s ...' % sfile )
                # remove file:
                if os.path.exists(sfile) : os.remove(sfile)

    # clean for undefined macro's:
    logging.info( '  remove files for undefined macro\'s ...' )
    for macr in macros_supported :
        # defined ? then next:
        if macr in macros_defined : continue
        # try to read list with files to be removed if macro is NOT defined:
        sfiles = rcf.get( 'build.configure.remove.ifndef.%s' % macr, default='' )
        # line was found ?
        if len(sfiles) > 0 :
            # loop over files to be removed:
            for sfile in sfiles.split() :
                # info ...
                logging.debug( '    %s ...' % sfile )
                # remove file:
                if os.path.exists(sfile) : os.remove(sfile)


    #
    # * user script
    #

    # call user script:
    pycus.Build_Configure( rcf )


    #
    # * deps
    #
    
    # create deps ?
    flag = rcf.get('build.configure.makedep','bool')
    if flag :
        
        # name of include file with dependencies to be written:
        makefile_deps = rcf.get('build.configure.makedep.includefile')
        # info ...
        logging.info( '  create %s ...' % makefile_deps )

        # name of executable to be build:
        exe = rcf.get( 'build.make.exec' )
        # file filter:
        files = rcf.get('build.configure.makedep.files')

        # command to generate a makefile:
        command = 'makedepf90 %s -o %s %s > %s' % (mk_defs, exe, files, makefile_deps)
        logging.info(command)
        os.system(command)
        # info ...
        # logging.info( '    run command: %s' % command )
        # import pdb; pdb.set_trace()
        # # run command:
        # try :
        #     # run as a shell command since the file list is probably '*.F90' :
        #     p = go.subprocess.call( command, shell=True )
        # except go.subprocess.CallingError, err :
        #     logging.error( err )
        #     raise Exception
        # except go.subprocess.StatusError, err :
        #     for line in err.stderr : logging.error(line)
        #     logging.error( err )
        #     raise Exception
        
        # # write result:
        # f = open( makefile_deps, 'w' )
        # for line in p.stdout : f.write(line+'\n')
        # f.close()
        # # add to log file:
        # logging.debug( '' )
        # logging.debug( '---[%s]------------------------------------------------------------' % makefile_deps )
        # for line in p.stdout : logging.debug(line)
        # logging.debug( '-------------------------------------------------------------------' )
        # logging.debug( '' )
        
    logging.debug( '  change back to %s ...' % current_dir )
    os.chdir( current_dir )
    

def Build_Make( rcf, pycus ) :

    # external:
    import os
    import logging
    
    # info ...
    logging.info( 'make executable ...' )
    
    # change directory:
    Make_dir = rcf.get('build.make.dir')
    logging.debug( '  change to %s ...' % Make_dir )
    os.chdir( Make_dir )

    # call user script:
    pycus.Build_Make( rcf )