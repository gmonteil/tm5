#! /usr/bin/env python
import sys
import os
import subprocess
import shutil
import re
import logging
import datetime

from pyshell.base.bin import rc
from pyshell.base.bin import pycasso_tools
from pyshell.base.bin import pycasso_user_scripts_tm5 as pycus


def main2(rcfile):
    rcf = rc.RcFile(rcfile)

    key = 'logfile'
    if rcf.has_key(key):
        logfile = rcf.get(key)

    if rcf.get('build.copy', 'bool'):
        Build_Copy(rcf, pycus)

    if rcf.get('build.configure', 'bool'):
        Build_Configure(rcf, pycus)

    if rcf.get('build.make', 'bool'):
        Build_Make(rcf, pycus)

    rundir = rcf.get('rundir')
    if len(rundir) > 0:
        pycasso_tools.CreateDirs(rundir)
        os.chdir(rundir)

    ifiles = rcf.get('install.copy')
    if len(ifiles) > 0:
        for ifile in ifiles.split():
            # if the file contains a colon ':' ...
            if ':' in ifile:
                # ... the name after the colon is the target name
                sour, targ = ifile.split(':')
            else:
                # ... otherwise, copy to current directory:
                sour, targ = ifile, os.path.curdir
            if not os.path.exists(sour):
                logging.error('source file not found : %s' % sour)
                raise IOError
            shutil.copy(sour, targ)

    newrc = rcf.get('install.rc')
    if len(newrc) > 0:
        rcf.WriteFile(newrc)

    submit_command = rcf.get('submit.command')

    build_prefix = rcf.get('build.prefix')
    if rcf.get('submit.relpaths', 'bool'):
        build_prefix = os.path.relpath(build_prefix, start=rundir)

    flag = rcf.get('submit.auto', 'bool')
    if flag:
        # change to run directory:
        os.chdir(rundir)
        retcode = subprocess.call(submit_command.split())
        if retcode != 0:
            logging.error('from submission of : %s' % submit_command)
            raise Exception


def Main( args, pycasso_user_scripts) :

    options, rcfile = ParseArguments( args, pycus )
    
    # initialise logging system:
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)

    # initialise messages to standard output:
    stdout_handler = Start_StdOut_Log(logger)
    
    logging.info( 'options and arguments ...' )
    logging.info( '  parsed options         : %s' % options )
    logging.info( '  parsed rcfile argument : %s' % rcfile )
    logging.info( 'read settings from %s ...' % rcfile )
    
    rcf = rc.RcFile(rcfile)
    
    if len(options) > 0 :
        logging.info( 'change rcfile settings given options ...' )
        for key,val in options.iteritems() :
            if rcf.has_key(key) :
                logging.info( '  reset "%s" to "%s" ...' % (key,str(val)) )
                rcf.replace( key, str(val) )
            else :
                logging.info( '  set "%s" to "%s" ...' % (key,str(val)) )
                rcf.add( key, str(val), comment='new key added after options passed to setup script' )
    
    key = 'logfile'
    if rcf.has_key(key) :
        logfile = rcf.get(key)
        logfile_handler = Start_File_Log( logger, logfile )
    else :
        logfile_handler = None

    # create a source ?
    if rcf.get('build.copy','bool'):
        Build_Copy( rcf, pycus )

    if rcf.get('build.configure','bool'):
        Build_Configure( rcf, pycus )

    if rcf.get('build.make','bool'):
        Build_Make( rcf, pycus )

    # change to destination directory ?  
    rundir = rcf.get('rundir')
    if len(rundir) > 0 :
        pycasso_tools.CreateDirs(rundir)
        os.chdir(rundir)

    # copy files ?
    ifiles = rcf.get('install.copy')
    if len(ifiles) > 0 :
        for ifile in ifiles.split() :
            # if the file contains a colon ':' ...
            if ':' in ifile :
                # ... the name after the colon is the target name
                sour, targ = ifile.split(':')
            else :
                # ... otherwise, copy to current directory:
                sour,targ = ifile,os.path.curdir
            if not os.path.exists(sour) :
                logging.error( 'source file not found : %s' % sour )
                raise IOError
            shutil.copy(sour, targ)

    newrc = rcf.get('install.rc')
    if len(newrc) > 0 :
        rcf.WriteFile(newrc)

    submit_command = rcf.get( 'submit.command' )

    build_prefix = rcf.get('build.prefix')
    if rcf.get('submit.relpaths','bool') :
        build_prefix = os.path.relpath(build_prefix,start=rundir)

    if rcf.get('submit.auto', 'bool'):
        os.chdir( rundir )
        logging.info( '  %s' % submit_command )
        retcode = subprocess.call( submit_command.split() )
        if retcode != 0 :
            logging.error( 'from submission of : %s' % submit_command )
            raise Exception

    if logfile_handler != None :
        logfile_handler.close()
    stdout_handler.close()
    

def ParseArguments( args, pycus ) :

    # external:
    import optparse
    
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
    

def Start_StdOut_Log(logger) :

    logformat = '[%(levelname)-8s] %(message)s'
    formatter = logging.Formatter(logformat)
    
    # handler for standard output:
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)

    return stdout_handler
    

def Start_File_Log(logger,logfile) :

    # set format of lines written to logging:
    if sys.version_info < (2, 5):
        logformat   = "%(asctime)s %(name)-10s, line %(lineno)4i  [%(levelname)-10s] %(message)s"
    else:
        logformat   = '%(lineno)4i %(filename)-30s -> %(funcName)-30s [%(levelname)-8s] %(message)s'
    formatter = logging.Formatter(logformat)
    
    logfile_handler = logging.FileHandler(filename=logfile,mode='w')
    logfile_handler.setLevel(logging.DEBUG)
    logfile_handler.setFormatter(formatter)
    logger.addHandler(logfile_handler)
    
    return logfile_handler


def Build_Copy( rcf, pycus ) :

    remove_existing_build = rcf.get('build.new','bool')
    prefix = rcf.get('build.prefix')

    if rcf.get('build.prefix.extend','bool') :
        prefix_ext = prefix
        flags = pycus.Build_FlagGroups(rcf)
        if len(flags) > 0 :
            for flag in flags : prefix_ext = prefix_ext+'_'+flag
        else :
            prefix_ext = prefix_ext+'_'
        logging.info('  build prefix extended : %s ' % prefix_ext)
        pycasso_tools.CreateDirs(prefix_ext, forceclean=remove_existing_build)
        if os.path.lexists(prefix) :
            os.remove(prefix)
        os.symlink(os.path.basename(prefix_ext), prefix)
    else :
        pycasso_tools.CreateDirs(prefix, forceclean=remove_existing_build)
    
    subdirs = rcf.get('build.copy.subdirs').split()
    for subdir in subdirs :
        sdir = os.path.join(prefix,subdir)
        pycasso_tools.CreateDirs(sdir)

    # loop over source directories:
    sourcedirs = rcf.get('build.copy.dirs')
    if len(sourcedirs) == 0 :
        logging.info( '  no source directories specified ...' )
        return

    # info ...
    logging.info( '  copy files from source directories...' )
    # some flags ...
    flag_remove__part = rcf.get('build.copy.remove.__part', 'bool')

    # Create a dict with files to copy
    # dict key is the filename, dict value is the path
    # if a same file is found in several projects, the value of the last scanned project is kept
    #
    # we copy only the files needed, and don't preserve their access time (so make doesn't consider them as up-to-date)
    files_to_copy = {}

    for sourcedir in sourcedirs.split():
        # info ...
        logging.info('    scanning %s ...' % sourcedir)
        # should be a directory ...
        if not os.path.isdir(sourcedir) :
            logging.error('specified source dir is not an existing directory : %s' % sourcedir)
            raise RuntimeError

        # empty ? then add something for current directory:
        if len(subdirs) == 0 :
            subdirs = os.path.curdir

        # loop over sub directories:
        for subdir in subdirs :
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

                            if stat in ['new', 'differ']:
                                # if file is new or different from what's already in the build dir, add it to the list
                                files_to_copy[os.path.join(subdir, outfile)] = sourcefile
                            elif os.path.join(subdir, outfile) in files_to_copy:
                                # if the file is identical to what is already in the build dir AND the file has been marked
                                # as file to copy from a previously scanned directory, then delete that previous marking
                                # e.g. if file x.f90 has a version in base and a version in proj A:
                                # - the version from proj A will be in the build dir (from a previous build with proj A),
                                # - when scanning the base, its status will be "differ", therefore it will be marked as file to copy
                                # - when scanning proj A, its status will be "...", therefore it won't be marked as file to copy
                                # - if we leave it like that, base/x.f90 will be copied to build, instead of proj/A/x.f90
                                # - if we delete it from the list of files to copy at this stage, we get the correct version (from proj A) in the build dir.
                                # - if a proj B comes after and has a third version of this file, it will be marked as file to copy normally.
                                del files_to_copy[os.path.join(subdir, outfile)]

    for k, v in files_to_copy.items():
        print("Copy %s to %s"%(v, os.path.join(prefix, k)))
        shutil.copy(v, os.path.join(prefix, k))


def Build_Configure( rcf, pycus ) :

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
    for flaggroup in flaggroups :
        fflags  =  fflags.strip()+' '+rcf.get('compiler.flags.'+flaggroup+'.fflags' )
        ldflags = ldflags.strip()+' '+rcf.get('compiler.flags.'+flaggroup+'.ldflags')

    flaggroups_basic = pycus.Build_FlagGroups( rcf, basic=True )
    for flaggroup in flaggroups_basic :
        fflags_basic = fflags_basic.strip()+' '+rcf.get('compiler.flags.'+flaggroup+'.fflags' )

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
        keybase = 'build.configure.macro.'+group
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
            for mdef in macros_def :
                if '=' in mdef :
                    mname,mval = mdef.split('=')
                    if mname in macros_all :
                        src.append( '#define %s %s\n' % (mname,mval) )
                elif mdef in macros_all :
                    src.append( '#define %s\n' % mdef )

            # write the source, replace existing file only if it was different:
            pycasso_tools.update_text_file( macros_hfile, src )

        for m in macros_all :
            if m not in macros_supported : macros_supported = macros_supported + macros_all
    
    # check for macro's that are not supported yet:
    any_error = False
    for macr in macros_defined :
        if macr not in macros_supported :
            if not any_error :
                # initial error message:
                logging.error( "one or more macro's have been defined that are not listed" )
                logging.error( "in any 'build.configure.macro.*.all' lists:" )
            logging.error( "  %s" % macr )
            any_error = True
    
    if any_error :
        raise Exception

    # create list of macro definitions as command line arguments, e.g. -Dwith_this_flag etc:
    fc_defs = ''   # for fortran compiler
    mk_defs = ''   # for makedepf90

    # add macro definitions ?
    define_D = rcf.get( 'build.configure.define.D', 'bool' )
    if define_D :
        fc_D = rcf.get('compiler.defineflag',default='-D')
        mk_D = '-D'
        for mdef in macros_defined :
            fc_defs = fc_defs.strip()+(' %s%s' % (fc_D,mdef) )
            mk_defs = mk_defs.strip()+(' %s%s' % (mk_D,mdef) )
    fflags = fflags.strip()+' '+fc_defs
    fflags_basic = fflags_basic.strip()+' '+fc_defs
    
    #
    # * libraries
    #
    
    # get list of default library names to be linked:
    libnames = rcf.get('build.configure.libs')
    
    # convert to list:
    libnames = libnames.split()

    # loop over defined macro's:
    for mdef in macros_defined :
        # read list of libraries that should be added if this macro is defined:
        libnames_ifdef = rcf.get( 'build.configure.libs.ifdef.%s' % mdef, default='' )
        if len(libnames_ifdef) > 0 :
            libnames = libnames+libnames_ifdef.split()
            logging.debug( '    %s    (macro `%s` defined)' % (libnames_ifdef,mdef) )

    # get list of all supported library names:
    libnames_all = rcf.get( 'build.configure.libs.all' ).split()
    
    # check if all libraries specified to be used are actually supported ...
    for libname in libnames :
        if libname not in libnames_all :
            logging.error( 'library name `%s` not in `build.configure.libs.all` list ...' % libname )
            raise Exception
    
    logging.debug( '  libraries linked (in this order):' )
    # now add compiler and linker flags ;
    # loop over all supported libraries (this specfifies the linking order!)
    for libname in libnames_all :
        if libname in libnames :
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

    logging.debug( '  change back to %s ...' % current_dir )
    os.chdir( current_dir )
    

def Build_Make( rcf, pycus ) :
    logging.info('make executable ...')
    Make_dir = rcf.get('build.make.dir')
    os.chdir(Make_dir)
    pycus.Build_Make(rcf)