import logging
import os
import go
import pycasso_tools


def DefineOptions(parser) :
    
    """
    Usage : DefineOptions( parser )
    Define options supported by user scripts.
    Arugments:
       parser   : optparse object 
    """
    
    # define flag to clean source code:
    parser.add_option( "-c", "--clean", 
                         help="""Remove high level object and module files before 
                            compilation. Low level objects from for example the 
                            'file_hdf' sources are not removed to speedup 
                            re-compilation of the model.
                            To have all objects removed, use the '--new' option
                            to create a complete new build; in older scripts
                            an option '-C' or '--clean-all' was defined for this.""",
                         dest="clean", action="store_true", default=False )


def StoreOptions( settings, values ) :
        
    """
    Add the parsed flags to a dictionairy.
    The values have data fields with names defined by 'dest'
    in the previous calls to 'parser.add_option' .
    """
    
    # translate options into a dictionairy if they should
    # replace rcfile values:
    if values.clean : settings['build.make.clean'    ] = True


def Build_FlagGroups( rcf, basic=False ) :

    """
    Return list of compiler flag groups to be used.
    """
    
    # default :
    flaggroups = ['default','real8']
    
    # add mpi ?
    if rcf.get('par.mpi','bool') :
        flaggroups.append('mpi')
    
    # include not standard flags ?
    if not basic :

        # add openmp ?
        if rcf.get('par.openmp','bool') :
            flaggroups.append('openmp')

        # read other:
        line = rcf.get('build.configure.flags')
        # add if necessary:
        if len(line) > 0 :
            flaggroups = flaggroups + line.split()

    #endif
    
    ## add adjoint ?
    #if rcf.get('var4d','bool') : flaggroups.append('adj')

    # ok
    return flaggroups

#enddef


# ***


def Build_Define( rcf, macro_group, mdefs ) :

    """
    Edit a list with macro's to be defined given
    the name of group of macro's and the other setings
    in the rcfile.
    """
    
    # external:
    import logging
    
    # info ...
    logging.info( '  user script Build_Define for macro group %s ...' % macro_group )
    
    # select on macro group name:
    if macro_group == 'tm5' :

        # macro to include checks on zooming;
        # required for zoom runs, but should be undefined 
        # for runs without zooming to speedup the program:
        mac = 'with_zoom'
        # regions defined via a list in the rcfile ?
        if len(rcf.get('regions',default='')) > 0 :
            # number of regions:
            nregions = len(rcf.get('regions').split())
        else :
            # the 'old' method defines the grid via source files;
            # key 'source.nregions' from the rcfile is used to set the actual number:
            nregions_text = rcf.get('source.nregions')
            # check for special value:
            if nregions_text == 'nregions_max' :
                nregions = 999
            else :
                nregions = int(nregions_text)
            #endif
        #endif
        # zoom regions ? then 
        if nregions > 1 :
            # define macro 'with_zoom' if not present yet:
            if mac not in mdefs :
                # add:
                mdefs.append( mac )
                # info ...
                logging.info( '    defined %s ...' % mac )
            #endif
        else :
            # no zooming, remove macro 'with_zoom' if present:
            if mac in mdefs :
                # remove:
                mdefs.remove( mac )
                # info ...
                logging.info( '    undefined %s ...' % mac )
            #endif
        #endif

        # macro's to enable MPI specific code:
        macs = ['MPI','with_mpi']
        # loop:
        for mac in macs :
            # MPI enabled ?
            if rcf.get('par.mpi','bool') :
                # define macro 'with_zoom' if not present yet:
                if mac not in mdefs :
                    # add:
                    mdefs.append( mac )
                    # info ...
                    logging.info( '    defined %s ...' % mac )
                #endif
            else :
                # MPI not enabled, remove if necessary:
                if mac in mdefs :
                    # add:
                    mdefs.remove( mac )
                    # info ...
                    logging.info( '    undefined %s ...' % mac )
                #endif
            #endif
        #endfor

    #endif
    
    # ok
    return mdefs
    
#enddef


# ***


def Build_Configure( rcf ) :

    """
    Configure a source code.
    This script is called from the source directory.
    Arguments:
      rcf  : dictionairy with settings from rc file
    """
    
    # external
    import logging
    
    # info ...
    logging.info( '  user script Build_Configure ...' )
    
    # write grid definitions:
    Build_Configure_Grid( rcf )
    
    # write tracer order include file:
    Build_Configure_TracerOrder( rcf )
    
    # check on depricated stuff ...
    Build_Configure_Check( rcf )

    # ok
    return
    
#enddef


# *


def Build_Configure_Grid( rcf ) :

    """
    Write source files for grid definition.
    """
    
    # external:
    import logging

    # info ...
    logging.info( '    user script Build_Configure_Grid ...' )
    
    # regions are defined either by a list of region names ('regions') ,
    # or the older method using 'source.nregions' and 'dims_grid.h' :
    if len(rcf.get('regions',default='')) > 0 :
        # write all source files:
        Build_Configure_Grid__regions( rcf )
    elif rcf.has_key('source.nregions') :
        # write only number of regions to include file:
        Build_Configure_Grid__nregions( rcf )
    else :
        logging.error( 'could not extract number of regions from rcfile settings;' )
        logging.error( 'either "regions" or "source.nregions" should be defined' )
        raise Exception
    #endif

    # ok
    return

#enddef


# *


def Build_Configure_Grid__nregions( rcf ) :

    """
    Write 'dims_grid.h' from rcfile settings.
    """
    
    # external:
    import logging
    import pycasso_tools

    # info ...
    logging.info( '      user script Build_Configure_Grid__nregions ...' )
    
    # which file to create ?    
    srcfile = 'dims_grid.h'

    # info ...
    logging.info( '    create %s ...' % srcfile )
    
    # number of regions; could be the name 'nregions_max', thus read as string:
    source_nregions = rcf.get('source.nregions')

    # fill lines:
    lines = []
    lines.append( '!\n' )
    lines.append( '! Define actual number of regions.\n' )
    lines.append( '! Included by \'dims_grid.F90\'.\n' )
    lines.append( '!\n' )
    lines.append( 'integer, parameter  ::  nregions = %s\n' % source_nregions )

    # update file if necessary ...
    pycasso_tools.update_text_file( srcfile, lines )
    
    # ok
    return

#enddef


# *


def Build_Configure_Grid__regions( rcf ) :

    """
    Write 'dims_grid.F90' from rcfile settings.
    """
    
    # which file to create ?
    srcfile = 'dims_grid.F90'
    
    # different versions could be made, depending on the release (former cycle);
    # get the target release, use ininite number ('latest' release) as default:
    build_release = rcf.get( 'build.configure.release', 'float', default=999.9 )

    # model regions:
    regions = rcf.get('regions').split()

    # actual number:
    nregions = len(regions)

    # maximum length of grid names:
    len_region_name = max(map(len,regions))

    # other ...
    maxref = rcf.get('region.maxref')
    dx = rcf.get('region.dx')
    dy = rcf.get('region.dy')
    dz = rcf.get('region.dz')

    # start with empty file:
    lines = []
    lines.append( '!#################################################################\n' )
    lines.append( '!\n' )
    lines.append( '! Grids.\n' )
    lines.append( '!\n' )
    lines.append( '!### macro\'s #####################################################\n' )
    lines.append( '!\n' )
    lines.append( 'module dims_grid\n' )
    lines.append( '\n' )
    lines.append( '  implicit none\n' )
    lines.append( '  \n' )
    lines.append( '  ! --- in/out ------------------------------\n' )
    lines.append( '  \n' )
    lines.append( '  public\n' )
    lines.append( '  \n' )
    lines.append( '  \n' )
    lines.append( '  ! --- const -------------------------------\n' )
    lines.append( '  \n' )
    lines.append( '  \n' )
    lines.append( '  ! Basic model definition: resolution etc. including some routines\n' )
    lines.append( '  ! to fill the data structure.\n' )
    lines.append( '\n' )
    lines.append( '  ! basic (coarsest) resolution in degrees for x and y (dz default 1.0)\n' )
    lines.append( '\n' )
    lines.append( '  real, parameter     ::  dx = %s\n' % dx )
    lines.append( '  real, parameter     ::  dy = %s\n' % dy )
    lines.append( '  real, parameter     ::  dz = %s\n' % dz )
    lines.append( '\n' )
    lines.append( '\n' )
    lines.append( '  ! Maximum number of zoom regions, \n' )
    lines.append( '  ! including the basic (coarsest grid) region;\n' )
    lines.append( '  ! arrays are allocated for each of these regions:\n' )
    lines.append( '  integer, parameter  ::  nregions_max = %i\n' % nregions )
    lines.append( '  \n' )
    lines.append( '  ! Actual number of zoom regions,\n' )
    lines.append( '  ! during testing this could be set to 1 to quickly run the model.\n' )
    lines.append( '  integer, parameter :: nregions = %s\n' % nregions )
    lines.append( '\n' )
    lines.append( '  ! region_name is used to recognise the METEO files\n' )
    lines.append( '  ! region_name is also used in the HDF output file name\n' )
    lines.append( '  ! region 1 should always be the global domain\n' )
    lines.append( '\n' )
    lines.append( '  integer, parameter  ::  len_region_name = %i\n' % len_region_name )
    lines.append( '  character(len=len_region_name), parameter  ::  region_name(1:nregions) = &\n' )
    line = '       (/ '
    for i in range(len(regions)) :
        if i > 0 : line = line + ', '
        fmt = "'%%-%is'" % len_region_name
        line = line + ( fmt % regions[i] )
    #endfor
    lines.append( line+'/)\n' )
    lines.append( '\n' )
    lines.append( '  ! coordinates (in degrees) for each region:\n' )
    lines.append( '  ! xcyc = 1 if the region has cyclic x-boundary conditions\n' )
    lines.append( '  ! touch_np = 1 if region touches the north pole\n' )
    lines.append( '  ! touch_sp = 1 if region touches the south pole\n' )
    lines.append( '  ! xbeg : the westmost border of the region\n' )
    lines.append( '  ! xend : the eastmost border of the region\n' )
    lines.append( '  ! ybeg : the southmost border of the region\n' )
    lines.append( '  ! yend : the northmost border of the region\n' )
    lines.append( '\n' )
    fields = ['xcyc','touch_np','touch_sp','xbeg','xend','ybeg','yend','im','jm']
    for ifield in range(len(fields)) :
        field = fields[ifield]
        line = '  integer, parameter  ::  %-8s(nregions) = (/ ' % field
        for iregion in range(len(regions)) :
            region = regions[iregion]
            if iregion > 0 : line = line + ', '
            val = rcf.get( 'region.%s.%s' % (region,field) )
            line = line + ( '%4i' % int(val) )
        #endfor
        lines.append( line+' /)\n' )
    #endfor
    lines.append( '\n' )
    lines.append( '\n' )
    lines.append( '  ! maximum refinement factor (can be arbitrary in principle):\n' )
    lines.append( '\n' )
    lines.append( '  integer, parameter :: maxref = %s\n' % maxref )
    lines.append( '\n' )
    lines.append( '  ! refinement factors for each region (<= maxref)\n' )
    lines.append( '  ! tref may differ from xref/yref. In the current \n' )
    lines.append( '  ! implementation it should be 1,2,4,6,...\n' )
    lines.append( '\n' )
    fields = ['xref','yref','zref','tref']
    for ifield in range(len(fields)) :
        field = fields[ifield]
        line = '  integer, parameter  :: %s(0:nregions) = (/ 1' % field
        for i in range(nregions) :
            #if i > 0 : line = line + ', '
            line = line + ', '
            val = rcf.get( 'region.%s.%s' % (regions[i],field) )
            line = line + ( '%4i' % int(val) )
        #endfor
        lines.append( line+' /)\n' )
    #endfor
    lines.append( '\n' )
    lines.append( '  ! Define the parent of each region. \n' )
    lines.append( '  ! Global region 1 should have parent 0 (globe single cell);\n' )
    lines.append( '  ! global surface region should have parent 1 (global region).\n' )
    line = '  integer, parameter  ::  parent(nregions) = (/ '
    for i in range(nregions) :
        if i > 0 : line = line + ', '
        val = rcf.get( 'region.%s.parent' % regions[i] )
        if val == 'globe' :
            ireg = 0
        else :
            ireg = regions.index(val) + 1
        #endif
        line = line + ( '%i' % ireg )
    #endfor
    lines.append( line+' /)\n' )
    lines.append( '\n' )
    lines.append( 'end module dims_grid\n' )

    # update file if necessary ...
    pycasso_tools.update_text_file( srcfile, lines )
    

def Build_Configure_TracerOrder( rcf ) :

    """
    Write 'chem_param.inc' from rcfile settings.
    """
    
    # macro defined ?
    if 'with_tracerorder' in rcf.get('build.configure.macro.define').split() :
    
        # external:
        import logging
        import pycasso_tools
        import traceback
        
        # info ...
        logging.info( '    user script Build_Configure_TracerOrder ...' )

        # special pycasso version of this module, might not exist yet:
        try :
            import pycasso_GenerateTracerOrder
        except :
            logging.error( 'Could not import "pycasso_GenerateTracerOrder" ;' )
            logging.error( 'not converted to Pycasso scripting yet ?' )
            logging.error( 'Procedure in module should look like :' )
            logging.error( '' )
            logging.error( '    def GenerateTracerOrder( rcf ) : ' )
            logging.error( '        ' )
            logging.error( '        # extract values from rcfile:' )
            logging.error( '        procs = rcf.get( \'par.ntask\', \'int\' )' )
            logging.error( '        ...' )
            logging.error( '        ' )
            logging.error( '        # fill lines (do not forget the newline character!) :' )
            logging.error( '        lines = []' )
            logging.error( '        lines.append( \'  integer, parameter :: ntrace = %i\\n\' % cnt )' )
            logging.error( '        ...' )
            logging.error( '        ' )
            logging.error( '        #ok' )
            logging.error( '        return lines' )
            logging.error( '        ' )
            logging.error( '    #enddef' )
            logging.error( '' )
            logging.error( traceback.format_exc() )
            raise Exception
        #endtry

        # which file to create ?    
        srcfile = 'chem_param.inc'

        # info ...
        logging.info( '    create %s ...' % srcfile )

        # return lines with file content:
        try :
            lines = pycasso_GenerateTracerOrder.GenerateTracerOrder( rcf )
        except :
            logging.error( 'Something wrong with call to GenerateTracerOrder.' )
            logging.error( traceback.format_exc() )
            raise Exception
        #endtry

        # update file if necessary ...
        pycasso_tools.update_text_file( srcfile, lines )
        
    #endif
    
    # ok
    return

#enddef

    
# *


def Build_Configure_Check( rcf ) :

    """
    Check source file for undesired features.
    """
    
    # external:
    import logging
    import os
    import fnmatch

    # info ...
    logging.info( '    user script Build_Configure_Check ...' )
    
    # keywords for checks to be performed:
    checknames = rcf.get( 'build.configure.checks', default='' )
    # empty ? then leave:
    if len(checknames) == 0 : return
    
    # error or just warnings ?
    with_error = rcf.get( 'build.configure.checks.error', 'bool', default=False )
    
    # set flag:
    any_warning = False
    
    # list files:
    srcfiles = os.listdir( os.curdir )
    srcfiles.sort()
    
    # loop over checks:
    for checkname in checknames.split() :
    
        # paterns:
        test_msg   = rcf.get( 'build.configure.check.%s.msg'   % checkname )
        test_files = rcf.get( 'build.configure.check.%s.files' % checkname ).split()
        test_skip  = rcf.get( 'build.configure.check.%s.skip'  % checkname ).split()
        test_line  = rcf.get( 'build.configure.check.%s.test'  % checkname )
        test_help  = rcf.get( 'build.configure.check.%s.help'  % checkname )

        # set flags:
        matching_files = False

        # loop over files:
        for srcfile in srcfiles :

            # match with patern ?
            match = False
            for pat in test_files :
                match = match or fnmatch.fnmatch(srcfile,pat)
                if match : break
            #endfor
            if not match : continue 

            # ... except if the name matches other patterns:
            match = False
            for pat in test_skip :
                match = match or fnmatch.fnmatch(srcfile,pat)
                if match : break
            #endfor
            if match : continue

            # read file:
            f = open( srcfile )
            lines = f.readlines()
            f.close()

            # search for something that is there, or something that is not there ...
            if test_line.startswith('not') :
                # by default no match:
                match = False
                # loop over lines:
                for line in lines :
                    # test on this line:
                    match = match or (not eval( test_line ))
                    # try next file after first match ...
                    if match : break
                #endfor
                # check next file if the requested code was found:
                if match : continue
                # revert:
                match = not match
            else :
                # by default no match:
                match = False
                # loop over lines:
                for line in lines[0:20] :
                    # test on this line:
                    match = match or eval( test_line )
                    # leave after first match ...
                    if match : break
                #endfor
            #endfor

            # found something ?
            if match :
                # info ...
                if not matching_files : logging.warning( '      %s : [found]' % test_msg )
                logging.warning( '        %s' % srcfile )
                # reset flags:
                matching_files = True
                any_warning = True
            #endif

        #endfor   # source files

        # info ...
        if matching_files :
            # display error message ?
            if with_error :
                # display help text; split at '\n' for newlines:
                for helpline in test_help.split('\\n') : logging.warning(helpline)
            #endif
        else :
            # no warnings for this test ...
            logging.info( '      %s [none ]' % test_msg )
        #endif

    #endfor  # checks
    
    # check for unknown macro's ?
    checkname = 'unknown_macro'
    flag = rcf.get( 'build.configure.check.%s' % checkname )
    if flag :

        # settings:    
        test_msg   = rcf.get( 'build.configure.check.%s.msg' % checkname )

        # names of macro groups:
        macgroups = rcf.get( 'build.configure.macro.groups' ).split()
        # collect all supported macro's:
        macall = []
        for macgroup in macgroups :
            macs = rcf.get( 'build.configure.macro.%s.all' % macgroup ).split()
            macall = macall + macs
        #endfor

        # flag ...
        logged_msg = False

        # loop over files:
        for srcfile in srcfiles :

            # read file (only if not *.mod, *.o or directory):
            if os.path.isdir(srcfile) : continue
            if fnmatch.fnmatch(srcfile,"*.mod") : continue
            if fnmatch.fnmatch(srcfile,"*.o") : continue
            f = open( srcfile )
            lines = f.readlines()
            f.close()
            
            # flags:
            logged_srcfile = False
            
            # loop over lines:
            for iline in range(len(lines)) :
                # current:
                line = lines[iline].strip()
                # macro test ?
                if line.startswith('#ifdef') or line.startswith('#ifndef') :
                    # second element of line is macro name:
                    mac = line.split()[1].strip()
                    # not supported ?
                    if mac not in macall :
                        # test description if not done yet:
                        if not logged_msg :
                            logging.info( '      %s' % test_msg )
                            logged_msg = True
                        #endif
                        # intro if necessary:
                        if not logged_srcfile :
                            logging.error( '        unsupported macro(s) in %s :' % srcfile )
                            logged_srcfile = True
                        #endif
                        # line number and content:
                        logging.error( '        %6i : %s' % (iline,line) )
                        # set flag:
                        any_warning = True
                    #endif   # unsuported macro
                #endif  # line with macro test  
            #endfor   # lines

        #endfor   # source files

        # jippy ...
        if not logged_msg :
            # no warnings for this test ...
            logging.info( '      %s [none ]' % test_msg )
        #endif

    #endif   # test on unsupported macro's

    # break ?
    if any_warning and with_error :
        logging.error( 'some source code checks failed; break' )
        logging.error( '(set "build.configure.checks.error : False" in the expert.rc to avoid this error)' )
        raise Exception
    #endif
    
    # ok
    return

#enddef

    
# ***


def Build_Compiler( rcf ) :

    """
    Set compiler and linker names.
    Usually it is enough to read the name from the rcfile,
    but some compiler families have aliases for compilation with
    MPI or OpenMP enabled.
    Arguments:
      rcf  : dictionairy with settings from rc file
    Return values:
      fc,linker
    """
    
    # external
    import logging
    
    # info ...
    logging.info( '  user script Build_Compiler ...' )
    
    # extract compiler name:
    fc = rcf.get('compiler.fc')
    # or supporting openmp ?
    if rcf.get('par.openmp','bool') : fc = rcf.get( 'compiler.fc.openmp' )

    # or with mpi support ?
    if rcf.get('par.mpi','bool') : 
        # extract compiler name:
        fc = rcf.get( 'mpi.compiler.fc' )
        # or supporting openmp ?
        if rcf.get('par.openmp','bool') : fc = rcf.get( 'mpi.compiler.fc.openmp' )
    #endif

    # f77 compiler, by default the same as fc:
    f77 = rcf.get( 'compiler.f77', default=fc )
    
    # assume linker is the same:
    linker = fc

    # info ...
    logging.debug( '    fortran compiler : %s' % fc )
    logging.debug( '    linker           : %s' % linker )
    
    # ok
    return fc,f77,linker
    
#enddef


# ***


def Build_Make( rcf ) :

    """
    Make and install an executable.
    This script is called from the source directory.
    Arguments:
      rcf  : dictionairy with settings from rc file
    """
    
    clean     = rcf.get('build.make.clean','bool',default=False)
    if clean :
        logging.debug( '  make clean ...' )
        for f in os.listdir(os.curdir) :
            if f.endswith('.o') or f.endswith('.mod') :
                # skip the most basic toolboxes ..
                if f.startswith('parray'   ) : continue
                if f.startswith('file_hdf' ) : continue
                if f.startswith('mdf'      ) : continue
                if f.startswith('file_grib') : continue
                if f.startswith('go'       ) : continue
                if f.startswith('binas'    ) : continue
                if f.startswith('num'      ) : continue
                if f.startswith('phys'     ) : continue
                if f.startswith('grid'     ) : continue
                if f.startswith('tmm'      ) : continue
                os.remove(f)

    # module dir ?
    mdir = rcf.get('compiler.mdir',default='None')
    if mdir != 'None' :
        if not os.path.exists( mdir ) :
            os.makedirs( mdir )

    # number of jobs available for make:
    build_jobs = rcf.get('build.jobs', default='')

    # get maker command; replace some keys:
    maker = rcf.get('maker').replace('%{build.jobs}',build_jobs)
    exe   = rcf.get('build.make.exec')
    command = maker.split() + ['-f', 'Makefile', exe]
    logging.info(str(command))
    go.subprocess.watch_call(command)

