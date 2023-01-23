
"""
Routine to setup the rcfile of a TM5 run that is actually read by the
executable.
"""

def WriteRcfile( rcf ):

    """
    Reads the settings from the input rcfile, and writes the rcfile that will
    be read by the exectuble.

    The name of new rcfile is returned.

    An import key in the input rcfile is the 'jobstep' number, this is used to
    decide whether a restart should be configured rather than an initial run.

    If no new run is to be started (for example because the end of the job
    chain is reached), the value None is returned.
    """

    # external:
    import os
    import sys
    import shutil
    import datetime
    import logging

    # info ...
    logging.info( 'write rcfile (to be read by executable) ...' )
    logging.info( '  read time period settings ...' )

    # time format in rcfiles:
    tfmt = '%Y-%m-%d %H:%M:%S'

    if sys.version_info[0:3] <= (2,4,3) :
        full_t1 = datetime.datetime( *map(int,rcf.get('timerange.start').replace('-',' ').replace(':',' ').split()) )
        full_t2 = datetime.datetime( *map(int,rcf.get('timerange.end'  ).replace('-',' ').replace(':',' ').split()) )
    else :
        full_t1 = datetime.datetime.strptime( rcf.get('timerange.start'), tfmt )
        full_t2 = datetime.datetime.strptime( rcf.get('timerange.end'  ), tfmt )
    logging.info( '    timerange start : %s' % str(full_t1) )
    logging.info( '    timerange end   : %s' % str(full_t2) )
    if full_t2 <= full_t1:
        logging.error( 'The start date (%s) is later than the end date (%s), please revise' % (full_t1.strftime(tfmt),full_t2.strftime(tfmt)) )
        raise ValueError

    # rcfile is valid for the previous step in the chain (or initial);
    # read the previous step number (or zero):
    prev_step = rcf.get( 'jobstep', 'int' )
    # check ...
    if prev_step < 0 :
        logging.error( 'step number in chain is %i, but should be 0 (initial) or 1,2,...' % step )
        raise Exception
    #endif
    # info ...
    logging.info( '  previous step   : %i' % prev_step )

    # split name of rcfile in base and extension:
    prev_bname,ext = os.path.splitext(rcf.filename)
    # new basename is initally the same:
    next_bname = prev_bname

    # copy settings into new rcfile:
    next_rcf = rcf

    # start of chain ?
    if prev_step == 0 :

        # fill first value as end of 'previous' step:
        prev_t2 = full_t1

    else :

        # continuation of chain ...

        # base name includes step indication, e.g. mytest_001 ;
        # remove the last 4 characters:
        next_bname = next_bname[0:-4]

        # 4D=var mode ?
        is_var4d = rcf.get('var4d','bool',default=False)

        # get runmode (1 = forward run):
        if is_var4d :
            # info ...
            logging.info( '    script options for 4D-var are enabled  ...' )
            logging.info( '    read runmode  ...' )
            # current runmode:
            runmode = rcf.get( '4DVAR.runmode', 'int' )
            # info ...
            logging.info( '      runmode = %i' % runmode )
        else :
            # forward run:
            runmode = 1
        #endif

        # setup for different run modes:
        #
        if runmode in [1,2,3,4,6,7,8] :    # no-iteration modes
        #

            # info ...
            logging.info( '  setup next forward run ...' )

            # create python datetime objects for previous time range:
            if sys.version_info[0:3] <= (2,4,3) :
                prev_t1 = datetime.datetime( *map(int,rcf.get('jobstep.timerange.start').replace('-',' ').replace(':',' ').split()) )
                prev_t2 = datetime.datetime( *map(int,rcf.get('jobstep.timerange.end'  ).replace('-',' ').replace(':',' ').split()) )
            else :
                prev_t1 = datetime.datetime.strptime( rcf.get('jobstep.timerange.start'), tfmt )
                prev_t2 = datetime.datetime.strptime( rcf.get('jobstep.timerange.end'  ), tfmt )
            #endif
            # info ...
            logging.info( '    previous timerange start : %s' % str(prev_t1) )
            logging.info( '    previous timerange end   : %s' % str(prev_t2) )
            # check ...
            if (prev_t2 < prev_t1) or (prev_t1 < full_t1) or (prev_t2 > full_t2) :
                logging.error( 'found strange previous range compared to full range' )
                raise Exception
            #endif
            # finished ?
            if prev_t2 == full_t2 :
                # no new file necessary ...
                next_rcf = None
                # info ...
                logging.info( '    end of full period reached; finish' )
            #endif

        #
        elif runmode == 5 : # 4D-var run
        #

            # info ...
            logging.info( '    setup next 4D-var run ...' )

            # file with 4d-var restart settings is rcfile with extension '.rs' instead of '.rc' :
            rsfile = prev_bname+'.rs'
            # check ...
            if not os.path.exists(rsfile) :
                logging.error( '4D-var restart file "%s" not found ...' % rsfile )
                raise Exception
            #endif
            # info ...
            logging.info( '      include settings from restart file in new rcfile:' )
            # read file:
            f = open(rsfile,'r')
            lines = f.readlines()
            f.close()
            # loop over lines; written in rcfile format
            for line in lines :
                # remove newlines etc:
                line = line.strip()
                # skip comment and empty lines:
                if line.startswith('!') : continue
                if len(line) == 0 : continue
                # value definitions, split at ':' :
                key,val = line.split(':')
                key = key.strip()
                val = val.strip()
                # this should be a key in the rcfile ...
                if not next_rcf.has_key(key) :
                    logging.error( 'restart file contains key which is not in rc file : %s' % key )
                    raise Exception
                #endif
                # replace in rcfile:
                next_rcf.replace( key, val )
                # info ...
                logging.info( '        %s   :  %s' % (key,val) )
            #endfor

            # finished ?
            val = next_rcf.get('outer_loop.finished','int')
            logging.info( 'check value of outer_loop.finished : %i' % val )
            if val == 0 :
                logging.info( '  --> VAR4D outer/inner loop to be performed or to be continued; submit next job' )
            elif val == 1 :
                logging.info( '  --> VAR4D final inner loop finished, final outer loop to be done; submit next job' )
            elif val == 2 :
                logging.info( '  --> VAR4D final outer loop finished; end of job chain' )
                next_rcf = None
            else :
                logging.error( '  --> unsupported value : %i' % val )
                raise Exception
            #endif

        #
        else :
        #
            # info ...
            logging.warning( "do not know how to setup next job for 4D-var runmode : %i" % runmode )
            logging.warning( "assume end of run" )
            # no new file necessary ...
            next_rcf = None

        #endif
        #

    #endif

    # new file to be written ?
    if next_rcf != None :

        # info ...
        logging.info( '  write new rcfile ...' )

        # increase counter:
        next_step = prev_step + 1
        # replace in rcfile:
        next_rcf.replace('jobstep',next_step)
        # info ...
        logging.info( '    next jobstep : %i' % next_step )

        # get time step
        td = next_rcf.get('jobstep.length', default='inf')
        # info ...
        logging.info( '    jobstep length : %s' % td )
        # single run or splitted ?
        if td in ['inf','infinite'] :
            # copy time range:
            next_t1 = full_t1
            next_t2 = full_t2
            # info ...
            logging.info( '    job will run complete period ...' )
        else :
            # start time of this period:
            next_t1 = prev_t2
            # set end time for this period:
            if td == 'month' :
                # break at end of this month
                if next_t1.month != 12:
                    next_t2 = datetime.datetime(next_t1.year  ,next_t1.month+1,01,00,00)
                else:
                    next_t2 = datetime.datetime(next_t1.year+1,              1,01,00,00)
            elif td == 'year':
                # break at the end of this year
                next_t2 = datetime.datetime(next_t1.year+1, 1, 1, 0, 0, 0)
            else :
                # assume that the interval specified is the number
                # of days to run forward before resubmitting
                next_t2 = datetime.datetime(next_t1.year,next_t1.month,next_t1.day) + datetime.timedelta(days=int(td))
            #endif
            # do not run beyond full range:
            if next_t2 > full_t2 : next_t2 = full_t2
            # info ...
            logging.info( '    job will run from %s to %s' % (str(next_t1),str(next_t2)) )
        #endif
        # replace time values:
        next_rcf.replace( 'jobstep.timerange.start', next_t1.strftime(tfmt) )
        next_rcf.replace( 'jobstep.timerange.end'  , next_t2.strftime(tfmt) )
        # name of previous output directory:
        prev_output_dir = rcf.get( 'output.dir' )

        #>>> not tested well, gave problems with N2O inversions;
        #    therefore commented until working version is made
        ## for 2nd and higer step the model might need to read data
        ## from the previous run ...
        #if next_step > 1 :
        #
        #    # replace previous output destination:
        #    next_rcf.replace( 'prev.output.dir', prev_output_dir )
        #
        #    # New istart : 33 if previous run saved restart files,
        #    # else use "save files".
        #    # Default restart.store is F as assumed in tm5_restart.F90
        #    restart_write = rcf.get('restart.write','bool',default='F')
        #    if restart_write :
        #        # restart from a 'restart' file:
        #        next_istart = '33'
        #    else :
        #        # restart from a 'save' file:
        #        next_istart = '3'
        #    #endif
        #
        #    # info ...
        #    logging.info( '    istart value : %s' % next_istart )
        #    # replace rcfile value:
        #    next_rcf.replace( 'istart', next_istart )
        #
        #    # restart files enabled ?
        #    if restart_write :
        #        # name of previous destination directory:
        #        prev_restart_dir = rcf.get( 'restart.write.dir' )
        #        # ensure that new restart file is read from the correct directory:
        #        next_rcf.replace( 'restart.read.dir', prev_restart_dir )
        #    #endif
        #
        ##endif  # next_step > 1
        #<<<

        # list of output directory extensions (could be empty):
        extensions = rcf.get( 'output.dir.extensions' ).split()
        # should be extended ?
        if len(extensions) > 0 :
            # get output directory without extensions:
            output_dir = next_rcf.get('output.dir.base')
            # loop over sub directories:
            for extension in extensions :
                # some special values:
                if extension == '<jobstep>' :
                    # 4-digit jobstep:
                    subdir = '%4.4i' % next_step
                elif extension == '<timerange>' :
                    # format for time: yyyymmdd_hhmnss
                    fmt = '%Y%m%d%H'
                    # start and end time in short format:
                    next_t1s = next_t1.strftime(fmt)
                    next_t2s = next_t2.strftime(fmt)
                    # fill name of subdirectory:
                    subdir = '%s-%s' % (next_t1s,next_t2s)
                elif next_rcf.has_key(extension) :
                    # read key from rcfile:
                    subdir = next_rcf.get(extension)
                else :
                    logging.error( 'unsupported output dir extension : %s' % extension )
                    raise Exception
                #endif
                # extend output path with subdirectory
                output_dir = os.path.join( output_dir, subdir )
            #endfor
            # replace original value in rcfile:
            next_rcf.replace('output.dir',output_dir)
        #endif

        # name of new rcfile including new step number:
        next_rcfile = '%s_%3.3i%s' % (next_bname,next_step,ext)
        # info ...
        logging.info( '    write to %s ...' % next_rcfile )
        # write :
        next_rcf.WriteFile( next_rcfile )

        # for backwards compatibility, write 'tm5_runtime.rc' too ?
        flag = next_rcf.get( 'jobstep.runtimerc.write', 'bool', default=False )
        if flag : Write_RuntimeRc( next_rcf )

    else :

        # dummy name ...
        next_rcfile = None

    #endif

    # ok
    return next_rcfile

#enddef

# *

def Write_RuntimeRc( rcf ) :

    """
    Write runtime rc file with time range.
    Necessary for backward compatibility with pre-pycasso scripts.
    """

    # external:
    import datetime

    # write rcfile to standard name expected by executable:
    rcf.WriteFile( 'tm5.rc' )

    # now write runtime rcfile with additional information ...

    # time format in rcfiles:
    tfmt = '%Y-%m-%d %H:%M:%S'

    # current values:
    step        = rcf.get( 'jobstep', 'int' )
    if sys.version_info[0:3] <= (2,4,3) :
        t1 = datetime.datetime( *map(int,rcf.get('jobstep.timerange.start').replace('-',' ').replace(':',' ').split()) )
        t2 = datetime.datetime( *map(int,rcf.get('jobstep.timerange.end'  ).replace('-',' ').replace(':',' ').split()) )
    else :
        t1 = datetime.datetime.strptime( rcf.get('jobstep.timerange.start'), tfmt )
        t2 = datetime.datetime.strptime( rcf.get('jobstep.timerange.end'  ), tfmt )
    #endif
    istart      = rcf.get( 'istart', 'int' )
    savedir     = rcf.get( 'savedir' )
    start_3_dir = rcf.get( 'start.3.dir', default=savedir )
    output_dir  = rcf.get( 'output.dir' )

    # start from save files after restart:
    if step > 1 : istart = 3

    # fill lines of file:
    lines = []
    lines.append( '\n' )
    lines.append( '!\n' )
    lines.append( '! TM5 Rcfile with runtime settings.\n' )
    lines.append( '! Created by : submit_tm5_setup_rcfile.py\n' )
    lines.append( '!\n' )
    lines.append( '\n' )
    lines.append( '! start time:\n' )
    lines.append( 'year1   : %i\n' % t1.year  )
    lines.append( 'month1  : %i\n' % t1.month )
    lines.append( 'day1    : %i\n' % t1.day   )
    lines.append( 'hour1   : %i\n' % t1.hour  )
    lines.append( 'minu1   : 0\n' )
    lines.append( 'sec1    : 0\n' )
    lines.append( '\n' )
    lines.append( '! end time:\n' )
    lines.append( 'year2   : %i\n' % t2.year  )
    lines.append( 'month2  : %i\n' % t2.month )
    lines.append( 'day2    : %i\n' % t2.day   )
    lines.append( 'hour2   : %i\n' % t2.hour  )
    lines.append( 'minu2   : 0\n' )
    lines.append( 'sec2    : 0\n' )
    lines.append( '\n' )
    lines.append( '! defines how to fill initial state:\n' )
    lines.append( 'istart  : %i\n' % istart )
    lines.append( '\n' )
    lines.append( '! directory with save files:\n' )
    lines.append( 'savedir : %s\n' % start_3_dir )
    lines.append( '\n' )
    lines.append( '! output directory:\n' )
    lines.append( 'outputdir  : %s\n' % output_dir )
    lines.append( '\n' )

    # write:
    f = open( 'tm5_runtime.rc', 'w' )
    for line in lines : f.write( line )
    f.close()

    # ok
    return

#enddef

