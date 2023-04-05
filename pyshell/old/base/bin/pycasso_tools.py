#!/usr/bin/env python

import os
import logging

"""
PYCASSO tools
"""


#-------------------------------------------------
# routines
#-------------------------------------------------


def write_text_file( fname, text ) :

    """
    Write a text file.
    Arguments:
      fname   : target file name
      text    : list of strings, a line should end with '\n'
    """
    
    # external:
    import os
    import logging
    
    # info ...
    logging.debug( '    write %s ...' % fname )

    # write new text:
    f = open( fname, 'w' )
    f.writelines( text )
    f.close()

    # ok
    return
    
#enddef


# ***


def diff_text_files( fname1, fname2 ) :
    with open(fname1, 'r') as f1:
        with open(fname2, 'r') as f2:
            return f1.readlines() == f2.readlines()


def update_text_file( fname, newtext ) :

    """
    Replace a file by a new text if the later differs
    from the current content.
    Arguments:
      fname      : target file name
      newtext    : list of strings, a line should end with '\n'
    """
    
    # external:
    import os
    import logging
    
    # file exists alread?
    if os.path.exists(fname) :
        # read current content:
        f = open( fname, 'r' )
        oldtext = f.readlines()
        f.close()
        # differences ?
        rewrite = newtext != oldtext
        ## info ...
        #for iline in range(len(oldtext)) :
        #    if iline < len(newtext) :
        #        if oldtext[iline] != newtext[iline] :
        #            logging.debug( '    first different lines:' )
        #            logging.debug( '      old: %s' % oldtext[iline] )
        #            logging.debug( '      new: %s' % newtext[iline] )
        #            break
        #        #endif
        #    #endif
        ##endfor
    else :
        # no file yet, always rewrite:
        rewrite = True
    #endif
    
    # write file ?
    if rewrite :
        # for info message:
        stat = 'replace'
        # write new text:
        f = open( fname, 'w' )
        f.writelines( newtext )
        f.close()
    else :
        # for info message:
        stat = 'keep'
    #endif
    
    # info ...
    logging.debug( '    %-8s %-40s' % (stat,fname) )

    # ok
    return
    
#enddef


# ***


def modify_text_file( fname, key, value ) :

    """
    Modify a text file by replacing a key by a new value.
    Arguments:
      fname      : text file name
      key        : value to be replaced
      value      : replacement value
    """
    
    # read file into list of lines:
    f = open( fname, 'r' )
    lines = f.readlines()
    f.close()

    # copy while replacing key with value:
    newlines = []
    for line in lines :
        # replace key with value:
        line = line.replace( key, value )
        # add:
        newlines.append(line)
    #endfor
    
    # write again:
    f = open( fname, 'w' )
    f.writelines(newlines)
    f.close()

    # ok
    return
    
#enddef


# ***


def CreateDirs( dirname, forceclean=False ) :

    """
    Create a directory and report success.
    """

    # external:
    import os
    import shutil
    import logging
    
    # already present ?
    if os.path.isdir( dirname ) :
        # remove existing directory ?
        if forceclean:
            # info ...
            logging.info( 'remove existing %s ...' % dirname )
            # use the chainsaw:
            shutil.rmtree( dirname )
        #endif
    #endif
    
    # not present (anymore) ?
    if not os.path.isdir( dirname ) :
        # info ...
        logging.info( 'create new directory %s ...' % dirname )
        # create directory including parents if necesary:
        os.makedirs(dirname)
    #endif
    
    # ok
    return None

#enddef

            
#-------------------------------------------------
# end
#-------------------------------------------------
