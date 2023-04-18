#! /usr/bin/env python
import re
import os
import sys
from loguru import logger
from pathlib import Path


class RcFile:
    def __init__(self, filename=None, marks=('${', '}')):
        self.outfile = []
        self.values = {}
        self.sources = {}
        self.curfile = None
        self.rootdir = None
        self.trace = []
        if filename is not None :
            self.readfile(filename, marks=marks)

    def readfile(self, filename, marks=('${', '}')):
        logger.debug(f'reading rcfile {filename}')
        self.curfile = Path(filename)
        self.rootdir = self.curfile.parent

        with open(filename, 'r') as f:
            inpfile = f.readlines()
            inpfile.extend(self.outfile)

        # create traceback info:
        inptrace = []
        for jline in range(len(inpfile)) :
            inptrace.append(f'"{filename}", line {jline+1}')

        # pass counter:
        ipass = 1

        # loop until all substitutions and inclusions are done:
        while True :
            self.outfile = []
            self.trace = []

            # init current line:
            line = ''
            something_done = False
            something_to_be_done = False
            unresolved_lines = []
            keys_with_unresolved_value = []
            undefined_keys = []

            # stack for conditional evaluation;
            # each element is a tuple with elements:
            #   resolved (boolean) true if the if-statement is evaluated
            #   flag     (boolean) true if the lines below the statement
            #              are to be included
            #   anyflag    (boolean) used to check if any of the 'if' or 'elif' conditions
            #              in this sequence evaluated to True
            #   line     (char) description of the line for messages
            ifstack = []

            # loop over lines in input file:
            iline = -1
            for inpline in inpfile :

                # line counter:
                iline = iline + 1

                # cut current traceback info from list:
                linetrace_curr = inptrace.pop(0)

                # set full traceback info:
                if line.endswith('\\') :
                    qfile, qlinenrs = linetrace.split(',')
                    qnrs = qlinenrs.replace('lines', '').replace('line', '')
                    if '-' in qnrs :
                        qnr1, qnr2 = qnrs.split('-')
                    else :
                        qnr1, qnr2 = qnrs, qnrs
                    linetrace = '%s, lines %-9s' % (qfile, '%i-%i' % (int(qnr1), int(qnr2)+1))
                else :
                    linetrace = linetrace_curr

                inpline = inpline.strip()

                if not inpline :
                    self.outfile.append('\n')
                    self.trace.append(linetrace)
                    line = ''
                    continue

                # skip comment:
                if inpline.startswith('!') :
                    self.outfile.append(f'{inpline}\n')
                    self.trace.append(linetrace)
                    line = ''
                    continue

                # continuation lines
                # current line has continuation mark '\' at the end ?
                # then add this input line:
                if line.endswith('\\') :
                    # remove continuation character, add input line:
                    line = line[:-1] + inpline
                else :
                    line = inpline

                # continuation mark ? then next line of input file:
                if line.endswith('\\') :
                    continue

                # conditional settings (1)

                # is this the begin of a new condition ?
                if line.startswith("#if"):
                    ifstack.append((False, True, False, linetrace))

                if line.startswith("#elif"):
                    if len(ifstack) == 0 :
                        logger.error(f'found orphin "#elif" in {linetrace}')
                        raise Exception
                    # remove current top from stack:
                    resolved, flag, anyflag, msg = ifstack.pop()
                    # did one of the previous #if/#elif evaluate to True already ?
                    if resolved and anyflag :
                        # push to stack that this line resolved to False :
                        ifstack.append((True, False, anyflag, linetrace))
                        self.outfile.append(f'{line}\n')
                        self.trace.append(linetrace)
                        continue
                    else :
                        # push temporary flag to stack, will be replaced after evaluation of condition:
                        ifstack.append((False, True, anyflag, linetrace))

                if line.startswith("#else") :
                    if len(ifstack) == 0 :
                        logger.error(f'found orphin "#else" in {linetrace}')
                        raise Exception
                    # remove current top from stack:
                    resolved, flag, anyflag, msg = ifstack.pop()
                    # get higher level settings:
                    if len(ifstack) > 0 :
                        resolved_prev, flag_prev, anyflag_prev, msg_prev = ifstack[-1]
                    else :
                        flag_prev = True
                    # should next lines be included ?
                    # reverse flag, take into acount higher level:
                    new_flag = (not flag) and (not anyflag) and flag_prev
                    # push to stack:
                    ifstack.append((resolved, new_flag, False, linetrace))
                    self.outfile.append(f"{line}\n")
                    self.trace.append(linetrace)
                    continue

                # is this the end of a condition ?
                if line.startswith("#endif") :
                    if len(ifstack) == 0 :
                        logger.error(f'found orphin "#endif" in {linetrace}')
                        raise Exception
                    _ = ifstack.pop()
                    self.outfile.append(f'{line}\n')
                    self.trace.append(linetrace)
                    continue

                # within if-statements ?
                if len(ifstack) > 0 :
                    resolved, flag, anyflag, msg = ifstack[-1]
                    # already resolved ? then check if this line should be skipped:
                    if resolved and (not flag) :
                        # skip this line, but include commented version in output:
                        self.outfile.append(f'!{line}\n')
                        self.trace.append(linetrace)
                        continue

                if line.startswith("#eval"):
                    raise NotImplementedError

                # ensure that common marks are evaluated correctly:
                start_mark = marks[0].replace('{', '\{').replace('<', '\<').replace('$', '\$')
                close_mark = marks[1].replace('}', '\}').replace('>', '\>')

                # set syntax of keywords to be matched, e.g. '${...}' :
                pattern = start_mark + '[A-Za-z0-9_.]+' + close_mark

                # make a regular expression that matches all variables:
                rc_varpat = re.compile(pattern)

                # search all matching paterns:
                pats = re.findall(rc_varpat, line)

                # counter for unexpanded substitutions:
                ntobedone = 0

                # loop over matches:
                for pat in pats :
                    key = pat.lstrip(start_mark).rstrip(close_mark)
                    if key in self.values:
                        val = self.values[key]
                        line = line.replace(pat, val)
                        something_done = True
                    elif key in os.environ:
                        val = os.environ[key]
                        line = line.replace(pat, val)
                        something_done = True
                    elif key == 'pid' :
                        val = str(int(os.getpid()))
                        line = line.replace(pat, val)
                        something_done = True
                    elif key == 'script' :
                        script, ext = os.path.splitext(os.path.basename(sys.argv[0]))
                        line = line.replace(pat, script)
                        something_done = True
                    else :
                        ntobedone = ntobedone + 1
                        if key not in undefined_keys :
                            undefined_keys.append(key)
                        continue

                if ntobedone > 0 :
                    self.outfile.append(line)
                    self.trace.append(linetrace)
                    something_to_be_done = True
                    unresolved_lines.append(f'{linetrace} | {line}')
                    if ':' in line :
                        qkey, qvalue = line.split(':', 1)
                        qkey = qkey.strip()
                        if (' ' not in qkey) and (start_mark not in qkey) and (not qkey.startswith('#')) :
                            if qkey not in keys_with_unresolved_value :
                                keys_with_unresolved_value.append(qkey)
                    continue

                if line.startswith("#include") :
                    inc_file = line.lstrip("#include").strip()
                    if not os.path.exists(inc_file) :
                        inc_file = os.path.join(self.rootdir, inc_file)
                        logger.debug(f'Added rootdir to requested include: {inc_file}')
                    if not os.path.exists(inc_file) :
                        logger.error('include file not found : {inc_file}')
                        logger.error(linetrace)
                        raise IOError
                    with open(inc_file, 'r') as inc_f:
                        inc_rc = inc_f.readlines()

                    self.outfile.append(f'! >>> {inc_file} >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n')
                    self.outfile.extend(inc_rc)
                    self.outfile.append(f'! <<< {inc_file} <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n')
                    self.trace.append(linetrace)
                    for jline in range(len(inc_rc)) :
                        self.trace.append(f'"{inc_file}", line {jline+1:.0f}')
                    self.trace.append(linetrace)
                    something_done = True
                    something_to_be_done = True
                    continue

                # evaluate conditional expressions:
                if line.startswith("#if") or line.startswith("#elif") :
                    # remove leading mark, what remains is logical expression:
                    expression = line.lstrip("#elif").strip()
                    if expression.endswith(':') :
                        expression = expression.rstrip(':').strip()
                    try :
                        flag = eval(expression)
                    except Exception as e:
                        logger.error(f'could not evaluate expression "{expression}" in {linetrace}')
                        raise e

                    # remove temporary top added before during this pass:
                    tmp_statement, tmp_flag, tmp_anyflag, tmp_msg = ifstack.pop()

                    # extract current top if necessary:
                    if len(ifstack) > 0 :
                        dummy_statement, prev_flag, dummy_anyflag, dummy_msg = ifstack[-1]
                    else :
                        prev_flag = True

                    # should next lines be included ?
                    new_flag = prev_flag and tmp_flag and flag
                    # any if/elif evaluated to true in this sequence ?
                    new_anyflag = tmp_anyflag or new_flag
                    # add to stack, now resolved, take into accout current flag:
                    ifstack.append((True, new_flag, new_anyflag, linetrace))
                    self.outfile.append(line + '\n')
                    self.trace.append(linetrace)
                    continue

                # error message
                if line.startswith("#error"):
                    raise NotImplementedError

                # checks
                # common mistake ...
                if line.startswith('#') :
                    logger.error(f'line {line} in rcfile starts with "#" but is not an "#include" or other special line;')
                    logger.error('if it is supposed to be comment, please start with "!" ...')
                    logger.error(linetrace)
                    raise RuntimeError

                # check ...
                if ':' not in line :
                    logger.error(f'key/value line should contain a ":"')
                    logger.error(linetrace)
                    raise RuntimeError

                # add to output
                self.outfile.append(line + '\n')
                self.trace.append(linetrace)

                # add key/value pair
                # not if inside an unresolved if-statement ...
                if len(ifstack) > 0 :
                    # get top values:
                    resolved, flag, anyflag, msg = ifstack[-1]
                    # not resolved yet ? then continue:
                    if not resolved :
                        continue

                # split in key and value;
                key, val = line.split(':', 1)

                # remove comment from value:
                if '!' in val:
                    if '\!' not in val:
                        val, comment = val.split('!')
                    val = val.replace('\!', '!')

                # remove spaces:
                key = key.strip()
                val = val.strip()

                # already defined ?
                if key in self.values:
                    # this will occure often after the first pass since
                    # the keys are resolved again and again ;
                    # therefore, only complain if this definition is read
                    # from a different line :

                    if linetrace != self.sources[key] and self.sources[key] != 'External':
                        logger.error(f'duplicated key \"{key}\" found:')
                        logger.error(f'first definition in {self.sources[key]} is:')
                        logger.error(f'  {key}  : {self.values[key]}')
                        logger.error(f'second definition in {linetrace.strip()} is:')
                        logger.error(f'  {key}  : {val}')
                        raise RuntimeError
                else :
                    self.values[key] = val
                    self.sources[key] = linetrace
                    something_done = True

            # check ...
            if len(ifstack) > 0 :
                logger.error('unterminated if-statement ; current stack:')
                for resolved, flag, anyflag, msg in ifstack:
                    logger.error(msg)
                raise RuntimeError

            # check ...
            if something_to_be_done :
                if not something_done :
                    logger.error('Could not resolve the following lines in rcfile(s):')
                    for uline in unresolved_lines :
                        logger.error(f'    {uline}')
                    logger.error('  Undefined key(s):')
                    for ukey in undefined_keys :
                        # do not list them if they are undefined because the value
                        # depends on other undefined keys:
                        if ukey not in keys_with_unresolved_value :
                            logger.error(f'    {ukey}')
                            # loop over unresolved lines to see in which the key is used:
                            for uline in unresolved_lines :
                                # search for  '${key}' pattern:
                                if marks[0] + ukey + marks[1] in uline :
                                    logger.error(f'      {uline}')
                    raise RuntimeError
            else :
                break

            # for safety ...
            if ipass == 100 :
                logger.error(f'resolving rc file has reached pass 100 ; something wrong ?')
                raise RecursionError

            # new pass:
            ipass += 1
            inpfile = self.outfile
            inptrace = self.trace

    def has_key(self, key) :
        return key in self.values

    def keys(self) :
        return self.values.keys()

    def get(self, key, totype='', default=None) :
        if key in self.values:
            value = self.values[key]
            if totype == 'bool' :
                if value in ['T', 'True', 'yes', '1'] :
                    return True
                elif value in ['F', 'False', 'no', '0'] :
                    return False
                else :
                    logger.error("value of key '%s' is not a boolean : %s" % (key, str(value)))
                    raise KeyError
            elif len(totype) > 0 :
                return eval('%s(%s)' % (totype, value))
            return value
        else :
            if default is not None :
                return default
            else :
                # something wrong ...
                logger.error("key '%s' not found in '%s' and no default specified" % (key, self.curfile))
                raise KeyError

    def replace(self, key, val) :
        found = False
        for iline in range(len(self.outfile)) :
            line = self.outfile[iline]
            if ':' not in line :
                continue
            k, v = line.split(':', 1)
            if k.strip() == key :
                self.outfile[iline] = '%s : %s\n' % (k.strip(), str(val))
                self.values[key] = val
                found = True
                break
        if not found :
            logger.error('could not replace key : %s' % key)
            raise KeyError

    def replace_add(self, key, val):
        if key in self.values:
            self.replace(key, val)
        else:
            self.add(key, val)

    def add(self, key, val, comment='') :
        if len(comment) > 0 :
            raise NotImplementedError
            # self.outfile.append('! %s\n' % comment)
        self.outfile.append('%s : %s\n' % (key.strip(), str(val)))

        # add to dictionairy:
        self.values[key.strip()] = val
        self.sources[key] = 'External'
        
    def substitute(self, line, marks=('${', '}')) :
        start_mark = marks[0].replace('{', '\{').replace('<', '\<').replace('$', '\$')
        close_mark = marks[1].replace('}', '\}').replace('>', '\>')
        pattern = start_mark + '[A-Za-z0-9_.]+' + close_mark
        rc_varpat = re.compile(pattern)
        pats = re.findall(rc_varpat, line)
        for pat in pats :
            key = pat.lstrip(start_mark).rstrip(close_mark)
            if key in self.values:
                val = self.values[key]
                line = line.replace(pat, val)
        return line

    def WriteFile(self, filename):
        with open(filename, 'w') as f:
            f.writelines(self.outfile)


def read(rcfilename) :
    return RcFile(rcfilename).values


def write(filename, rcdict):
    with open(filename, 'w') as f:
        for k, v in rcdict.items():
            f.write('%-20s:%s\n' % (k, v))
