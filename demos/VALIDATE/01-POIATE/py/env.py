
import os, sys, io

DL = int ( os.getenv("CHDEBUG", 1) )

#
# Fromat string in lines.
#
def _print( *obj, sep=' ', end='\n', pre="" ) :
    s = sep.join(map(str,obj))
    s = s.splitlines()
    ret = ""
    for l in s :
        if len(pre) : ret += "# " + pre + ": ";
        ret += l + end;
    print(ret, end='')
    sys.stdout.flush()


#
# Log functions
#
def dlog( dl, *obj, sep=' ', end='\n', file=sys.stdout, flush=False, pre=True ) :
    if DL >= dl :
        pt = f"Deb[{dl}]" if pre else "";
        _print(*obj, sep=sep, end=end, pre=pt )

def ilog( *obj, sep=' ', end='\n', file=sys.stdout, flush=False, pre=True ) :
    pt = "Info  " if pre else "";
    _print(*obj, sep=sep, end=end, pre=pt )

def wlog( *obj, sep=' ', end='\n', file=sys.stdout, flush=False, pre=True ) :
    pt = "Warn  " if pre else "";
    _print(*obj, sep=sep, end=end, pre=pt )

def elog( *obj, sep=' ', end='\n', file=sys.stdout, flush=False, pre=True ) :
    pt = "Error " if pre else "";
    _print(*obj, sep=sep, end=end, pre=pt )

def flog( *obj, sep=' ', end='\n', file=sys.stdout, flush=False, pre=True ) :
    pt = "Fail  " if pre else "";
    _print(*obj, sep=sep, end=end, pre=pt )
    exit(-1)

##
#
# Creates a directory and do not throw error if it already exists.
#
def mkdir(d) :
    dlog(1, f"Creating directory '{d}'...", end="")
    try:
        os.mkdir(d)
        dlog(1," [Done]", pre=False)
    except: 
        dlog(1," [Existed]", pre=False)
        pass

##
#
# Delete file if exists
#
def file_remove(d) :
    dlog(1, f"Deleting file'{d}'...", end="")

    # First try to get rid of the file
    try:
        import shutil
        shutil.move(path, "/tmp/__trash_")
        os.remove("/tmp/__trash_")
        dlog(1," [Done]", pre=False)
    except: 
        dlog(1," [Did not exist]", pre=False)
        pass


#
# The program dies if infiles are older than the outfiles
#
def die_if_uptodate( in_files, out_files ) :
    INM  = [ os.path.getmtime(x) for x in in_files ]

    may_skip = True
    for f in out_files : 
        if not os.path.isfile(f) : may_skip = False

    if may_skip :
        OUTM = [ os.path.getmtime(x) for x in out_files ]
        if max(INM) < min(OUTM)  :
          ilog()
          ilog("########################################")
          ilog("Files are up to date. Exitting ...")
          ilog("########################################")
          exit(0)

