#!/usr/bin/env -S python3 

import os,time

PROC = None

#
#
#
def connect( SERVER ) :
    global PROC
    from subprocess import Popen, PIPE
    ssh_cmd = [ "ssh", "-T", SERVER ]
    PROC = Popen(ssh_cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE,  universal_newlines=True)
    os.set_blocking(PROC.stdout.fileno(), False)
    os.set_blocking(PROC.stderr.fileno(), False)
    print(f"#SSH# Connected to {SERVER}.")
    return PROC

#
#
#
def cmd( cmd, quiet=1, timeout=5, deb=0 ) :
    global PROC
    #
    # PROCEDURE : Run command in remote host
    #
    print(f"#SSH# Running command: $ \"{cmd}\" ...")
    PROC.stdin.write(f"{cmd}\n")
    PROC.stdin.flush()
    
    ENDTAG = "THIS IS AN END TAG. THIS IS AN END TAG. THIS IS AN END TAG."
    # Writes in the error pipe to avoir missing errors from the command.
    PROC.stdin.write(f"echo {ENDTAG} 1>&2\n")
    PROC.stdin.flush()

    stdout = ""
    stderr = ""

    #
    # PROCEDURE : Wait for the proces to start outputting
    #
    i=0
    while True :
        stderr_line = PROC.stderr.readline()
        stdout_line = PROC.stdout.readline()
        # Keep going if we have seen somethin in STDOUT or STDERR
        if stderr_line : break
        if stdout_line : break

        i+=1 
        if i > 100 : 
            print("#SSH# E: Process did not output as expected after 10 seconds.")
            return "", ""
        time.sleep(0.1)

    #
    # PROCEDURE : Collect all the output - break with ENDTAG or 1s idle
    #
    i=0
    while True :   # Still working
        # Collect stdout
        if stdout_line :
            if not quiet : print(stdout_line.strip())

            stdout += stdout_line
            i=0
            stdout_line = PROC.stdout.readline()
            continue  ## STDOUT has the higher priority. Collect all those before checking for errors

        # Collect stderr
        if stderr_line :
            # The command is over and we have seen the endtag.
            if stderr_line.strip() == ENDTAG :
                if deb : print("#SSH#debug# Found end tag.")
                break 
            if not quiet : print(stderr_line.strip())

            stderr += stderr_line
            stderr_line = PROC.stderr.readline()
            i=0
            continue

        # Empty. Give one second
        i+=1 
        if i > timeout/0.1 : 
            print(f"#SSH# {timeout} second without output - ssh command TIMEOUT.")
            break
        time.sleep(0.1)


    # DONE.
    return stdout, stderr

## FACILITIES TO CONTROL THE CLUSTER SUBMISSIONS
#
# Return the number of processes running currently
#
def n_jobs_running( deb = 0 ) :
    # Flush outputs
    n_tries = 0
    while n_tries < 10 :
        try :
            o,e = cmd( "squeue -u bfq9 -h -t pending,running -r | wc -l", deb=deb )
            print(f"#SSH# Currently running jobs for BFQ9: {int(o)}.")
            n = int(o)
            break
        except :
            print(f"#SSH# Failed to fetch the number of jobs running (n_tries={n_tries}). STDOUT='{o}'")
            time.sleep(0.5)
            n_tries += 1
    
    if n_tries >= 10 :
        print(f"#SSH# Failed to get the number of jobs running. Returning default value 1000 to keep things going...")
        n=1000

    return n


## EXAMPLE OF USAGE

# # MAIN
# p = chimas_ssh.connect()
# msg,err = chimas_ssh.cmd( p, "ls", 1 )
# print( msg )

# msg,err = ssh_cmd( p, "cd temp ; ls",1 )
# print( msg )
# print( err )

# msg,err = ssh_cmd( p, "echo AAAA 1>&2" )
# print( err )

