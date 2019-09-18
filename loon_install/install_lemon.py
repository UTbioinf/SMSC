#!/usr/bin/env python

import os
import subprocess
import multiprocessing

def command_path(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def install_lemon(lemon_dir, libexec):
    old_dir = os.getcwd()
    os.chdir( lemon_dir )

    subprocess.check_call( ["tar", "xf", "lemon-1.3.1.tar.gz"] )
    subprocess.check_call( ["mkdir", "-p", "build_lemon"] )
    os.chdir( "build_lemon" )
    install_dir = os.path.join(libexec, "lemon")
    subprocess.check_call( ["cmake", "../lemon-1.3.1",
            "-DCMAKE_INSTALL_PREFIX={}".format(install_dir),
            "-DCMAKE_BUILD_TYPE=Release",
            "-DLEMON_ENABLE_COIN=NO",
            "-DLEMON_ENABLE_ILOG=NO",
            "-DLEMON_ENABLE_GLPK=NO",
            "-DLEMON_DOC_SOURCE_BROWSER=NO",
            "-DTEST_WITH_VALGRIND=NO",
            "-DLEMON_DOC_USE_MATHJAX=NO"] )
    n_cpu = str(multiprocessing.cpu_count())
    subprocess.check_call( ["make", "-j", n_cpu] )
    subprocess.check_call( ["make", "-j", n_cpu, "check"] )
    subprocess.check_call( ["make", "install"])
    os.chdir( old_dir )
    return os.path.join(install_dir, "share", "lemon", "cmake")

def check_lemon(libexec):
    lemon_path = command_path("dimacs-solver")
    if lemon_path:
        lemon_path = os.path.dirname( lemon_path )
        lemon_path = os.path.normpath( os.path.join(lemon_path, "..", "share", "lemon", "cmake") )
        if os.path.isfile( os.path.join(lemon_path, "LEMONConfig.cmake") ):
            return lemon_path
    lemon_path = os.path.normpath( os.path.join(libexec, "lemon", "share", "lemon", "cmake") )
    if os.path.isfile( os.path.join(lemon_path, "LEMONConfig.cmake") ):
        return lemon_path
    else:
        return ""

if __name__ == "__main__":
    print check_lemon()
    print install_lemon( os.getcwd() )
