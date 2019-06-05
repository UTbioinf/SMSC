#!/usr/bin/env python2

from loon_install import multitool_install
from loon_install import util
import shutil
import os
import subprocess

developer_mode = False

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

def add_lemon_parameters(obj, cmd_args):
    obj.logger.info("Check lemon")
    # lemon will be installed into "<prefix>/libexec/lemon/"
    from loon_install import install_lemon
    lemon_path = install_lemon.check_lemon( obj.install["libexec"])
    if not lemon_path:
        obj.logger.info("Install lemon")
        lemon_path = install_lemon.install_lemon( os.path.join(obj.source_dir, "depends"), obj.install["libexec"] )
    else:
        obj.logger.info("Lemon has been installed")
    cmd_args.append("-DLEMON_DIR={}".format( lemon_path )) # LEMON_DIR will be used by cmake directly to find the LEMONConfig.cmake file

def check_blasr(obj, cmd_args):
    obj.logger.info("Check Blasr")
    if not command_path("blasr"):
        obj.logger.warning("BLASR is not available! Please install it before running smsc")


def install_mummer(obj, cmd_args, mummer_dirname):
    # mummer will be installed into "<prefix>/libexec/mummer/<mummer_dirname>"
    obj.logger.info("Check Mummer")
    mummer_dir = os.path.join(obj.install["libexec"], "mummer")
    if command_path( os.path.join(mummer_dir, mummer_dirname, "mummer") ):
        obj.logger.info("Mummer has been installed")
        return

    obj.logger.info("Install Mummer")
    util.mkdir_p( mummer_dir )
    shutil.copy( os.path.join(obj.source_dir, "depends", "{}.tar.gz".format( mummer_dirname )), mummer_dir )
    
    old_dir = os.getcwd()
    os.chdir( mummer_dir )
    
    subprocess.check_call( ["tar", "xf", "{}.tar.gz".format( mummer_dirname )] )
    os.chdir( mummer_dirname )
    util.mkdir_p( "aux_bin" )
    subprocess.check_call( ["make", "check"] )
    subprocess.check_call( ["make", "install"] )

    os.chdir( old_dir )

def install_mhap(obj, cmd_args, mhap_jar):
    # mhap will be installed into "<prefix>/libexec/mhap/"
    obj.logger.info("Check mhap")
    mhap_dir = os.path.join(obj.install["libexec"], "mhap")

    if os.path.isfile( os.path.join(mhap_dir, mhap_jar) ):
        obj.logger.info("Mhap has been installed")
        return

    obj.logger.info("Install mhap")
    util.mkdir_p( mhap_dir )
    shutil.copy( os.path.join(obj.source_dir, "depends", "{}.gz".format( mhap_jar )), mhap_dir )

    old_dir = os.getcwd()
    os.chdir( mhap_dir )

    subprocess.check_call( ["gzip", "-d", "{}.gz".format( mhap_jar )] )

    os.chdir( old_dir )

def install_muscle(obj, cmd_args, muscle_exe):
    # muscle will be installed into "<prefix>/libexec/muscle/"
    obj.logger.info("Check muscle: (linux64 version only. Please contact the author if you need other versions)")
    muscle_dir = os.path.join( obj.install["libexec"], "muscle" )

    if command_path(os.path.join(muscle_dir, muscle_exe)):
        obj.logger.info("Muscle has been installed")
        return

    obj.logger.info("Install Muscle")
    util.mkdir_p( muscle_dir )
    shutil.copy( os.path.join(obj.source_dir, "depends", "muscle", "{}.tar.gz".format( muscle_exe )), muscle_dir )

    old_dir = os.getcwd()
    os.chdir( muscle_dir )
    
    subprocess.check_call( ["tar", "xf", "{}.tar.gz".format( muscle_exe )] )

    os.chdir( old_dir )
    
def add_config(obj, cmd_args, mummer_dirname, mhap_jar, muscle_exe):
    obj.logger.info("Generate config.sh")
    mummer_path = os.path.join( obj.install["libexec"], "mummer", mummer_dirname)
    mhap_path = os.path.join( obj.install["libexec"], "mhap", mhap_jar )
    muscle_path = os.path.join( obj.install["libexec"], "muscle", muscle_exe )
    with open(os.path.join(obj.source["cpp"], "smsc", "config.sh"), "w") as fout:
        fout.write("""#!/bin/bash

# configure mummer root directory, full path is a must
mummer_path={}

# configure mhaps, the executable jar is included. Full path is a must
mhap_path={}
mhap_mem=32g

# configure muscle, the executable is included. Full path is a must
MUSCLE_EXEC={}
max_iters=5
        """.format(mummer_path, mhap_path, muscle_path))

    with open(os.path.join(obj.source["cpp"], "smsc", "version.h"), "w") as fout:
        fout.write("""#ifndef __SMSC_VERSION_H
#define __SMSC_VERSION_H

#define SMSC_VERSION "{}"

#endif
        """.format( obj.get_version().strip() ))

def main():
    default_config = {
            "cpp_commands_extra": []
        }

    if developer_mode:
        try:
            import developer
            print "===> developer_mode = ON"
            default_config = developer.config
        except ImportError:
            print "===> developer_mode is forced off"

    mummer_dirname = "MUMmer3.23"
    mhap_jar = "mhap-2.1.jar"
    muscle_exe = "muscle3.8.31_i86linux64"
    multitool_install.MultiToolInstall(
            project = "smsc",
            license = "LICENSE",
            source_dir = os.path.dirname( __file__ ),

            version = "2.2",
            start_hex = "7994245c2697eb09bac1d8b7dcab4c190e17bcf1",
            cpp = {
                    "commands": [
                        (check_blasr,),
                        (add_lemon_parameters, ),
                        (install_mummer, mummer_dirname),
                        (install_mhap, mhap_jar),
                        (install_muscle, muscle_exe),
                        (add_config, mummer_dirname, mhap_jar, muscle_exe)
                        ]
                }
        ).run()
    

if __name__ == "__main__":
    main()

