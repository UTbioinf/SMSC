#!/usr/bin/env python2

from loon_install import multitool_install
import os

developer_mode = False

def add_lemon_parameters(obj, cmd_args):
    import sys
    sys.path.append(obj.source["cpp"])
    import install_lemon
    lemon_path = install_lemon.check_lemon()
    if not lemon_path:
        lemon_path = install_lemon.install_lemon( obj.source["cpp"] )
    cmd_args.append("-DLEMON_DIR={}".format( lemon_path )) # LEMON_DIR will be used by cmake directly to find the LEMONConfig.cmake file

def main():
    default_config = {
            "LOGGER_LEVEL": "3",
            "cpp_commands_extra": []
        }

    if developer_mode:
        try:
            import developer
            print "===> developer_mode = ON"
            default_config = developer.config
        except ImportError:
            print "===> developer_mode is forced off"

    multitool_install.MultiToolInstall(
            project = "smsc",
            license = "LICENSE.txt",
            source_dir = os.path.dirname( __file__ ),

            version = "2.2",
            start_hex = "7994245c2697eb09bac1d8b7dcab4c190e17bcf1",
            cpp = {
                    "commands": [(add_lemon_parameters, )] + default_config["cpp_commands_extra"]
                }
        ).run()
    

if __name__ == "__main__":
    main()

