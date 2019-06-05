#!/usr/bin/env python2

import subprocess
import errno
import os
import shutil
import argparse
import distutils.spawn
import logging
import util

class MultiToolInstall(object):
    def __init__(self, project, version, source_dir, license = None, start_hex=None, python=None, cpp=None, makefile=None):
        self.project = project
        self.version = version
        self.start_hex = start_hex

        self.py = python
        self.cpp = cpp
        self.makefile = makefile

        # parse cmd parameters
        self.args = self.parse_args()
        if self.args.verbose:
            logging.basicConfig(format="[%(asctime)s] [%(levelname)s]: %(message)s", level=logging.INFO)
        else:
            logging.basicConfig(format="[%(asctime)s] [%(levelname)s]: %(message)s", level=logging.ERROR)

        # init logger
        self.logger = logging.getLogger()

        # init directories
        self.source_dir = os.path.abspath( os.path.realpath(source_dir) )
        if self.cpp and self.args.build_debug:
            self.build_dir = os.path.join(self.source_dir, "build", "Debug")
        else:
            self.build_dir = os.path.join(self.source_dir, "build", "Release")
        if self.args.prefix:
            self.install_dir = os.path.abspath( os.path.realpath( self.args.prefix ) )
        else:
            self.install_dir = os.path.join( self.source_dir, "build", "local" )
            
        self.install = {}
        install_structure = {"bin": ["bin"],
                             "lib": ["lib"],
                             "libexec": ["libexec"],
                             "libpython": ["lib", "python", "{}_lib".format(self.project)],
                             "libexec_bin": ["libexec", "bin"],
                             "libexec_mk": ["libexec", "makefiles"]
                            }
        for key, val in install_structure.iteritems():
            self.install[ key ] = os.path.join(self.install_dir, *val)        

        self.source = {}
        source_structure = {"python_root": ["python"],
                            "python": ["python", "src"],
                            "pylib":  ["python", "lib"],
                            "cpp": ["cpp"],
                            "makefile": ["makefiles"]}
        for key, val in source_structure.iteritems():
            self.source[ key ] = os.path.join(self.source_dir, "src", *val)

        # Other parameters
        self.brew = None

        # license
        self.license = os.path.join(self.source_dir, license)

    def parse_args(self):
        parser = argparse.ArgumentParser(description = "Installing script")
        parser.add_argument("--prefix", help="Installing prefix (default: <repo dir>/build/local)")
        if self.cpp != None:
            parser.add_argument("--build-debug", action="store_true", help="Set CPP build type as Debug rather than Release")
            cpp_parameters = self.cpp.get("parameters", []) # cmake parameters. e.g., abc = 3 will be set as -Dabc=3 when invoking the smake command
            for each in cpp_parameters:
                if len(each) == 1:
                    parser.add_argument("--{}".format(each[0]), action="store_true", help="CPP parameter")
                else:
                    parser.add_argument("--{}".format(each[0]), help="CPP parameter")
        parser.add_argument("-v", "--verbose",  action="store_true", help="Verbose")
        return parser.parse_args()        

    def get_version(self):
        if self.start_hex == None:
            return self.version
        if os.path.isdir( os.path.join(self.source_dir, ".git") ):
            return "{}.{}".format(self.version, util.git_count(self.start_hex))
        else:
            return self.version

    def install_version(self):
        with open(os.path.join(self.install["bin"], "version.txt"), "w") as fout:
            this_version = self.get_version()
            self.logger.info("Install version {}".format( this_version ))
            fout.write(this_version)

    def install_license(self):
        if self.license:
            self.logger.info("Install license")
            shutil.copy(self.license, self.install["bin"])
            

    # Install python packages
    def install_python_preinstall(self):
        pass

    def install_python_postinstall(self):
        pass

    def install_python(self):
        if self.py == None: return

        self.install_python_preinstall()

        # install dependencies
        # format: (str) a file name, say, "requirements.txt"
        pyitem = self.py.get("dependencies")
        if pyitem:
            self.logger.info("Install python dependencies")
            subprocess.check_call(["pip2", "install", "-r", os.path.join(self.source["python_root"], pyitem)])        

        # install python bin
        # format: [(<source fname>, <(optional) target_fname>), ...]
        pyitem = self.py.get("bin")
        if pyitem:
            self.logger.info("Install python bin")
            util.mkdir_p( self.install["bin"] )
            for each in pyitem:
                if 1 <= len(each) <= 2:
                    target_file = os.path.join(self.install["bin"], each[-1])
                    shutil.copyfile(os.path.join(self.source["python"], each[0]), target_file)
                    subprocess.check_call(["chmod", "+x", target_file])
                else:
                    raise ValueError("Unsupported tuple length: {}".format(str(each)))

        # install lib
        # format: [<lib_file1>, <lib_file2>, ...]
        # N.B. the "__init__.py" must be provided, since there might be a customized file of it
        pyitem = self.py.get("lib")
        if pyitem:
            self.logger.info("Install python lib")
            util.mkdir_p( self.install["libpython"] )
            for each in pyitem:
                shutil.copyfile(os.path.join(self.source["pylib"], each), os.path.join(self.install["libpython"], each))        

        # install libexec_bin
        # format: same as python bin; [(<source_fname>, <(optional) target_fname>), ...]
        pyitem = self.py.get("libexec_bin")
        if pyitem:
            self.logger.info("Install python libexec/bin")
            util.mkdir_p( self.install["libexec_bin"] )
            for each in pyitem:
                if 1 <= len(each) <= 2:
                    target_file = os.path.join(self.install["libexec_bin"], each[-1])
                    shutil.copyfile(os.path.join( self.source["python"], each[0]), target_file)
                    subprocess.check_call(["chmod", "+x", target_file])
                else:
                    raise ValueError("Unsupported tuple length: {}".format(str(each)))        

        self.install_python_postinstall()

    # Install CPP packages
    def install_cpp_preinstall(self):
        pass
    def install_cpp_postinstall(self):
        pass
    def install_cpp(self):
        if self.cpp == None: return

        self.install_cpp_preinstall()

        self.logger.info("Install cpp")
        util.mkdir_p( self.build_dir )
        # detect if homebrew is installed
        # format: True | False
        if self.cpp.get("detect_brew"):
            if distutils.spawn.find_executable("brew") != None:
                self.brew = subprocess.check_output(["brew", "--prefix"]).strip()

        cmd_args = ["cmake", self.source["cpp"], "-DCMAKE_INSTALL_PREFIX={}".format(self.install_dir)]
        if self.args.build_debug:
            cmd_args.append("-DCMAKE_BUILD_TYPE=Debug")
        else:
            cmd_args.append("-DCMAKE_BUILD_TYPE=Release")

        # cmake parameters
        cpp_parameters = self.cpp.get("parameters", [])
        for each in cpp_parameters:
            if len(each) == 1: # optional parameter
                value = getattr(self.args, each[0])
                if value != None:
                    cmd_args.append("-D{}".format(each[0]))
            elif len(each) == 2:
                value = getattr(self.args, each[0])
                if value != None:
                    cmd_args.append("-D{}={}".format(each[0], value))
                elif each[1] == False: # this one is optional
                    pass
                elif each[1] != None: # use the default value
                    cmd_args.append("-D{}={}".format(each[0], each[1]))
                else:
                    raise ValueError("CPP parameter is not set for {}".format(str(each)))
            elif len(each) == 3:
                value = getattr(self.args, each[0])
                if value != None:
                    cmd_args.append("-D{}={}".format(each[0], value))
                elif each[1] == "brew" and self.brew != None:
                    cmd_args.append("-D{}={}".format(each[0], self.brew))
                elif each[2] == False: # this one is optional
                    pass
                elif each[2] != None: # use the default value
                    cmd_args.append("-D{}={}".format(each[0], each[2]))
                else:
                    raise ValueError("CPP parameter is not set for {}".format(str(each)))
            else:
                raise ValueError("CPP parameter {} has illegal number of options".format(each))

        # add complicated cmake parameters
        cpp_commands = self.cpp.get("commands", [])
        for each in cpp_commands:
            each[0](self, cmd_args, *each[1:])

        os.chdir( self.build_dir )
        self.logger.info("run cmake: {}".format(str(cmd_args)))
        subprocess.check_call(cmd_args) # call cmake
        subprocess.check_call(["make", "install"])        

        self.install_cpp_postinstall()


    # Install makefiles
    def install_makefile_preinstall(self):
        pass
    def install_makefile_postinstall(self):
        pass

    def install_makefile(self):
        if self.makefile == None: return
        self.install_makefile_preinstall()

        self.logger.info("Install makefiles")
        util.mkdir_p( self.install["libexec_mk"] )
        for each in self.makefile:
            shutil.copyfile(os.path.join(self.source["makefile"], each), os.path.join(self.install["libexec_mk"], each))        

        self.install_makefile_postinstall()

    def prerun(self):
        pass

    def postrun(self):
        pass

    def run(self):
        self.prerun()

        util.mkdir_p( self.install_dir )
        self.logger.info("Install to the path: [{}]".format(self.install_dir))
        
        util.mkdir_p( self.install["bin"] )
        os.chdir( self.source_dir )

        self.install_version()
        self.install_license()
        self.install_python()
        self.install_cpp()
        self.install_makefile()

        self.postrun()

    

if __name__ == "__main__":
    def add_lemon_parameters(obj, cmd_args):
        obj.logger.info("add_lemon_parameters() called")

    MultiToolInstall(
            project = "oppac",
            version = "0.2",
            source_dir = os.path.dirname( __file__ ),
            start_hex = "ebe39546bc21fe9b299c857b4a69d81fe04952a0",
            python = {
                    "bin": [("oppac_alt.py", "oppac_alt"), ("oppac_dpc.py", "oppac_dpc")],
                    "lib": ["__init__.py", "util.py"],
                    "dependencies": "requirements.txt"
                },
            cpp = {
                    "detect_brew": True,
                    #"parameters": [("LOONLIB_ROOT_DIR", "brew", None)]
                    "parameters": [("LOGGER_LEVEL", "0")],

                    # self and cmd_args are the default parameters
                    "commands": [(add_lemon_parameters, )]
                }
        )    
