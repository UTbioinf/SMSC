#!/usr/bin/env python2

import os
import errno
import subprocess

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as err:
        if err.errno != errno.EEXIST:
            raise

def git_count(start_hex, end_hex = "HEAD"):
    return subprocess.check_output(["git", "rev-list", "--count", "{}..{}".format(start_hex, end_hex)]).rstrip()

