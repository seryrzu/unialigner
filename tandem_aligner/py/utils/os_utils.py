# (c) 2020 by Authors
# This file is a part of the SD program.
# see LICENSE file

import os
import errno


def smart_mkdir(dirname):
    try:
        os.mkdir(dirname)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc


def smart_makedirs(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc


def expandpath(path):
    expanduser_path = os.path.expanduser(path)
    real_path = os.path.realpath(expanduser_path)
    abs_path = os.path.abspath(real_path)
    return abs_path


def cat(infns, outfile):
    with open(outfile, "w") as o:
        for fn in infns:
            with open(fn) as i:
                for line in i:
                    o.write(line)
