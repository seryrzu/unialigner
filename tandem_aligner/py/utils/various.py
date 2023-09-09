# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

import numpy as np


def list2str(lst, sep=" "):
    return sep.join(str(e) for e in lst)


def listEls2str(lst):
    return [str(e) for e in lst]


def fst_iterable(iterable):
    return next(iter(iterable))


def running_mean(data, window_size):
    cumsum = np.cumsum(np.insert(data, 0, 0))
    return (cumsum[window_size:] - cumsum[:-window_size]) / window_size
