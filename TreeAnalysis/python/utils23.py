#!/usr/bin/env python

############################################
# Utilties which do not require ROOT       #
#                                          #
# Author: A. Mecca (alberto.mecca@cern.ch) #
############################################

import sys
if sys.version_info.major <= 2:
    from collections import Mapping
else:
    from collections.abc import Mapping


def deep_update(orig, new):
    for k, v in new.items():
        if isinstance(v, Mapping):
            orig[k] = deep_update(orig.get(k, {}), v)
        else:
            orig[k] = v
    return orig


def _test_deep_update():
    d1 = {'A': {'a1': 1}, 'B': {'b1': 1}}
    d2 = {'A': {'a2': 2}, 'B': {'b1': 3, 'b2': 2}}
    target = {'A': {'a1': 1, 'a2':2}, 'B': {'b1':3, 'b2': 2}}
    assert deep_update(d1, d2) == target, 'deep_update test failed'


def main():
    '''
    Run the tests
    '''
    _test_deep_update()

    return 0


if __name__ == '__main__':
    exit(main())
