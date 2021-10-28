#! /usr/bin/env python2

import commands
commands.getstatusoutput("make clean")
commands.getstatusoutput("rm CMakeCache.txt cmake_install.cmake Makefile")
commands.getstatusoutput("rm -r CMakeFiles/")
