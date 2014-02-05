#! /usr/bin/env python

import commands
commands.getstatusoutput("make clean")
commands.getstatusoutput("rm CMakeCache.txt cmake_install.cmake Makefile")
commands.getstatusoutput("rm -r CMakeFiles/")
