#!/usr/bin/env python

import subprocess
subprocess.call(["make", "clean"])
subprocess.call(["rm", "CMakeCache.txt", "cmake_install.cmake", "Makefile"])
subprocess.call(["rm", "-r", "CMakeFiles/"])
