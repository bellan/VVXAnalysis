#!/usr/bin/env python

import subprocess
subprocess.run(["make", "clean"])
subprocess.run(["rm", "CMakeCache.txt", "cmake_install.cmake", "Makefile"])
subprocess.run(["rm", "-r", "CMakeFiles/"])
