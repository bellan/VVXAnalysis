#! /usr/bin/env python

from __future__ import print_function
import ROOT, sys

fname = sys.argv[1]
print(fname)

fin  = ROOT.TFile(fname)
tree = fin.Get("treePlanter/ElderTree")
tree.Print("toponly")
