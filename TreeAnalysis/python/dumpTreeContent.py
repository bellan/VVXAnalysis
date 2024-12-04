#!/usr/bin/env python2

import ROOT, sys

fname = sys.argv[1]
print fname

fin  = ROOT.TFile(fname)
tree = fin.Get("treePlanter/ElderTree")
tree.Print("toponly")
