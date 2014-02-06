#! /usr/bin/env python

import ROOT, sys

fname = sys.argv[1]
print fname

fin  = ROOT.TFile(fname)
tree = fin.Get("treePlanter/ElderTree")
tree.Print("toponly")
