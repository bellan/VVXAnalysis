#! /usr/bin/env python

import ROOT
import sys, os, commands, math, re, string

#inputdir = 'samples/CR2P2F_HZZ/'
inputdir = 'samples/'

dataset = ['DoubleMuA','DoubleMuB','DoubleMuC','DoubleMuD',
           'DoubleEleA','DoubleEleB','DoubleEleC','DoubleEleD',
           'MuEGA', 'MuEGB', 'MuEGC', 'MuEGD']

files = [inputdir+ds+'.root' for ds in dataset]


tree = ROOT.TChain('treePlanter/ElderTree')
for f in files:
    tree.Add(f)

entries = tree.GetEntries()

events = {}

for jentry in xrange(entries):
    # get the next tree in the chain and verify
    ientry = tree.LoadTree(jentry)
    if ientry < 0: break
     # copy next entry into memory and verify
    nb = tree.GetEntry(jentry)
    if nb<=0: continue
    # use the values directly from the tree
    # print tree.run, tree.event
    if tree.event in events: 
        print events[tree.event], tree.run, tree.event
    events[tree.event] = tree.run
